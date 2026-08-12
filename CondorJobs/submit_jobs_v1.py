#!/usr/bin/env python3

import argparse
import os
import re
import subprocess
import sys
from pathlib import Path


# ============================================================
# Global paths
# ============================================================

CONDOR_DIR = Path(__file__).resolve().parent
PACKAGE_DIR = CONDOR_DIR.parent

PACKAGE_NAME = "SSBNanoAODANCode"
TARBALL_PATH = CONDOR_DIR / f"{PACKAGE_NAME}.tar.gz"

RUN_SCRIPT = CONDOR_DIR / "run_condor_v1.sh"

LOG_BASE = CONDOR_DIR / "condorLog"
SUBMIT_BASE = CONDOR_DIR / "condorSubmit"


# ============================================================
# External InputList locations
# ============================================================

INPUT_LIST_BASE = {
    "UL2018": Path(
        "/u/user/sha/Develop/CPviolation/SSB/AnalysisCode/"
        "NanoAODNtuple_v1/NtupleList_v1/"
        "2018_v6-FromGuks/FileList/InputList"
    ),

    # Later:
    # "UL2017": Path("..."),
    # "UL2016PreVFP": Path("..."),
    # "UL2016PostVFP": Path("..."),
}


# ============================================================
# Data filtering
# ============================================================

def is_data_sample(sample: str) -> bool:
    return sample.startswith("Data_")


def is_sample_for_channel(sample: str, channel: str, run_period: str) -> bool:
    """
    Keep only data datasets relevant for the selected channel.
    MC samples always pass.
    """

    if not is_data_sample(sample):
        return True

    dimuon_list = ["SingleMuon", "DoubleMuon"]
    dielec_list = ["SingleElectron", "DoubleEG"]
    muelec_list = ["SingleMuon", "SingleElectron", "MuonEG"]

    if run_period == "UL2018":
        # Keep compatibility with possible EGamma-style samples.
        dielec_list = ["EGamma", "SingleElectron", "DoubleEG"]
        muelec_list = ["SingleMuon", "EGamma", "SingleElectron", "MuonEG"]

    if channel == "MuMu":
        allowed = dimuon_list
    elif channel == "ElEl":
        allowed = dielec_list
    elif channel == "MuEl":
        allowed = muelec_list
    else:
        return False

    return any(token in sample for token in allowed)


# ============================================================
# Config / branch-list mapping
# ============================================================

def get_job_config(run_period: str, channel: str, sample: str):
    """
    Determine config and branch-list files.
    Start from the mapping used in the previous unified submitter.
    """

    config_map = {
        "MuMu": "dimuon.config",
        "ElEl": "dielec.config",
        "MuEl": "muelec.config",
    }

    config_file = config_map[channel]

    # Current NanoAOD v15 branch list
    branch_list = f"{run_period}/branch_list_v15.txt"

    # Preserve room for special year/data mappings later.
    #
    # Example:
    #
    # if run_period == "UL2017" and sample.startswith("Data_"):
    #     ...
    #

    return config_file, branch_list


# ============================================================
# Create analysis tarball
# ============================================================

def create_package_tarball():
    """
    Create SSBNanoAODANCode.tar.gz immediately before submission.

    Excluded:
      - input
      - output
      - CondorJobs
      - .git
      - local build products
      - common temporary files

    Everything else is included by default so that runtime
    correction/config/data dependencies are not accidentally omitted.
    """

    print("=" * 70)
    print("[PACKAGE] Creating analysis tarball")
    print("=" * 70)

    if TARBALL_PATH.exists():
        print(f"[INFO] Removing old tarball: {TARBALL_PATH}")
        TARBALL_PATH.unlink()

    cmd = [
        "tar",
        "-czf",
        str(TARBALL_PATH),

        f"--exclude={PACKAGE_NAME}/input",
        f"--exclude={PACKAGE_NAME}/output",
        f"--exclude={PACKAGE_NAME}/CondorJobs",
        f"--exclude={PACKAGE_NAME}/.git",

        f"--exclude={PACKAGE_NAME}/*.o",
        f"--exclude={PACKAGE_NAME}/ssb_analysis",

        f"--exclude={PACKAGE_NAME}/.DS_Store",
        f"--exclude={PACKAGE_NAME}/__pycache__",
        f"--exclude={PACKAGE_NAME}/*.pyc",

        PACKAGE_NAME,
    ]

    subprocess.run(
        cmd,
        cwd=PACKAGE_DIR.parent,
        check=True,
    )

    size_mb = TARBALL_PATH.stat().st_size / (1024 * 1024)

    print(f"[INFO] Tarball : {TARBALL_PATH}")
    print(f"[INFO] Size    : {size_mb:.1f} MB")
    print("=" * 70)


# ============================================================
# Collect list files
# ============================================================

def collect_input_lists(sample_dir: Path, sample: str):
    """
    Collect files matching:

        SAMPLE/SAMPLE_<number>.list

    Returns:
        [(number, full_path), ...]
    """

    jobs = []

    if not sample_dir.is_dir():
        return jobs

    pattern = re.compile(
        rf"^{re.escape(sample)}_(\d+)\.list$"
    )

    for path in sorted(sample_dir.glob("*.list")):
        match = pattern.match(path.name)

        if not match:
            continue

        job_number = int(match.group(1))
        jobs.append((job_number, path))

    return jobs


# ============================================================
# Condor submission
# ============================================================

def submit_condor(jdl_path: Path):
    cmd = ["condor_submit", str(jdl_path)]

    print(f"[INFO] Running: {' '.join(cmd)}")

    subprocess.run(
        cmd,
        check=True,
    )


# ============================================================
# Generate one JDL per sample
# ============================================================

def generate_sample_jdl(
    args,
    sample,
    jobs,
    config_file,
    branch_list,
):
    submit_dir = (
        SUBMIT_BASE
        / args.study
        / args.run_period
        / args.channel
    )

    log_dir = (
        LOG_BASE
        / args.study
        / args.run_period
        / args.channel
        / sample
    )

    submit_dir.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)

    jdl_path = submit_dir / f"{sample}.sub"

    #
    # Queue file contains:
    #
    # JobId InputListPath InputListName
    #
    queue_file = submit_dir / f"{sample}.queue"

    with queue_file.open("w") as f:
        for job_number, list_path in jobs:
            f.write(
                f"{job_number} "
                f"{list_path} "
                f"{list_path.name}\n"
            )

    with jdl_path.open("w") as f:

        f.write("Universe = vanilla\n")
        f.write(f"Executable = {RUN_SCRIPT}\n")

        f.write(
            f"Log = {log_dir}/{sample}_$(JobId).log\n"
        )
        f.write(
            f"Output = {log_dir}/{sample}_$(JobId).out\n"
        )
        f.write(
            f"Error = {log_dir}/{sample}_$(JobId).err\n"
        )

        f.write("\n")

        # KNU resources: start conservatively.
        f.write("RequestCpus = 1\n")
        f.write("RequestMemory = 4 GB\n")
        f.write("RequestDisk = 10 GB\n")

        f.write('+JobType = "short"\n')

        f.write("\n")

        f.write("should_transfer_files = YES\n")
        f.write("when_to_transfer_output = ON_EXIT\n")
        f.write("use_x509userproxy = true\n")

        #
        # Important:
        #
        # Do NOT transfer the output ROOT file back to the UI.
        # run_condor_v1.sh performs direct xrdcp stage-out.
        #
        f.write("transfer_output_files = \"\"\n")

        f.write("\n")

        #
        # Per-job transfer:
        #
        #   SSBNanoAODANCode.tar.gz
        #   one .list file
        #
        f.write(
            "transfer_input_files = "
            f"{TARBALL_PATH},$(InputListPath)\n"
        )

        f.write("\n")

        #
        # run_condor_v1.sh arguments
        #
        f.write(
            "Arguments = "
            f"\"{args.run_period} "
            f"{args.study} "
            f"{args.channel} "
            f"{sample} "
            f"{config_file} "
            f"{branch_list} "
            f"{args.max_events} "
            f"$(InputListName)\"\n"
        )

        f.write("\n")

        f.write(
            "Queue JobId, InputListPath, InputListName "
            f"from {queue_file}\n"
        )

    return jdl_path


# ============================================================
# Main workflow
# ============================================================

def run(args):

    if args.run_period not in INPUT_LIST_BASE:
        print(
            f"[ERROR] InputList base path is not configured "
            f"for {args.run_period}"
        )
        sys.exit(1)

    input_base = INPUT_LIST_BASE[args.run_period]

    if not input_base.is_dir():
        print(f"[ERROR] InputList directory not found:")
        print(f"        {input_base}")
        sys.exit(1)

    if not RUN_SCRIPT.is_file():
        print(f"[ERROR] run_condor_v1.sh not found:")
        print(f"        {RUN_SCRIPT}")
        sys.exit(1)

    #
    # Create package once per invocation.
    #
    create_package_tarball()

    #
    # Determine requested samples
    #
    all_samples = sorted(
        p.name
        for p in input_base.iterdir()
        if p.is_dir()
    )

    if "all" in args.samples:
        samples = all_samples
    else:
        samples = args.samples

    print()
    print("=" * 70)
    print("[SUBMISSION CONFIGURATION]")
    print("=" * 70)
    print(f"Study       : {args.study}")
    print(f"Run period  : {args.run_period}")
    print(f"Channel     : {args.channel}")
    print(f"Max events  : {args.max_events}")
    print(f"Input base  : {input_base}")
    print(f"Samples     : {', '.join(samples)}")
    print("=" * 70)

    total_jobs = 0
    jdl_files = []

    for sample in samples:

        sample_dir = input_base / sample

        if not sample_dir.is_dir():
            print(
                f"[WARNING] Sample directory does not exist: "
                f"{sample_dir}"
            )
            continue

        if not is_sample_for_channel(
            sample,
            args.channel,
            args.run_period,
        ):
            print(
                f"[INFO] Skip data sample not used by "
                f"{args.channel}: {sample}"
            )
            continue

        config_file, branch_list = get_job_config(
            args.run_period,
            args.channel,
            sample,
        )

        jobs = collect_input_lists(
            sample_dir,
            sample,
        )

        if not jobs:
            print(
                f"[WARNING] No valid .list files found: {sample}"
            )
            continue

        #
        # Useful for initial Condor validation:
        #
        if args.test:
            jobs = jobs[:1]

        jdl_path = generate_sample_jdl(
            args,
            sample,
            jobs,
            config_file,
            branch_list,
        )

        total_jobs += len(jobs)
        jdl_files.append(jdl_path)

        print(
            f"[READY] {sample:<45} "
            f"{len(jobs):>5} jobs"
        )

    print()
    print("=" * 70)
    print(f"[SUMMARY] Total jobs : {total_jobs}")
    print(f"[SUMMARY] JDL files  : {len(jdl_files)}")
    print("=" * 70)

    if args.dry_run:
        print("[INFO] --dry-run: no jobs submitted.")
        return

    for jdl in jdl_files:
        submit_condor(jdl)


# ============================================================
# CLI
# ============================================================

def parse_args():

    parser = argparse.ArgumentParser(
        description="SSBNanoAOD analysis HTCondor submitter"
    )

    parser.add_argument(
        "--study",
        required=True,
        help="Study name, e.g. Testv1, AN_v7, AN_v7_JESUp",
    )

    parser.add_argument(
        "--run-period",
        required=True,
        choices=[
            "UL2016PreVFP",
            "UL2016PostVFP",
            "UL2017",
            "UL2018",
        ],
    )

    parser.add_argument(
        "--channel",
        required=True,
        choices=[
            "MuMu",
            "ElEl",
            "MuEl",
        ],
    )

    parser.add_argument(
        "--samples",
        nargs="+",
        required=True,
        help=(
            "Samples to submit. "
            "Use '--samples all' for all available samples."
        ),
    )

    parser.add_argument(
        "--max-events",
        default="-1",
        help="Maximum events per job. Default: -1",
    )

    parser.add_argument(
        "--test",
        action="store_true",
        help="Submit only the first list file of each selected sample.",
    )

    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Generate tarball/JDLs but do not run condor_submit.",
    )

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    run(args)
