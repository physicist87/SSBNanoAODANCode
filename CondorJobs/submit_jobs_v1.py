#!/usr/bin/env python3

import argparse
import os
import re
import subprocess
import sys
from pathlib import Path
from collections import defaultdict


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
# User
# ============================================================

USER_ID = os.environ.get("USER") or os.getlogin()


# ============================================================
# External InputList locations
# ============================================================

INPUT_LIST_BASE = {

    "UL2018": Path(
        f"/u/user/{USER_ID}/Develop/CPviolation/SSB/AnalysisCode/"
        "NanoAODNtuple_v1/NtupleList_v1/"
        "2018_v6-FromGuks/FileList/InputList"
    ),

    # Add later:
    #
    # "UL2017": Path("..."),
    # "UL2016PreVFP": Path("..."),
    # "UL2016PostVFP": Path("..."),
}


# ============================================================
# Data utilities
# ============================================================

def is_data_sample(sample: str) -> bool:
    return sample.startswith("Data_")


def is_sample_for_channel(
    sample: str,
    channel: str,
    run_period: str,
) -> bool:

    if not is_data_sample(sample):
        return True

    dimuon_list = [
        "SingleMuon",
        "DoubleMuon",
    ]

    dielec_list = [
        "SingleElectron",
        "DoubleEG",
    ]

    muelec_list = [
        "SingleMuon",
        "SingleElectron",
        "MuonEG",
    ]

    if run_period == "UL2018":

        dielec_list = [
            "EGamma",
            "SingleElectron",
            "DoubleEG",
        ]

        muelec_list = [
            "SingleMuon",
            "EGamma",
            "SingleElectron",
            "MuonEG",
        ]

    if channel == "MuMu":
        allowed = dimuon_list

    elif channel == "ElEl":
        allowed = dielec_list

    elif channel == "MuEl":
        allowed = muelec_list

    else:
        return False

    return any(
        token in sample
        for token in allowed
    )


# ============================================================
# Config / branch-list mapping
# ============================================================

def get_job_config(
    run_period: str,
    channel: str,
    sample: str,
):

    config_map = {
        "MuMu": "dimuon.config",
        "ElEl": "dielec.config",
        "MuEl": "muelec.config",
    }

    config_file = config_map[channel]

    branch_list = (
        f"{run_period}/branch_list_v15.txt"
    )

    return config_file, branch_list


# ============================================================
# Create analysis tarball
# ============================================================

def create_package_tarball():

    print("=" * 70)
    print("[PACKAGE] Creating analysis tarball")
    print("=" * 70)

    if TARBALL_PATH.exists():

        print(
            f"[INFO] Removing old tarball: "
            f"{TARBALL_PATH}"
        )

        TARBALL_PATH.unlink()

    cmd = [
        "tar",
        "-czf",
        str(TARBALL_PATH),

        # Local/runtime directories not needed on worker
        f"--exclude={PACKAGE_NAME}/input",
        f"--exclude={PACKAGE_NAME}/output",
        f"--exclude={PACKAGE_NAME}/CondorJobs",

        # Git
        f"--exclude={PACKAGE_NAME}/.git",

        # Local build products
        f"--exclude={PACKAGE_NAME}/*.o",
        f"--exclude={PACKAGE_NAME}/ssb_analysis",

        # Temporary files
        f"--exclude={PACKAGE_NAME}/.DS_Store",
        f"--exclude={PACKAGE_NAME}/__pycache__",
        f"--exclude={PACKAGE_NAME}/*.pyc",

        PACKAGE_NAME,
    ]

    print(
        f"[INFO] Package directory : "
        f"{PACKAGE_DIR}"
    )

    print(
        f"[INFO] Output tarball    : "
        f"{TARBALL_PATH}"
    )

    subprocess.run(
        cmd,
        cwd=PACKAGE_DIR.parent,
        check=True,
    )

    size_mb = (
        TARBALL_PATH.stat().st_size
        / (1024 * 1024)
    )

    print(
        f"[INFO] Tarball created successfully: "
        f"{size_mb:.1f} MB"
    )

    print("=" * 70)


# ============================================================
# Collect normal jobs from InputList
# ============================================================

def collect_input_lists(
    sample_dir: Path,
    sample: str,
):

    jobs = []

    if not sample_dir.is_dir():
        return jobs

    pattern = re.compile(
        rf"^{re.escape(sample)}_(\d+)\.list$"
    )

    for path in sorted(
        sample_dir.glob("*.list")
    ):

        match = pattern.match(
            path.name
        )

        if not match:
            continue

        job_number = int(
            match.group(1)
        )

        jobs.append(
            (
                job_number,
                path.resolve(),
            )
        )

    return sorted(
        jobs,
        key=lambda x: x[0],
    )


# ============================================================
# Read bad-job file
#
# Expected format from check_jobs_v1.py:
#
# SAMPLE JOB_NUMBER STATUS INPUT_LIST_PATH
#
# Example:
#
# TTZToLLNuNu 27 RUNNING /.../TTZToLLNuNu_27.list
# TTbar_AllHadronic 141 RUNNING /.../TTbar_AllHadronic_141.list
# ============================================================

def read_bad_jobs(
    bad_jobs_file: Path,
):

    if not bad_jobs_file.is_file():

        print(
            "[ERROR] Bad-job file not found:"
        )

        print(
            f"        {bad_jobs_file}"
        )

        sys.exit(1)

    jobs_by_sample = defaultdict(list)

    print("=" * 70)
    print("[BAD-JOB MODE]")
    print("=" * 70)

    print(
        f"[INFO] Reading bad-job file:"
    )

    print(
        f"       {bad_jobs_file}"
    )

    with bad_jobs_file.open() as f:

        for line_number, line in enumerate(
            f,
            start=1,
        ):

            line = line.strip()

            if not line:
                continue

            if line.startswith("#"):
                continue

            parts = line.split(
                maxsplit=3
            )

            if len(parts) != 4:

                print(
                    f"[WARNING] Invalid line "
                    f"{line_number}:"
                )

                print(
                    f"          {line}"
                )

                continue

            sample = parts[0]
            job_number_text = parts[1]
            status = parts[2]
            input_list_text = parts[3]

            try:

                job_number = int(
                    job_number_text
                )

            except ValueError:

                print(
                    f"[WARNING] Invalid job number "
                    f"at line {line_number}: "
                    f"{job_number_text}"
                )

                continue

            input_list = Path(
                input_list_text
            )

            if not input_list.is_file():

                print(
                    f"[ERROR] Input list from "
                    f"bad-job file does not exist:"
                )

                print(
                    f"        {input_list}"
                )

                sys.exit(1)

            jobs_by_sample[
                sample
            ].append(
                (
                    job_number,
                    input_list.resolve(),
                )
            )

            print(
                f"[BAD] {sample:<45} "
                f"job={job_number:<5} "
                f"previous_status={status}"
            )

    # Remove duplicate job numbers
    # and sort each sample.

    cleaned = {}

    for sample, jobs in jobs_by_sample.items():

        unique_jobs = {}

        for job_number, input_list in jobs:

            unique_jobs[
                job_number
            ] = input_list

        cleaned[sample] = sorted(
            unique_jobs.items(),
            key=lambda x: x[0],
        )

    total = sum(
        len(jobs)
        for jobs in cleaned.values()
    )

    print("-" * 70)

    print(
        f"[INFO] Samples to resubmit : "
        f"{len(cleaned)}"
    )

    print(
        f"[INFO] Jobs to resubmit    : "
        f"{total}"
    )

    print("=" * 70)

    if total == 0:

        print(
            "[INFO] No bad jobs found. "
            "Nothing to submit."
        )

        sys.exit(0)

    return cleaned


# ============================================================
# Condor submission
# ============================================================

def submit_condor(
    jdl_path: Path,
):

    cmd = [
        "condor_submit",
        str(jdl_path),
    ]

    print(
        f"[INFO] Running: "
        f"{' '.join(cmd)}"
    )

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
    is_resubmit=False,
):

    # --------------------------------------------------------
    # Directory layout
    # --------------------------------------------------------

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

    submit_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    log_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    # --------------------------------------------------------
    # Filename
    # --------------------------------------------------------

    if is_resubmit:

        jdl_path = (
            submit_dir
            / f"{sample}_resubmit.sub"
        )

        queue_file = (
            submit_dir
            / f"{sample}_resubmit.queue"
        )

    else:

        jdl_path = (
            submit_dir
            / f"{sample}.sub"
        )

        queue_file = (
            submit_dir
            / f"{sample}.queue"
        )

    # --------------------------------------------------------
    # Queue file
    #
    # JobId InputListPath InputListName
    # --------------------------------------------------------

    with queue_file.open("w") as f:

        for job_number, list_path in jobs:

            f.write(
                f"{job_number} "
                f"{list_path} "
                f"{list_path.name}\n"
            )

    # --------------------------------------------------------
    # JDL
    # --------------------------------------------------------

    with jdl_path.open("w") as f:

        f.write(
            "Universe = vanilla\n"
        )

        f.write(
            f"Executable = {RUN_SCRIPT}\n"
        )

        f.write("\n")

        # ----------------------------------------------------
        # Use the same log filenames as original jobs.
        #
        # check_jobs_v1.py uses the LAST Condor termination
        # record, so resubmissions are handled correctly.
        # ----------------------------------------------------

        f.write(
            f"Log = "
            f"{log_dir}/{sample}_$(JobId).log\n"
        )

        f.write(
            f"Output = "
            f"{log_dir}/{sample}_$(JobId).out\n"
        )

        f.write(
            f"Error = "
            f"{log_dir}/{sample}_$(JobId).err\n"
        )

        f.write("\n")

        # ----------------------------------------------------
        # Resources
        # ----------------------------------------------------

        f.write(
            "RequestCpus = 1\n"
        )

        f.write(
            "RequestMemory = 4 GB\n"
        )

        f.write(
            "RequestDisk = 10 GB\n"
        )

        f.write(
            '+JobType = "short"\n'
        )

        f.write("\n")

        # ----------------------------------------------------
        # Condor file transfer
        # ----------------------------------------------------

        f.write(
            "should_transfer_files = YES\n"
        )

        f.write(
            "when_to_transfer_output = ON_EXIT\n"
        )

        f.write(
            "use_x509userproxy = true\n"
        )

        # Output ROOT is sent directly from worker to SE
        # using xrdcp in run_condor_v1.sh.

        f.write(
            'transfer_output_files = ""\n'
        )

        f.write("\n")

        # ----------------------------------------------------
        # Input transfer
        #
        # 1. Analysis package tarball
        # 2. One input .list file
        # ----------------------------------------------------

        f.write(
            "transfer_input_files = "
            f"{TARBALL_PATH},$(InputListPath)\n"
        )

        f.write("\n")

        # ----------------------------------------------------
        # Worker arguments
        # ----------------------------------------------------

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

        # ----------------------------------------------------
        # Queue
        # ----------------------------------------------------

        f.write(
            "Queue "
            "JobId, InputListPath, InputListName "
            f"from {queue_file}\n"
        )

    return jdl_path, queue_file


# ============================================================
# Normal sample selection
# ============================================================

def select_samples(
    args,
    input_base: Path,
):

    all_samples = sorted(
        p.name
        for p in input_base.iterdir()
        if p.is_dir()
    )

    if "all" in args.samples:

        samples = all_samples

    elif "data" in args.samples:

        samples = [
            sample
            for sample in all_samples
            if is_data_sample(sample)
        ]

    elif "mc" in args.samples:

        samples = [
            sample
            for sample in all_samples
            if not is_data_sample(sample)
        ]

    else:

        samples = args.samples

    samples = [
        sample
        for sample in samples
        if is_sample_for_channel(
            sample,
            args.channel,
            args.run_period,
        )
    ]

    return samples


# ============================================================
# Normal submission mode
# ============================================================

def prepare_normal_jobs(
    args,
    input_base: Path,
):

    samples = select_samples(
        args,
        input_base,
    )

    jobs_by_sample = {}

    for sample in samples:

        sample_dir = (
            input_base
            / sample
        )

        if not sample_dir.is_dir():

            print(
                f"[WARNING] Sample directory "
                f"does not exist:"
            )

            print(
                f"          {sample_dir}"
            )

            continue

        jobs = collect_input_lists(
            sample_dir,
            sample,
        )

        if not jobs:

            print(
                f"[WARNING] No valid .list files "
                f"found: {sample}"
            )

            continue

        if args.test:
            jobs = jobs[:1]

        jobs_by_sample[
            sample
        ] = jobs

    return jobs_by_sample


# ============================================================
# Main workflow
# ============================================================

def run(args):

    # --------------------------------------------------------
    # Validate run period
    # --------------------------------------------------------

    if (
        args.run_period
        not in INPUT_LIST_BASE
    ):

        print(
            "[ERROR] InputList base path "
            f"is not configured for "
            f"{args.run_period}"
        )

        sys.exit(1)

    input_base = (
        INPUT_LIST_BASE[
            args.run_period
        ]
    )

    if not input_base.is_dir():

        print(
            "[ERROR] InputList directory "
            "not found:"
        )

        print(
            f"        {input_base}"
        )

        sys.exit(1)

    # --------------------------------------------------------
    # Worker script
    # --------------------------------------------------------

    if not RUN_SCRIPT.is_file():

        print(
            "[ERROR] run_condor_v1.sh "
            "not found:"
        )

        print(
            f"        {RUN_SCRIPT}"
        )

        sys.exit(1)

    # --------------------------------------------------------
    # Determine submission mode
    # --------------------------------------------------------

    is_resubmit = (
        args.bad_jobs is not None
    )

    # --------------------------------------------------------
    # Build list of jobs
    # --------------------------------------------------------

    if is_resubmit:

        bad_jobs_file = Path(
            args.bad_jobs
        )

        if not bad_jobs_file.is_absolute():

            bad_jobs_file = (
                Path.cwd()
                / bad_jobs_file
            )

        jobs_by_sample = read_bad_jobs(
            bad_jobs_file.resolve()
        )

    else:

        jobs_by_sample = prepare_normal_jobs(
            args,
            input_base,
        )

    if not jobs_by_sample:

        print(
            "[INFO] No jobs selected."
        )

        return

    # --------------------------------------------------------
    # Create analysis tarball once
    # --------------------------------------------------------

    create_package_tarball()

    # --------------------------------------------------------
    # Submission information
    # --------------------------------------------------------

    print()
    print("=" * 70)
    print("[SUBMISSION CONFIGURATION]")
    print("=" * 70)

    print(
        f"Mode         : "
        f"{'BAD-JOB RESUBMIT' if is_resubmit else 'NORMAL'}"
    )

    print(
        f"Study        : {args.study}"
    )

    print(
        f"Run period   : {args.run_period}"
    )

    print(
        f"Channel      : {args.channel}"
    )

    print(
        f"Max events   : {args.max_events}"
    )

    print(
        f"Input base   : {input_base}"
    )

    print(
        f"Samples      : "
        f"{len(jobs_by_sample)}"
    )

    print("=" * 70)

    total_jobs = 0
    jdl_files = []

    # --------------------------------------------------------
    # Generate JDLs
    # --------------------------------------------------------

    for sample in sorted(
        jobs_by_sample.keys()
    ):

        jobs = (
            jobs_by_sample[
                sample
            ]
        )

        if not is_sample_for_channel(
            sample,
            args.channel,
            args.run_period,
        ):

            print(
                f"[WARNING] Skip sample not "
                f"valid for {args.channel}: "
                f"{sample}"
            )

            continue

        config_file, branch_list = (
            get_job_config(
                args.run_period,
                args.channel,
                sample,
            )
        )

        jdl_path, queue_file = (
            generate_sample_jdl(
                args,
                sample,
                jobs,
                config_file,
                branch_list,
                is_resubmit=is_resubmit,
            )
        )

        total_jobs += len(
            jobs
        )

        jdl_files.append(
            jdl_path
        )

        if is_resubmit:

            print(
                f"[RESUBMIT] "
                f"{sample:<45} "
                f"{len(jobs):>5} jobs"
            )

            print(
                "           Job numbers: "
                + ", ".join(
                    str(job_number)
                    for job_number, _
                    in jobs
                )
            )

        else:

            print(
                f"[READY] "
                f"{sample:<45} "
                f"{len(jobs):>5} jobs"
            )

    # --------------------------------------------------------
    # Summary
    # --------------------------------------------------------

    print()
    print("=" * 70)

    if is_resubmit:

        print(
            f"[SUMMARY] Jobs to resubmit : "
            f"{total_jobs}"
        )

    else:

        print(
            f"[SUMMARY] Total jobs       : "
            f"{total_jobs}"
        )

    print(
        f"[SUMMARY] JDL files        : "
        f"{len(jdl_files)}"
    )

    print("=" * 70)

    # --------------------------------------------------------
    # Dry run
    # --------------------------------------------------------

    if args.dry_run:

        print(
            "[INFO] --dry-run: "
            "no jobs submitted."
        )

        return

    # --------------------------------------------------------
    # Submit
    # --------------------------------------------------------

    for jdl in jdl_files:

        submit_condor(
            jdl
        )


# ============================================================
# CLI
# ============================================================

def parse_args():

    parser = argparse.ArgumentParser(
        description=(
            "SSBNanoAOD analysis "
            "HTCondor submitter"
        )
    )

    parser.add_argument(
        "--study",
        required=True,
        help=(
            "Study name, e.g. "
            "NanoAODv15_Testv1"
        ),
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

    # --------------------------------------------------------
    # Selection mode
    #
    # Exactly one of --samples / --bad-jobs is required.
    # --------------------------------------------------------

    selection_group = (
        parser.add_mutually_exclusive_group(
            required=True
        )
    )

    selection_group.add_argument(
        "--samples",
        nargs="+",
        help=(
            "Samples to submit. "
            "Supported shortcuts: "
            "all, data, mc."
        ),
    )

    selection_group.add_argument(
        "--bad-jobs",
        help=(
            "Bad-job list produced by "
            "check_jobs_v1.py --write-bad. "
            "Only jobs listed in this file "
            "will be resubmitted."
        ),
    )

    parser.add_argument(
        "--max-events",
        default="-1",
        help=(
            "Maximum events per job. "
            "Default: -1"
        ),
    )

    parser.add_argument(
        "--test",
        action="store_true",
        help=(
            "Normal submission mode only: "
            "use the first list file of "
            "each selected sample."
        ),
    )

    parser.add_argument(
        "--dry-run",
        action="store_true",
        help=(
            "Generate tarball/JDL/queue files "
            "but do not run condor_submit."
        ),
    )

    args = parser.parse_args()

    # --test does not make sense in bad-job mode.
    if (
        args.bad_jobs
        and args.test
    ):

        parser.error(
            "--test cannot be used together "
            "with --bad-jobs"
        )

    return args


# ============================================================
# Entry point
# ============================================================

if __name__ == "__main__":

    args = parse_args()

    run(args)
