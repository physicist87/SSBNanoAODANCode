#!/usr/bin/env python3

import argparse
import os
import re
import subprocess
import sys
from pathlib import Path
from collections import Counter


# ============================================================
# Paths
# ============================================================

CONDOR_DIR = Path(__file__).resolve().parent

LOG_BASE = CONDOR_DIR / "condorLog"

USER_ID = os.environ.get("USER") or os.getlogin()

SE_HOST = "root://cluster142.knu.ac.kr"
SE_BASE = f"/store/user/{USER_ID}/CPV_Run2/ULSummer20"


# ============================================================
# External InputList locations
#
# Keep this consistent with submit_jobs_v1.py
# ============================================================

INPUT_LIST_BASE = {
    "UL2018": Path(
        "/u/user/sha/Develop/CPviolation/SSB/AnalysisCode/"
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
# Sample utilities
# ============================================================

def is_data_sample(sample: str) -> bool:
    return sample.startswith("Data_")


def is_sample_for_channel(sample: str, channel: str, run_period: str) -> bool:

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
# Collect expected jobs from InputList
# ============================================================

def collect_input_lists(sample_dir: Path, sample: str):

    jobs = []

    pattern = re.compile(
        rf"^{re.escape(sample)}_(\d+)\.list$"
    )

    if not sample_dir.is_dir():
        return jobs

    for path in sample_dir.glob("*.list"):

        match = pattern.match(path.name)

        if not match:
            continue

        job_number = int(match.group(1))

        jobs.append(
            {
                "job_number": job_number,
                "list_path": path,
                "list_name": path.name,
                "root_name": f"{path.stem}.root",
            }
        )

    return sorted(
        jobs,
        key=lambda x: x["job_number"]
    )


# ============================================================
# Safe text reader
# ============================================================

def read_text(path: Path):

    if not path.is_file():
        return ""

    try:
        return path.read_text(
            errors="replace"
        )

    except Exception:
        return ""


# ============================================================
# Check Condor event log
#
# IMPORTANT:
#
# The same .log filename may contain records from multiple
# submissions/resubmissions.
#
# Therefore always use the LAST termination record.
# ============================================================

def check_condor_log(log_file: Path):

    if not log_file.is_file():

        return {
            "exists": False,
            "finished": False,
            "exit_code": None,
        }

    content = read_text(log_file)

    exit_codes = []

    # --------------------------------------------------------
    # Format 1:
    #
    # Normal termination (return value 0)
    # --------------------------------------------------------

    matches = re.findall(
        r"Normal termination\s+\(return value\s+(-?\d+)\)",
        content,
    )

    exit_codes.extend(
        int(code)
        for code in matches
    )

    # --------------------------------------------------------
    # Format 2:
    #
    # Job terminated ... with exit-code 0
    #
    # Use this only if Format 1 was not found.
    # --------------------------------------------------------

    if not exit_codes:

        matches = re.findall(
            r"exit-code\s+(-?\d+)",
            content,
        )

        exit_codes.extend(
            int(code)
            for code in matches
        )

    # --------------------------------------------------------
    # Use the LAST exit code.
    # --------------------------------------------------------

    exit_code = (
        exit_codes[-1]
        if exit_codes
        else None
    )

    finished = (
        exit_code is not None
    )

    return {
        "exists": True,
        "finished": finished,
        "exit_code": exit_code,
    }


# ============================================================
# Check worker stdout
# ============================================================

def check_stdout(out_file: Path):

    if not out_file.is_file():

        return {
            "exists": False,
            "analysis_done": False,
            "worker_done": False,
            "stageout_done": False,
        }

    content = read_text(out_file)

    return {
        "exists": True,

        "analysis_done":
            "Analysis destructor completed." in content,

        "worker_done":
            "=== [WORKER DONE]" in content,

        "stageout_done":
            "=== [STAGE-OUT CHECK]" in content,
    }


# ============================================================
# Check stderr
#
# Warnings are NOT treated as failures.
# ============================================================

def check_stderr(err_file: Path):

    if not err_file.is_file():

        return {
            "exists": False,
            "has_error": False,
            "serious_lines": [],
        }

    content = read_text(err_file)

    serious_patterns = [
        r"segmentation violation",
        r"segmentation fault",
        r"Traceback \(most recent call last\)",
        r"\[ERROR\]",
        r"fatal error",
        r"std::exception",
        r"terminate called",
        r"Aborted",
        r"core dumped",
    ]

    serious_lines = []

    for line in content.splitlines():

        for pattern in serious_patterns:

            if re.search(
                pattern,
                line,
                re.IGNORECASE,
            ):
                serious_lines.append(line)
                break

    return {
        "exists": True,
        "has_error": bool(serious_lines),
        "serious_lines": serious_lines,
    }


# ============================================================
# Get ROOT outputs on SE
#
# One xrdfs call per sample.
# ============================================================

def get_se_root_files(
    study: str,
    run_period: str,
    channel: str,
    sample: str,
):

    se_dir = (
        f"{SE_BASE}/"
        f"{study}/"
        f"{run_period}/"
        f"{channel}/"
        f"{sample}"
    )

    cmd = [
        "xrdfs",
        SE_HOST,
        "ls",
        se_dir,
    ]

    result = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )

    if result.returncode != 0:
        return set(), se_dir

    files = set()

    for line in result.stdout.splitlines():

        name = Path(
            line.strip()
        ).name

        if name.endswith(".root"):
            files.add(name)

    return files, se_dir


# ============================================================
# Determine job status
# ============================================================

def determine_status(
    condor_info,
    out_info,
    err_info,
    output_exists,
):

    # --------------------------------------------------------
    # Strong success condition
    #
    # 1. Analysis completed
    # 2. Worker wrapper reached the end
    # 3. Output exists on SE
    # 4. Latest Condor exit code is 0
    # --------------------------------------------------------

    if (
        out_info["analysis_done"]
        and out_info["worker_done"]
        and output_exists
        and condor_info["exit_code"] == 0
    ):
        return "OK"

    # --------------------------------------------------------
    # Explicit non-zero Condor exit
    # --------------------------------------------------------

    if (
        condor_info["exit_code"] is not None
        and condor_info["exit_code"] != 0
    ):
        return "FAILED"

    # --------------------------------------------------------
    # Serious stderr error
    # --------------------------------------------------------

    if err_info["has_error"]:
        return "FAILED"

    # --------------------------------------------------------
    # Analysis finished but wrapper did not finish.
    #
    # Usually stage-out / SE problem.
    # --------------------------------------------------------

    if (
        out_info["analysis_done"]
        and not out_info["worker_done"]
    ):
        return "STAGEOUT_FAIL"

    # --------------------------------------------------------
    # ROOT exists on SE, but local logs are incomplete.
    # --------------------------------------------------------

    if output_exists:
        return "OUTPUT_ONLY"

    # --------------------------------------------------------
    # Condor log exists but no termination record yet.
    # --------------------------------------------------------

    if (
        condor_info["exists"]
        and not condor_info["finished"]
    ):
        return "RUNNING"

    # --------------------------------------------------------
    # Partial logs exist.
    # --------------------------------------------------------

    if (
        condor_info["exists"]
        or out_info["exists"]
        or err_info["exists"]
    ):
        return "INCOMPLETE"

    return "MISSING"


# ============================================================
# Check one sample
# ============================================================

def check_sample(
    args,
    input_base: Path,
    sample: str,
):

    sample_dir = (
        input_base
        / sample
    )

    jobs = collect_input_lists(
        sample_dir,
        sample,
    )

    if not jobs:

        print(
            f"[WARNING] No input lists found for {sample}"
        )

        return None

    log_dir = (
        LOG_BASE
        / args.study
        / args.run_period
        / args.channel
        / sample
    )

    # --------------------------------------------------------
    # Fetch SE contents once per sample.
    # --------------------------------------------------------

    se_files, se_dir = get_se_root_files(
        args.study,
        args.run_period,
        args.channel,
        sample,
    )

    results = []

    for job in jobs:

        job_number = (
            job["job_number"]
        )

        prefix = (
            f"{sample}_{job_number}"
        )

        log_file = (
            log_dir
            / f"{prefix}.log"
        )

        out_file = (
            log_dir
            / f"{prefix}.out"
        )

        err_file = (
            log_dir
            / f"{prefix}.err"
        )

        condor_info = check_condor_log(
            log_file
        )

        out_info = check_stdout(
            out_file
        )

        err_info = check_stderr(
            err_file
        )

        output_exists = (
            job["root_name"]
            in se_files
        )

        status = determine_status(
            condor_info,
            out_info,
            err_info,
            output_exists,
        )

        results.append(
            {
                **job,

                "status":
                    status,

                "exit_code":
                    condor_info["exit_code"],

                "analysis_done":
                    out_info["analysis_done"],

                "worker_done":
                    out_info["worker_done"],

                "stageout_done":
                    out_info["stageout_done"],

                "output_exists":
                    output_exists,

                "errors":
                    err_info["serious_lines"],
            }
        )

    counts = Counter(
        item["status"]
        for item in results
    )

    ok = counts.get(
        "OK",
        0,
    )

    total = len(results)

    bad = (
        total - ok
    )

    print(
        f"{sample:<50} "
        f"{ok:>5}/{total:<5} "
        f"OK"
        + (
            f"   BAD={bad}"
            if bad
            else ""
        )
    )

    # --------------------------------------------------------
    # Show bad jobs automatically.
    #
    # --details shows successful jobs too.
    # --------------------------------------------------------

    if args.details or bad:

        for item in results:

            if (
                not args.details
                and item["status"] == "OK"
            ):
                continue

            exit_code = (
                item["exit_code"]
            )

            exit_string = (
                "-"
                if exit_code is None
                else str(exit_code)
            )

            print(
                f"    "
                f"{item['job_number']:>5} "
                f"{item['status']:<14} "
                f"exit={exit_string:<4} "
                f"analysis={str(item['analysis_done']):<5} "
                f"worker={str(item['worker_done']):<5} "
                f"SE={str(item['output_exists']):<5}"
            )

            if (
                args.show_errors
                and item["errors"]
            ):

                for line in item["errors"][:5]:

                    print(
                        f"          ERR: {line}"
                    )

    return {
        "sample":
            sample,

        "results":
            results,

        "counts":
            counts,

        "se_dir":
            se_dir,
    }


# ============================================================
# Main
# ============================================================

def run(args):

    # --------------------------------------------------------
    # InputList configuration
    # --------------------------------------------------------

    if (
        args.run_period
        not in INPUT_LIST_BASE
    ):

        print(
            "[ERROR] InputList base path "
            f"is not configured for {args.run_period}"
        )

        sys.exit(1)

    input_base = (
        INPUT_LIST_BASE[
            args.run_period
        ]
    )

    if not input_base.is_dir():

        print(
            "[ERROR] InputList directory not found:"
        )

        print(
            f"        {input_base}"
        )

        sys.exit(1)

    # --------------------------------------------------------
    # All available samples
    # --------------------------------------------------------

    all_samples = sorted(
        p.name
        for p in input_base.iterdir()
        if p.is_dir()
    )

    # --------------------------------------------------------
    # Sample selector
    # --------------------------------------------------------

    if "all" in args.samples:

        samples = (
            all_samples
        )

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

        samples = (
            args.samples
        )

    # --------------------------------------------------------
    # Channel-specific Data filtering
    # --------------------------------------------------------

    samples = [
        sample
        for sample in samples
        if is_sample_for_channel(
            sample,
            args.channel,
            args.run_period,
        )
    ]

    print("=" * 78)
    print("SSB CONDOR JOB CHECK")
    print("=" * 78)

    print(
        f"Study       : {args.study}"
    )

    print(
        f"Run period  : {args.run_period}"
    )

    print(
        f"Channel     : {args.channel}"
    )

    print(
        f"User        : {USER_ID}"
    )

    print(
        f"Input base  : {input_base}"
    )

    print(
        f"SE base     : {SE_HOST}/{SE_BASE}"
    )

    print(
        f"Samples     : {len(samples)}"
    )

    print("=" * 78)

    all_results = []

    # --------------------------------------------------------
    # Check samples
    # --------------------------------------------------------

    for sample in samples:

        if not (
            input_base
            / sample
        ).is_dir():

            print(
                f"[WARNING] Sample not found: {sample}"
            )

            continue

        result = check_sample(
            args,
            input_base,
            sample,
        )

        if result:
            all_results.append(
                result
            )

    # --------------------------------------------------------
    # Global summary
    # --------------------------------------------------------

    global_counts = Counter()

    total_jobs = 0

    bad_jobs = []

    for sample_result in all_results:

        for status, count in (
            sample_result["counts"].items()
        ):

            global_counts[
                status
            ] += count

        for item in (
            sample_result["results"]
        ):

            total_jobs += 1

            if item["status"] != "OK":

                bad_jobs.append(
                    (
                        sample_result["sample"],
                        item,
                    )
                )

    print()
    print("=" * 78)
    print("SUMMARY")
    print("=" * 78)

    print(
        f"Total expected jobs : {total_jobs}"
    )

    print(
        f"OK                  : {global_counts['OK']}"
    )

    for status in [
        "RUNNING",
        "FAILED",
        "STAGEOUT_FAIL",
        "INCOMPLETE",
        "OUTPUT_ONLY",
        "MISSING",
    ]:

        count = (
            global_counts.get(
                status,
                0,
            )
        )

        if count:

            print(
                f"{status:<20}: {count}"
            )

    print("-" * 78)

    print(
        f"Bad / unfinished    : "
        f"{len(bad_jobs)}"
    )

    print("=" * 78)

    # --------------------------------------------------------
    # Optional bad-job list
    # --------------------------------------------------------

    if args.write_bad:

        output_file = (
            CONDOR_DIR
            / (
                f"bad_jobs_"
                f"{args.study}_"
                f"{args.run_period}_"
                f"{args.channel}.txt"
            )
        )

        with output_file.open(
            "w"
        ) as f:

            for sample, item in bad_jobs:

                f.write(
                    f"{sample} "
                    f"{item['job_number']} "
                    f"{item['status']} "
                    f"{item['list_path']}\n"
                )

        print()
        print(
            "[INFO] Bad-job list written to:"
        )

        print(
            f"       {output_file}"
        )

    # --------------------------------------------------------
    # Exit status
    #
    # 0 = all jobs complete
    # 1 = one or more bad / unfinished jobs
    # --------------------------------------------------------

    if bad_jobs:
        sys.exit(1)

    sys.exit(0)


# ============================================================
# CLI
# ============================================================

def parse_args():

    parser = argparse.ArgumentParser(
        description=(
            "Check SSBNanoAOD HTCondor jobs "
            "against InputList, Condor logs, "
            "and SE ROOT outputs."
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

    parser.add_argument(
        "--samples",
        nargs="+",
        default=["all"],
        help=(
            "all, data, mc, or explicit sample names. "
            "Default: all"
        ),
    )

    parser.add_argument(
        "--details",
        action="store_true",
        help=(
            "Show status of every job, "
            "including successful jobs."
        ),
    )

    parser.add_argument(
        "--show-errors",
        action="store_true",
        help=(
            "Print serious stderr lines "
            "for failed jobs."
        ),
    )

    parser.add_argument(
        "--write-bad",
        action="store_true",
        help=(
            "Write bad/unfinished jobs "
            "to a text file."
        ),
    )

    return parser.parse_args()


# ============================================================
# Entry point
# ============================================================

if __name__ == "__main__":

    args = parse_args()

    run(args)
