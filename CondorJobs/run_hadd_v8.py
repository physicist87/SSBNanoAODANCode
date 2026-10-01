#!/usr/bin/env python3

import argparse
import os
import re
import subprocess
import sys
from pathlib import Path


# ============================================================
# User / Paths
# ============================================================

USER_ID = os.environ.get("USER") or os.getlogin()

# Output SE, NFS view
OUTPUT_BASE = Path(
    f"/pnfs/knu.ac.kr/data/cms/store/user/{USER_ID}/"
    "CPV_Run2/ULSummer20"
)


# ============================================================
# Sample utilities
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
# Collect expected chunk outputs from InputList
# ============================================================

def collect_expected_outputs(
    input_sample_dir: Path,
    output_sample_dir: Path,
    sample: str,
):

    pattern = re.compile(
        rf"^{re.escape(sample)}_(\d+)\.list$"
    )

    expected = []

    if not input_sample_dir.is_dir():
        return expected

    for list_file in input_sample_dir.glob("*.list"):

        match = pattern.match(
            list_file.name
        )

        if not match:
            continue

        job_number = int(
            match.group(1)
        )

        root_file = (
            output_sample_dir
            / f"{sample}_{job_number}.root"
        )

        expected.append(
            {
                "job_number": job_number,
                "list_file": list_file,
                "root_file": root_file,
            }
        )

    return sorted(
        expected,
        key=lambda x: x["job_number"],
    )


# ============================================================
# Check whether all ROOT files in an InputList contain 0 events
#
# NtupleForge may legitimately produce a zero-event ROOT file
# when all events are rejected by Golden JSON filtering.
#
# This function is called ONLY when the corresponding analysis
# output ROOT file is missing.
#
# Return:
#   True  -> every ROOT file was readable and total events == 0
#   False -> input contains events OR could not be checked safely
#
# A failed/ambiguous edmFileUtil check therefore remains a real
# missing output and will NOT be silently ignored.
# ============================================================

def input_list_has_zero_events(list_path: Path):

    try:
        content = list_path.read_text()
    except OSError:
        return False

    root_files = [
        line.strip()
        for line in content.splitlines()
        if line.strip()
        and not line.strip().startswith("#")
    ]

    if not root_files:
        return False

    total_events = 0

    for root_file in root_files:

        try:
            result = subprocess.run(
                [
                    "edmFileUtil",
                    root_file,
                ],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                timeout=120,
            )

        except (
            subprocess.TimeoutExpired,
            OSError,
        ):
            return False

        if result.returncode != 0:
            return False

        text = (
            result.stdout
            + "\n"
            + result.stderr
        )

        matches = re.findall(
            r"(\d+)\s+events\b",
            text,
        )

        if not matches:
            return False

        total_events += int(
            matches[-1]
        )

    return total_events == 0


# ============================================================
# Check one sample
# ============================================================

def check_sample(
    input_base: Path,
    output_base: Path,
    sample: str,
):

    input_sample_dir = (
        input_base
        / sample
    )

    output_sample_dir = (
        output_base
        / sample
    )

    expected = collect_expected_outputs(
        input_sample_dir,
        output_sample_dir,
        sample,
    )

    if not expected:

        return {
            "sample": sample,
            "expected": [],
            "found": [],
            "empty_input": [],
            "missing": [],
            "output_dir": output_sample_dir,
        }

    found = []
    empty_input = []
    missing = []

    for item in expected:

        if (
            item["root_file"].is_file()
            and item["root_file"].stat().st_size > 0
        ):
            # Analysis output exists.
            #
            # Even if the analysis selection itself produced
            # zero selected events, this is still a legitimate
            # analysis output and must be included in hadd.
            found.append(item)

        else:

            # Analysis output is absent.
            #
            # Check whether the corresponding upstream
            # NtupleForge input itself contains zero events.
            # If so, no analysis output is expected.
            if input_list_has_zero_events(
                item["list_file"]
            ):
                empty_input.append(item)

            else:
                missing.append(item)

    return {
        "sample": sample,
        "expected": expected,
        "found": found,
        "empty_input": empty_input,
        "missing": missing,
        "output_dir": output_sample_dir,
    }


# ============================================================
# Write missing-file report
# ============================================================

def write_missing_report(
    result,
):

    sample = result["sample"]
    output_dir = result["output_dir"]

    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    report_file = (
        output_dir
        / f"{sample}_missing_check.log"
    )

    with report_file.open("w") as f:

        f.write(
            f"Sample   : {sample}\n"
        )

        f.write(
            f"Expected : {len(result['expected'])}\n"
        )

        f.write(
            f"Found    : {len(result['found'])}\n"
        )

        f.write(
            f"Missing  : {len(result['missing'])}\n"
        )

        f.write("\n")

        for item in result["missing"]:

            f.write(
                f"job={item['job_number']} "
                f"list={item['list_file']} "
                f"root={item['root_file']}\n"
            )

    return report_file


# ============================================================
# Run hadd
# ============================================================

def run_hadd(
    result,
    recreate=False,
    dry_run=False,
):

    sample = result["sample"]

    input_files = [
        item["root_file"]
        for item in result["found"]
    ]

    output_dir = (
        result["output_dir"]
    )

    output_file = (
        output_dir
        / f"{sample}.root"
    )

    log_file = (
        output_dir
        / f"{sample}_hadd.log"
    )

    cmd = [
        "hadd",
    ]

    if recreate:
        cmd.append("-f")

    cmd.append(
        str(output_file)
    )

    cmd.extend(
        str(path)
        for path in input_files
    )

    print()
    print(
        f"[HADD] {sample}"
    )

    print(
        f"       Inputs : {len(input_files)}"
    )

    print(
        f"       Output : {output_file}"
    )

    if dry_run:

        print(
            "       [DRY-RUN] hadd not executed."
        )

        return True

    if output_file.exists() and not recreate:

        print(
            "[SKIP] Merged output already exists:"
        )

        print(
            f"       {output_file}"
        )

        print(
            "       Use --recreate to overwrite it."
        )

        return True

    print(
        f"       Log    : {log_file}"
    )

    try:

        with log_file.open("w") as log:

            subprocess.run(
                cmd,
                stdout=log,
                stderr=subprocess.STDOUT,
                check=True,
            )

    except subprocess.CalledProcessError as exc:

        print(
            f"[ERROR] hadd failed for {sample}"
        )

        print(
            f"        return code = {exc.returncode}"
        )

        print(
            f"        log = {log_file}"
        )

        return False

    # --------------------------------------------------------
    # Validate merged output
    # --------------------------------------------------------

    if (
        not output_file.is_file()
        or output_file.stat().st_size == 0
    ):

        print(
            "[ERROR] Merged ROOT file was not created "
            "or is empty:"
        )

        print(
            f"        {output_file}"
        )

        return False

    size_mb = (
        output_file.stat().st_size
        / 1024.0
        / 1024.0
    )

    print(
        f"[OK] Merged successfully: "
        f"{size_mb:.1f} MB"
    )

    return True


# ============================================================
# Select samples
# ============================================================

def select_samples(
    args,
    input_base: Path,
):

    all_samples = sorted(
        path.name
        for path in input_base.iterdir()
        if path.is_dir()
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

    # Filter channel-dependent data samples
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
# Main
# ============================================================

def run(args):

    # --------------------------------------------------------
    # InputList base directory
    # --------------------------------------------------------

    input_base = Path(
        args.input_base
    ).expanduser()

    if not input_base.is_absolute():
        input_base = (
            Path.cwd()
            / input_base
        )

    input_base = input_base.resolve()

    if not input_base.is_dir():

        print(
            "[ERROR] InputList directory not found:"
        )

        print(
            f"        {input_base}"
        )

        sys.exit(1)

    # --------------------------------------------------------
    # Study output directory
    # --------------------------------------------------------

    output_base = (
        OUTPUT_BASE
        / args.study
        / args.run_period
        / args.channel
    )

    if not output_base.is_dir():

        print(
            "[ERROR] Analysis output directory "
            "does not exist:"
        )

        print(
            f"        {output_base}"
        )

        sys.exit(1)

    # --------------------------------------------------------
    # Select samples
    # --------------------------------------------------------

    samples = select_samples(
        args,
        input_base,
    )

    print("=" * 78)
    print("SSB ROOT HADD")
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
        f"Output base : {output_base}"
    )

    print(
        f"Samples     : {len(samples)}"
    )

    print(
        f"Check only  : {args.check_only}"
    )

    print(
        f"Dry run     : {args.dry_run}"
    )

    print("=" * 78)

    total_expected = 0
    total_found = 0
    total_empty_input = 0
    total_missing = 0

    merge_ok = 0
    merge_failed = 0
    merge_skipped = 0

    # --------------------------------------------------------
    # Process samples
    # --------------------------------------------------------

    for sample in samples:

        if not (
            input_base
            / sample
        ).is_dir():

            print(
                f"[WARNING] Input sample not found: "
                f"{sample}"
            )

            continue

        result = check_sample(
            input_base,
            output_base,
            sample,
        )

        expected_count = len(
            result["expected"]
        )

        found_count = len(
            result["found"]
        )

        empty_input_count = len(
            result["empty_input"]
        )

        missing_count = len(
            result["missing"]
        )

        total_expected += expected_count
        total_found += found_count
        total_empty_input += empty_input_count
        total_missing += missing_count

        if expected_count == 0:

            print(
                f"{sample:<50} "
                "NO INPUT LIST"
            )

            continue

        print(
            f"{sample:<50} "
            f"{found_count:>5}/{expected_count:<5} "
            f"FOUND",
            end="",
        )

        if empty_input_count:

            print(
                f"   EMPTY_INPUT={empty_input_count}",
                end="",
            )

        if missing_count:

            print(
                f"   MISSING={missing_count}"
            )

            report_file = (
                write_missing_report(
                    result
                )
            )

            print(
                f"    -> skip hadd"
            )

            print(
                f"    -> report: "
                f"{report_file}"
            )

            merge_skipped += 1

            continue

        if empty_input_count:
            print("   OK")
        else:
            print("   OK")

        # ----------------------------------------------------
        # Check-only mode
        # ----------------------------------------------------

        if args.check_only:

            continue

        # ----------------------------------------------------
        # Hadd
        # ----------------------------------------------------

        success = run_hadd(
            result,
            recreate=args.recreate,
            dry_run=args.dry_run,
        )

        if success:
            merge_ok += 1
        else:
            merge_failed += 1

    # --------------------------------------------------------
    # Summary
    # --------------------------------------------------------

    print()
    print("=" * 78)
    print("SUMMARY")
    print("=" * 78)

    print(
        f"Expected chunk files : {total_expected}"
    )

    print(
        f"Found chunk files    : {total_found}"
    )

    print(
        f"Empty input files    : {total_empty_input}"
    )

    print(
        f"Missing chunk files  : {total_missing}"
    )

    if not args.check_only:

        print(
            f"Merged / ready       : {merge_ok}"
        )

        print(
            f"Merge failed         : {merge_failed}"
        )

        print(
            f"Merge skipped        : {merge_skipped}"
        )

    print("=" * 78)

    if (
        total_missing > 0
        or merge_failed > 0
    ):
        sys.exit(1)

    sys.exit(0)


# ============================================================
# CLI
# ============================================================

def parse_args():

    parser = argparse.ArgumentParser(
        description=(
            "Check and hadd SSBNanoAOD Condor outputs."
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
        "--input-base",
        required=True,
        help=(
            "Base directory containing the sample "
            "directories with per-job InputList .list files."
        ),
    )

    parser.add_argument(
        "--samples",
        nargs="+",
        default=["all"],
        help=(
            "Samples to merge. "
            "Use all, data, mc, or explicit sample names. "
            "Default: all"
        ),
    )

    parser.add_argument(
        "--check-only",
        action="store_true",
        help=(
            "Check completeness only; "
            "do not run hadd."
        ),
    )

    parser.add_argument(
        "--dry-run",
        action="store_true",
        help=(
            "Show which samples would be merged, "
            "but do not execute hadd."
        ),
    )

    parser.add_argument(
        "--recreate",
        action="store_true",
        help=(
            "Overwrite an existing merged "
            "<sample>.root file using hadd -f."
        ),
    )

    return parser.parse_args()


# ============================================================
# Entry point
# ============================================================

if __name__ == "__main__":

    args = parse_args()

    run(args)
