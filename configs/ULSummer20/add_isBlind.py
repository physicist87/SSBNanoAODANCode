#!/usr/bin/env python3
# Inserts a Blinding section into all .config files under UL* subdirectories.
# Placement: immediately after the Luminosities line.
# Skips files that already have an isBlind entry.

import os
import glob

BLIND_SECTION = "\n####### Blinding #######\nisBlind : \"True\" ### after unblinding False\n"
KEY = "isBlind"

BASE_DIRS = [
    "UL2016PreVFP",
    "UL2016PostVFP",
    "UL2017",
    "UL2018",
]

def add_isblind(config_path):
    with open(config_path, "r") as f:
        lines = f.readlines()

    # Skip if isBlind already present
    for line in lines:
        if line.strip().startswith(KEY):
            print(f"  [skip] already has isBlind: {config_path}")
            return

    # Find insertion point: last uncommented Luminosities line
    insert_after = -1
    for i, line in enumerate(lines):
        stripped = line.strip()
        if not stripped.startswith("#") and stripped.startswith("Luminosities"):
            insert_after = i

    # Fallback: last uncommented RunRange line
    if insert_after == -1:
        for i, line in enumerate(lines):
            stripped = line.strip()
            if not stripped.startswith("#") and "RunRange" in stripped:
                insert_after = i

    if insert_after == -1:
        with open(config_path, "a") as f:
            f.write(BLIND_SECTION)
        print(f"  [done/appended] {config_path}")
        return

    lines.insert(insert_after + 1, BLIND_SECTION)

    with open(config_path, "w") as f:
        f.writelines(lines)
    print(f"  [done] {config_path}  (after line {insert_after + 1}: {lines[insert_after].rstrip()})")

def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    for base in BASE_DIRS:
        pattern = os.path.join(script_dir, base, "*.config")
        files = sorted(glob.glob(pattern))
        if not files:
            print(f"[warn] no .config files found in {base}/")
            continue
        print(f"\n[{base}]")
        for f in files:
            add_isblind(f)

if __name__ == "__main__":
    main()
