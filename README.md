# SSB NanoAOD Analysis Code

This repository contains code for analyzing CMS NanoAOD data using the SSB framework. The analysis processes event data and extracts relevant physics observables.

## Installation and Setup

## Correctionlib Integration

This project uses [correctionlib](https://cms-nanoaod.github.io/correctionlib/) for JEC, JER, and lepton scale factor corrections.

The current packages utilize correctionlib, and if you are using ROOT through CMSSW, the following procedures are not necessary.
You can simply install it by using a release like CMSSW_13_x_x.

The Makefile is configured to automatically detect the appropriate correctionlib installation using the following priority:

1. **Environment Variables**:  
   If `CORRECTION_PATH` and `CORRECTION_LIBPATH` are defined, they are used for headers and libraries.
2. **CVMFS Installation**:  
   Falls back to a known working CVMFS release:  
   `/cvmfs/sft.cern.ch/lcg/releases/correctionlib/2.2.2-13853/x86_64-centos7-gcc11-opt`
3. **pip-installed fallback** (e.g. on macOS):  
   Uses `correction config --cflags` and `--ldflags --rpath` if no other paths are available.

### Example for Local Build
If you have built correctionlib locally or installed via pip:
```bash
export CORRECTION_PATH=$HOME/mycorrectionlib/include
export CORRECTION_LIBPATH=$HOME/mycorrectionlib/lib
make -f Makefile_ssb
```
### Cloning the Repository
To get started, clone the repository from GitHub and check out the required branch:
```sh
git clone --branch Run2_ULSummer20_v7 https://github.com/physicist87/SSBNanoAODANCode.git
cd SSBNanoAODANCode
```

If this repository is forked to your account:
```sh
git clone --branch Run2_ULSummer20_v7 https://github.com/<username>/SSBNanoAODANCode.git
cd SSBNanoAODANCode
```

### Compilation Instructions
To compile the analysis program, use the provided `Makefile`:
```sh
make -f Makefile_ssb
```
This will generate the executable **`ssb_analysis`**.

## Repository Structure and Code Summary

### **Main Components**
- **`main_ssb.cpp`**: The main execution file that initializes the analysis and handles input data.
- **`src/Analysis.cpp`**: Implements core analysis functions, including event selection and variable calculations.
- **`interface/Analysis.h`**: Defines the analysis class, including member functions and variables used for event processing.
- **`interface/NanoAODBranchReader.h` / `src/NanoAODBranchReader.cpp`**: Owns the dynamic branch-name -> `TTreeReader` bindings for a NanoAOD `TChain` and reads them back out in a way that's resilient to storage-type changes and renames between NanoAOD versions (added for the v15 migration - see below). `Analysis` calls through this for all branch I/O instead of touching `TTreeReader` maps directly.
- **`interface/JetID.h` / `src/JetID.cpp`**: The NanoAODv15 PUPPI jet ID (tight / tight-lepton-veto), reimplemented directly from jet energy fractions since the `Jet_jetId` bitmask no longer exists for Puppi jets. Extracted out of `Analysis` as part of an ongoing object-oriented refactor (see below).
- **`interface/Jet.h`**: Small value-object structs (`RawJet`, `CorrT1Jet`) bundling a jet's per-event raw NanoAOD fields together, used by `Analysis::MakeJetCollection()` to read the `Jet_`/`CorrT1METJet_` collections in one pass instead of scattered per-field reads.
- **`interface/SSBCorrections.h` / `src/SSBCorrections.cpp`**: Loads and applies all `correctionlib`-based corrections - JEC/JER/JVM, Type-1 MET, b-tagging SF, PU jet ID SF, muon/electron/trigger SF, PU weight, MET-XY correction.
- **`CommonTools.cpp/hpp`**: Provides utility functions used across multiple parts of the analysis.
- **`TextReader/TextReader.cpp/hpp`**: Handles reading text-based inputs for configurations or dataset lists.
- **`configs/`**: Contains various configuration files for different datasets and analysis setups.
- **`input/`**: Stores input dataset lists, specifying ROOT files to be analyzed.
- **`output/`**: Contains processed ROOT files after running the analysis.
- **`xsecAndsample/`**: Holds cross-section and sample information for different datasets.
- **`branchlist/`**: Contains lists of branches that are used in the analysis, ensuring compatibility with different versions of NanoAOD.

## Running the Analysis
### Input Files
Input file lists are stored in the `input/` directory. Example:
- `input/UL2016PreVFP/TTbar_Signal_1.list`

### Execution Format
Based on `run_ssb_check.sh`, the program should be executed with the following arguments:
```sh
./ssb_analysis <runPeriod> <StudyName> <Channels> <input_list>
```
Where:
- **`runPeriod`**: Specifies the data-taking period (e.g., `UL2016PreVFP`, `UL2017`, `UL2018`).
- **`StudyName`**: Defines the study version (e.g., `Testv1`).
- **`Channels`**: The physics channel being analyzed:
  - `MuMu` → Dimuon channel -> Configfile : dimuon.config
  - `ElEl` → Dielectron channel -> Configfile: dielec.config
  - `MuEl` → Muon-electron channel -> Configfile: muelec.config
- **`input_list`**: The dataset file (e.g., `input/UL2016PreVFP/TTbar_Signal_1.list`).

### Example Execution Commands
For a dimuon analysis in UL2016 PreVFP:
```sh
./ssb_analysis "${runPeriod}/TTbar_Signal/${i}.list" "${StudyName}/${runPeriod}/${Channels}/TTbar_Signal/${i}.root" "ULSummer20/UL2016PreVFP/dimuon.config" "None" ${runPeriod} -1
```
For a dielectron analysis in UL2017:
```sh
./ssb_analysis UL2017 StudyX ElEl input/UL2016PreVFP/TTbar_Signal_1.list
./ssb_analysis "${runPeriod}/TTbar_Signal/${i}.list" "${StudyName}/${runPeriod}/${Channels}/TTbar_Signal/${i}.root" "ULSummer20/UL2016PreVFP/dielec.config" "None" ${runPeriod} -1
```
For a muon-electron analysis in UL2018:
```sh
./ssb_analysis "${runPeriod}/TTbar_Signal/${i}.list" "${StudyName}/${runPeriod}/${Channels}/TTbar_Signal/${i}.root" "ULSummer20/UL2016PreVFP/muelec.config" "None" ${runPeriod} -1

```

### Running a Check
To validate the setup and check if the program is correctly configured, use:
```sh
bash run_ssb_check.sh
```
Alternatively, execute it directly:
```sh
./run_ssb_check.sh
```

### Checking the Output
The processed ROOT files will be stored in the `output/` directory. Example:
```sh
ls output/Testv1/UL2016PreVFP/MuMu/
```
Expected output:
```
UL2016PreVFP_TTbar_Signal_1.root
TTbar_Signal_1.root
```

## Branchlist Directory and Its Role
The `branchlist/` directory contains text files that list the specific branches read from NanoAOD files. These files ensure compatibility with different NanoAOD versions by specifying which branches to extract. Each line in `branchlist/branch_list.txt` follows this format:
```
<Branch Name>, <Category>, <Data Type>, <Structure>
```
Where:
- **Branch Name**: The name of the branch in the ROOT file.
- **Category**: Indicates the type of physics object (e.g., `trigger`, `muon`, `electron`).
- **Data Type**: The type of data stored in the branch (e.g., `Bool_t`, `UInt_t`, `Float_t`).
- **Structure**: Specifies whether it is a `single` value or a `vector`.

### **Example Entries**
```
HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ, trigger, Bool_t, single
HLT_IsoMu24, trigger, Bool_t, single
nMuon, muon, UInt_t, single
Muon_pfRelIso03_all, muon, Float_t, vector
Muon_eta, muon, Float_t, vector
Muon_mass, muon, Float_t, vector
```
This structure ensures that only the necessary branches are accessed, optimizing performance and reducing memory usage.

For NanoAODv15 running, use `branchlist/<era>/branch_list_v15.txt` (e.g. `branchlist/UL2018/branch_list_v15.txt`) instead of the v9 `branch_list.txt` - see the next section.

## NanoAODv15 Migration

This branch adds support for **NanoAODv15** input (the dilepton `MuMu`/`ElEl`/`MuEl` channels) alongside the existing v9 code path, without requiring a separate copy of the analysis code. It has been verified end-to-end against a real UL2018 `TTbar_Signal` NanoAODv15 file: clean compile, full event-loop run with no crashes, all corrections loading and evaluating successfully.

The single biggest driver of change is that NanoAODv15's default AK4 jet collection switched from PF+CHS (`AK4PFchs`) to **Puppi** (`AK4PFPuppi`), which cascades into JEC naming, JER, b-tagging, and the whole MET propagation chain. Highlights:

- **Version-agnostic branch reading**: a new `NanoAODBranchReader` class asks the `TChain` directly for each branch's real type/structure instead of trusting a branch list's declared type, so the same code handles both v9 and v15 storage-type changes (e.g. several `Int_t`/`UInt_t` branches narrowed to `Short_t`/`UChar_t`) without per-version edits.
- **Jet ID**: `Jet_jetId` no longer exists for Puppi jets - the POG-recommended tight/tight-lepton-veto working points are reimplemented directly from jet energy fractions (`Jet_neHEF`/`neEmEF`/`chHEF`/`chEmEF`/`muEF`/multiplicities).
- **b-tagging**: switched from DeepCSV/DeepJet to **UParT** (Unified Particle Transformer), including automatic derivation of the numeric discriminant cut from the working-point letter alone (`Jet_btag = "UParTM"` is enough - no separate `BTagDiscCut` config value required).
- **JEC/JER/Type-1 MET**: rebuilt to match the CMS JERC ApplicationTutorial term-for-term, including starting Type-1 MET from the genuinely uncorrected `RawMET_pt`/`RawPuppiMET_pt` (not the already-Type1-corrected `MET_pt`/`PuppiMET_pt`), muon-subtracting jets before Type-1 propagation, and including the low-pT `CorrT1METJet_` collection.
- **cvmfs paths**: default correction-file distribution switched to CMS's new "CAT" layout (`/cvmfs/cms-griddata.cern.ch/cat/metadata/...`), with per-POG campaign-folder differences documented (not every POG needed a new v15-era folder).

All four Run-2 UL run periods now have v15 configs/branch lists: `UL2018`, `UL2017`, `UL2016PreVFP`, `UL2016PostVFP`. See `v15_migration/README.md` for the full "what changed, by area" writeup and `v15_migration/NOTES.md` for the complete, dated account of every bug found and fixed along the way (including two real physics bugs: a JEC scale factor being squared into jet mass, and Type-1 MET being built from an already-corrected MET branch instead of the genuinely raw one).

### Per-era JEC/JER/JetVeto tags (confirmed against actual `jet_jerc.json.gz` / `jetvetomaps.json.gz`)

| Era | JECName | JERName | JetVetoName |
|---|---|---|---|
| UL2018 | `Summer20UL18NanoV15_V1` | `Summer19UL18_JRV3_MC_ScaleFactor_AK4PFPuppi` | `Summer19UL18_V1` |
| UL2017 | `Summer20UL17NanoV15_V1` | `Summer19UL17_JRV4_MC_ScaleFactor_AK4PFPuppi` | `Summer19UL17_V1` |
| UL2016PreVFP | `Summer20UL16APVNanoV15_V1` | `Summer20UL16APV_JRV5_MC_ScaleFactor_AK4PFPuppi` | `Summer19UL16_V1` |
| UL2016PostVFP | `Summer20UL16NanoV15_V1` | `Summer20UL16_JRV5_MC_ScaleFactor_AK4PFPuppi` | `Summer19UL16_V1` |

For UL2017, the JME twiki's published "Recommended maps" table lists a `V2` jet-veto map; the actual `jetvetomaps.json.gz` shipped in this copy of `jsonpog-integration` only contains `V1`, so the config is deliberately set to `V1` with a comment explaining the discrepancy. Re-check if a newer `jetvetomaps.json.gz` becomes available.

`Jet_btag` is set to `UParTM` (UParT Medium WP) everywhere; the numeric cut is not hardcoded in the configs - it's looked up automatically at runtime from `btagging.json.gz`'s `UParTAK4_wp_values` correction (confirmed e.g. `0.161` for UL2018 MC in a real job log). `BTagSFType` is `comb` for all eras, since UParT only ships a combined heavy-flavor SF (no separate `mujets` variant like the old DeepJet WP had).

PU/JMAR/Muon/Electron SF `*Path` keys are intentionally left on their old v9-style paths for every era (not yet verified against a v15-restructured directory) - flagged with an inline comment in each config. JEC/JER/JetVeto paths, by contrast, are confirmed against real correction files.

### Era-specific trigger / branch-list splits

Some run periods need more than one branch list + config because the HLT menu changed mid-era:

- **UL2017**: split into `Run2017B` (no `*_Mass3p8` dimuon trigger; also needs `HLT_IsoMu24_eta2p1` as the single-muon trigger) and `Run2017C-F` (has `*_Mass3p8`). Branch lists: `branchlist/UL2017/branch_list_Run2017B_v15.txt`, `branchlist/UL2017/branch_list_Run2017CtoF_v15.txt` (MC/other-data use `branchlist/UL2017/branch_list_v15.txt`, the union of both trigger sets).
- **UL2016PostVFP**: split into `Run2016H` (needs DZ-variant dimuon/dielectron/muon-electron triggers, and drops `HLT_DoubleEle33_CaloIdL_GsfTrkIdVL`) and the base/other eras. Branch lists: `branchlist/UL2016PostVFP/branch_list_RunH_v15.txt` vs `branchlist/UL2016PostVFP/branch_list_v15.txt`.
- **UL2016PreVFP**: no additional era split needed - one branch list/config set covers the whole era.

`CondorJobs/submit_jobs_v1.py`'s `get_job_config()` and the local interactive `run_v4.sh` both implement this era routing at submission/test time (matching on the run-era substring of the data sample name, e.g. `Run2017B`, `Run2016H`); MC samples always use the base/union branch list for their run period.

## Object-Oriented Refactor (In Progress)

Alongside the v15 migration, the codebase is being incrementally refactored toward clearer separation of concerns - each step compiled and run against real data before moving to the next, to keep risk low against an already-verified physics pipeline. Completed so far:

1. **Stringly-typed dispatch -> enum + shared helper**: a `BTagAlgo` enum (`interface/SSBCorrections.h`) replaced two independently-maintained string-parsing blocks (in `SSBCorrections.cpp` and `Analysis.cpp`) that had already drifted out of sync once in practice. A generic `LoadOptionalCorrection<T>()` template replaced 6 near-identical try/catch blocks for optional `correctionlib` loads.
2. **`JetID` class** (`interface/JetID.h`, `src/JetID.cpp`): the NanoAODv15 PUPPI jet ID logic, pulled out of `Analysis` (see above).
3. **`Jet` value object** (`interface/Jet.h`): `RawJet`/`CorrT1Jet` structs consolidating `Analysis::MakeJetCollection()`'s per-jet raw reads into a single pass, instead of scattered/duplicated bounds-checked array access. `SSBCorrections`' own jet-correction/MET function signatures were deliberately left unchanged in this step, since that code is the most carefully tutorial-verified part of the migration.

Planned next: splitting `SSBCorrections` into per-domain classes (JEC/JER, Type-1 MET, b-tag SF, PU jet ID, lepton SF) behind a facade that preserves its current public interface, and parametrizing systematic variation (nominal/up/down) instead of ad hoc flags. See `v15_migration/NOTES.md` items 26-28 for the full rationale behind each step.

## HTCondor Workflow

The current Condor workflow consists of three main utilities under `CondorJobs/`:

1. `submit_jobs_v1.py` — submit analysis jobs
2. `check_jobs_v1.py` — validate job completion and identify failed jobs
3. `run_hadd_v8.py` — check output completeness and merge ROOT files

```text
InputList
   |
   v
submit_jobs_v1.py
   |
   v
HTCondor jobs
   |
   v
check_jobs_v1.py
   |
   +---- bad jobs ----> resubmit
   |                       |
   |<----------------------+
   |
   v
run_hadd_v8.py --check-only
   |
   v
run_hadd_v8.py
   |
   v
<sample>.root
```

### 1. Job submission (`submit_jobs_v1.py`)

Creates the analysis package, JDL files, queue files, and submits jobs to HTCondor. Two modes: normal submission, and bad-job resubmission.

```bash
python3 CondorJobs/submit_jobs_v1.py \
    --study <STUDY> \
    --run-period <RUN_PERIOD> \
    --channel <CHANNEL> \
    --input-base <INPUT_LIST_BASE> \
    --samples <SAMPLES>
```

- Run periods: `UL2016PreVFP`, `UL2016PostVFP`, `UL2017`, `UL2018`
- Channels: `ElEl`, `MuEl`, `MuMu`
- Sample selectors: `all`, `data`, `mc`, or explicit sample names

Example (UL2018 MuMu, all samples):

```bash
INPUT2018=/path/to/InputList/2018
python3 CondorJobs/submit_jobs_v1.py \
    --study NanoAODv15_v3 --run-period UL2018 --channel MuMu \
    --input-base "$INPUT2018" --samples all
```

Useful flags:
- `--test` — submit only the first InputList file of each selected sample
- `--dry-run` — build the tarball/JDL/queue files but skip `condor_submit` (good for sanity-checking job config before submission)

Before submission the script builds `CondorJobs/SSBNanoAODANCode.tar.gz` from the current package, excluding `input/`, `output/`, `CondorJobs/`, `.git/`, `*.o`, `ssb_analysis`, `.DS_Store`, `__pycache__/`, `*.pyc`. The tarball and the per-job InputList are transferred to the worker node.

**Condor resources per job** (current):

```text
RequestCpus   = 1
RequestMemory = 4 GB
RequestDisk   = 10 GB
+JobType      = "short"
```

ROOT output files are *not* returned via Condor's output-transfer mechanism — `run_condor_v1.sh` stages them directly to the storage element.

JDLs/queue files land under `CondorJobs/condorSubmit/<study>/<run-period>/<channel>/`; logs under `CondorJobs/condorLog/<study>/<run-period>/<channel>/<sample>/` as `<sample>_<job>.{log,out,err}`. Resubmission reuses the same log filenames, so `check_jobs_v1.py` always looks at the *latest* Condor termination record for a job's status.

**Data-sample filtering by channel:**

| Channel | 2016/2017 | 2018 |
|---|---|---|
| `MuMu` | `SingleMuon`, `DoubleMuon` | same |
| `ElEl` | `SingleElectron`, `DoubleEG` | `EGamma`, `SingleElectron`, `DoubleEG` |
| `MuEl` | `SingleMuon`, `SingleElectron`, `MuonEG` | `SingleMuon`, `EGamma`, `SingleElectron`, `MuonEG` |

MC samples are not channel-filtered. See the "Era-specific trigger / branch-list splits" subsection above for how `UL2017`/`UL2016PostVFP` data route to the right config + branch list.

### 2. Checking jobs (`check_jobs_v1.py`)

Cross-checks expected jobs (from the InputList) against Condor logs, worker stdout/stderr, and the ROOT outputs on the SE.

```bash
python3 CondorJobs/check_jobs_v1.py \
    --study NanoAODv15_v3 --run-period UL2018 --channel MuMu --samples all
```

Job statuses:

- **`OK`** — analysis destructor completed, worker wrapper completed, ROOT output exists on the SE, latest Condor exit code is 0.
- **`FAILED`** — non-zero Condor exit code or a serious stderr pattern (`segmentation violation/fault`, `Traceback`, `[ERROR]`, `fatal error`, `std::exception`, `terminate called`, `Aborted`, `core dumped`). Ordinary warnings don't count.
- **`STAGEOUT_FAIL`** — analysis completed but the worker wrapper didn't reach normal completion (worker/stage-out problem).
- **`OUTPUT_ONLY`** — ROOT output exists on the SE but local Condor logs are incomplete.
- **`RUNNING`** — Condor log exists but has no termination record yet.
- **`INCOMPLETE`** — some local logs exist but the job can't be classified as complete.
- **`MISSING`** — no valid output and no useful Condor job info.
- **`EMPTY_INPUT`** — the expected output is missing, but every ROOT file in that job's InputList was successfully checked with `edmFileUtil` and the total input event count is exactly zero. **Not treated as a failure.**

`EMPTY_INPUT` exists because NtupleForge can legitimately produce a zero-event output (e.g. every lumisection in an input file falls outside the Golden JSON). A job is only classified `EMPTY_INPUT` if *every* input file in its list can be read and the total is exactly zero; if any file can't be checked, the job stays a genuine bad/unfinished job — this conservative behavior stops real input-access problems from being hidden as `EMPTY_INPUT`.

Useful flags: `--details` (print every job, not just bad ones), `--show-errors` (print the serious stderr lines found), `--write-bad` (write `CondorJobs/bad_jobs_<study>_<run-period>_<channel>.txt` with `SAMPLE JOB_NUMBER STATUS INPUT_LIST_PATH` per line; `OK` and `EMPTY_INPUT` jobs are excluded).

### 3. Resubmitting failed jobs

Feed the bad-job file straight back into `submit_jobs_v1.py`:

```bash
python3 CondorJobs/submit_jobs_v1.py \
    --study NanoAODv15_v3 --run-period UL2018 --channel MuMu \
    --bad-jobs CondorJobs/bad_jobs_NanoAODv15_v3_UL2018_MuMu.txt
```

`--input-base` isn't needed here since the bad-job file already has the full InputList path. Re-run `check_jobs_v1.py` afterward (it reuses the latest termination record automatically) and repeat until `Bad / unfinished : 0`.

### 4. Checking and merging outputs (`run_hadd_v8.py`)

```bash
python3 CondorJobs/run_hadd_v8.py \
    --study NanoAODv15_v3 --run-period UL2018 --channel MuMu \
    --input-base "$INPUT2018" --samples all --check-only
```

For every `<sample>_<job>.list` in the InputList, expects a matching `<sample>_<job>.root` in the output directory (`/pnfs/knu.ac.kr/data/cms/store/user/${USER}/CPV_Run2/ULSummer20/<study>/<run-period>/<channel>/<sample>/`). Missing outputs are re-checked with `edmFileUtil`: if every corresponding input file has zero events, the chunk is `EMPTY_INPUT` and excluded from hadd (no output expected); otherwise it stays `MISSING` and hadd is skipped for that sample. An existing output file is always treated as a valid chunk, even if it contains zero selected events after the analysis cuts.

Missing chunks (genuine) get a report written to `<sample>_missing_check.log` next to the output, with expected/found/missing counts and the missing job numbers/paths.

Once the check is clean, drop `--check-only` to actually run `hadd`:

```bash
python3 CondorJobs/run_hadd_v8.py \
    --study NanoAODv15_v3 --run-period UL2018 --channel MuMu \
    --input-base "$INPUT2018" --samples all
```

Output: `<sample>/<sample>.root` + `<sample>/<sample>_hadd.log`. Already-merged samples are skipped unless `--recreate` is given (runs `hadd -f`). `--dry-run` shows what would be merged without running `hadd`.

### 5. Recommended production loop

```text
submit (all) -> check (--write-bad) -> resubmit bad jobs -> check again -> ... -> 0 bad
   -> hadd --check-only (confirm 0 genuinely missing) -> hadd
```

### Notes on the three scripts

- `submit_jobs_v1.py` and `run_hadd_v8.py` take the InputList base via `--input-base`; `check_jobs_v1.py` currently has its InputList location hardcoded in the script itself (not yet unified to `--input-base` - a possible future cleanup).
- `EMPTY_INPUT` is a legitimate upstream zero-event condition, not an analysis failure, and is handled consistently by both `check_jobs_v1.py` and `run_hadd_v8.py`.
- `run_hadd_v8.py` never silently assumes a missing output is a zero-event case - if it can't verify that with `edmFileUtil`, the chunk stays `MISSING` and blocks hadd for that sample.

## Recent Fixes

- **`SetUpKINObs()` degenerate-solution guard was killing every good KinReco solution** (`src/Analysis.cpp`): the `isKinSol = false;` line sat outside the `if (Top.Energy() < 0.01)` braces, so it ran unconditionally whenever `isKinSol` had just been set `true` - discarding every successful kinematic reconstruction, not just genuinely degenerate (near-zero-energy) ones. This made `h_Lep1pt_8` (AfterTopReconstruction) and the CP-observable histogram `h_Reco_CPO1_ReRange` empty for every sample. Fixed by requiring `isKinSol && Top.Energy() < 0.01` (and the `AnTop` equivalent) before resetting. Verified on a local UL2018 `TTbar_Signal` test: `h_Lep1pt_8` went from 0 to 15985 entries, consistent with the expected KinReco efficiency relative to `h_Lep1pt_5` (AfterBTagMultiplicity, 17689 entries).
- **UL2017 `JetVetoName` reverted `V2 -> V1`**: an earlier change (based on the JME twiki's published "Recommended maps" table) set `Summer19UL17_V2`, but the actual `jetvetomaps.json.gz` in this copy of `jsonpog-integration` only defines `Summer19UL17_V1`. Reverted across all 11 UL2017 configs with a comment explaining the discrepancy.

## Notes
- Ensure that ROOT is properly installed and configured before running the analysis.
- Modify `Makefile` if necessary to match your system's compiler settings.
- Use the configuration files in `configs/` to customize your analysis workflow.

## Contact
For questions or contributions, please open an issue or contact the maintainers.

---
This README provides essential instructions for setting up and running the SSB NanoAOD analysis. If you need additional details, feel free to modify and expand it!
