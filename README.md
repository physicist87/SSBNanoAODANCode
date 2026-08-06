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
git clone --branch Run2_ULSummer20_v6p4 https://github.com/physicist87/SSBNanoAODANCode.git
cd SSBNanoAODANCode
```

If this repository is forked to your account:
```sh
git clone --branch Run2_ULSummer20_v6p4 https://github.com/<username>/SSBNanoAODANCode.git
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

Only UL2018 configs/branch lists exist so far; UL2016PreVFP/UL2016PostVFP/UL2017 are not yet migrated. See `v15_migration/README.md` for the full "what changed, by area" writeup and `v15_migration/NOTES.md` for the complete, dated account of every bug found and fixed along the way (including two real physics bugs: a JEC scale factor being squared into jet mass, and Type-1 MET being built from an already-corrected MET branch instead of the genuinely raw one).

## Object-Oriented Refactor (In Progress)

Alongside the v15 migration, the codebase is being incrementally refactored toward clearer separation of concerns - each step compiled and run against real data before moving to the next, to keep risk low against an already-verified physics pipeline. Completed so far:

1. **Stringly-typed dispatch -> enum + shared helper**: a `BTagAlgo` enum (`interface/SSBCorrections.h`) replaced two independently-maintained string-parsing blocks (in `SSBCorrections.cpp` and `Analysis.cpp`) that had already drifted out of sync once in practice. A generic `LoadOptionalCorrection<T>()` template replaced 6 near-identical try/catch blocks for optional `correctionlib` loads.
2. **`JetID` class** (`interface/JetID.h`, `src/JetID.cpp`): the NanoAODv15 PUPPI jet ID logic, pulled out of `Analysis` (see above).
3. **`Jet` value object** (`interface/Jet.h`): `RawJet`/`CorrT1Jet` structs consolidating `Analysis::MakeJetCollection()`'s per-jet raw reads into a single pass, instead of scattered/duplicated bounds-checked array access. `SSBCorrections`' own jet-correction/MET function signatures were deliberately left unchanged in this step, since that code is the most carefully tutorial-verified part of the migration.

Planned next: splitting `SSBCorrections` into per-domain classes (JEC/JER, Type-1 MET, b-tag SF, PU jet ID, lepton SF) behind a facade that preserves its current public interface, and parametrizing systematic variation (nominal/up/down) instead of ad hoc flags. See `v15_migration/NOTES.md` items 26-28 for the full rationale behind each step.

## Notes
- Ensure that ROOT is properly installed and configured before running the analysis.
- Modify `Makefile` if necessary to match your system's compiler settings.
- Use the configuration files in `configs/` to customize your analysis workflow.

## Contact
For questions or contributions, please open an issue or contact the maintainers.

---
This README provides essential instructions for setting up and running the SSB NanoAOD analysis. If you need additional details, feel free to modify and expand it!
