# DY Estimation (R_out/in) Workflow

End-to-end recipe for the data-driven Drell-Yan normalization with the
`DYEstimation` mode of this branch: Condor production of Data + DY MC,
job checking and resubmission, hadd, and the R_out/in calculation with
`CalDYSF_Run2UL_v15.C`.

What the `DYEstimation` flag changes in `Analysis` is described in the main
`README.md` ("DY Estimation Mode"). In short: no Z-veto before the
jet/MET/b-tag steps, MET > 40 GeV in all channels (incl. `muel`), and the
extra histograms `h_cf_DiLepMass_<step>`, `h_NoMETCut_*_<step>` and
`h_cf_NoMETCut_DiLepMass_<step>`.

Branch: `Run2_ULSummer20_v7_DYEstimation` (branched from `Run2_ULSummer20_v7`).

## Files

| File | Purpose |
| --- | --- |
| `configs/ULSummer20/<era>/*_DY.config` | Copies of the normal configs with `DYEstimation : "True"` |
| `CondorJobs/submit_jobs_v1.py --dy-estimation` | Same era routing, picks the `*_DY.config` twin |
| `DYEstimation/CalDYSF_Run2UL_v15.C` | R_out/in, k_ll and DY SF for all Run 2 UL eras |

## 1. Environment (cluster)

```bash
git clone --branch Run2_ULSummer20_v7_DYEstimation https://github.com/physicist87/SSBNanoAODANCode.git
cd SSBNanoAODANCode                   # the directory name must be exactly SSBNanoAODANCode
cmsenv
make -f Makefile_ssb
BASE=/u/user/sha/Develop/CPviolation/SSB/AnalysisCode/NanoAODNtuple_v1/NtupleList_v1/NanoAODv15_v2/InputList
STUDY=NanoAODv15_DYEst_v1
```

`submit_jobs_v1.py` tars the package as `SSBNanoAODANCode` and the worker
expects that name, so a clone under another directory name fails at the
tarball step.

## 2. Submit Data + DY only

Explicit sample names are used because `--samples data` ignores any other
name given with it. Data samples are still filtered per channel by the
script.

```bash
for ERA in 2016PreVFP 2016PostVFP 2017 2018; do
  SAMPLES=$(ls $BASE/$ERA | grep -E '^(Data_|DYJetsToLL_)')
  for CH in MuMu ElEl MuEl; do
    python3 CondorJobs/submit_jobs_v1.py --study $STUDY --run-period UL$ERA --channel $CH \
        --input-base $BASE/$ERA --samples $SAMPLES --dy-estimation
  done
done
```

Run once with `--dry-run` first and check that every JDL uses a `_DY` config:

```bash
grep -l "Arguments" -r CondorJobs/condorSubmit/$STUDY | xargs grep -L "_DY.config"   # must print nothing
```

Jobs of this study in the queue (the JDLs use `Arguments`, so `-af Args` is always empty):

```bash
condor_q $USER -af Arguments | grep -c "$STUDY"
```

## 3. Check jobs

Only after the queue for this study is empty (`--write-bad` also lists
jobs that are still idle or running):

```bash
for ERA in 2016PreVFP 2016PostVFP 2017 2018; do
  SAMPLES=$(ls $BASE/$ERA | grep -E '^(Data_|DYJetsToLL_)')
  for CH in MuMu ElEl MuEl; do
    echo "=== UL$ERA $CH"
    python3 CondorJobs/check_jobs_v1.py --study $STUDY --run-period UL$ERA --channel $CH \
        --samples $SAMPLES --write-bad | grep -E '^[A-Za-z_/ ]+ +: [0-9]+$'
  done
done
```

Pass the same sample list as at submission; `--samples all` would report
every MC sample that was never submitted as MISSING.

### Failure modes seen in production

| Symptom | Cause | Action |
| --- | --- | --- |
| `exit=134`, `std::runtime_error`, every job of one era/channel | config points to a correction file that does not exist (fixed for UL2016 dimuon) | fix the config, resubmit the era/channel |
| `exit=54`, `analysis=True worker=False` | xrootd stage-out failure (`xrdcp`/`xrdfs`) | resubmit |
| `Missing required NanoAOD branch(es): run, ...` | all input files failed to open (xrootd auth/timeout) | resubmit |
| `RUNNING` but nothing in the queue | worker disconnected, job aborted (no termination record) | resubmit |
| status `OK` but `TNetXNGFile::Open` in `.err` | some input files of the list failed to open; TChain skipped them silently, so the job lost events | resubmit, even though it is OK |

## 4. Resubmit

Collect the bad jobs plus the OK jobs with input-open errors into one list
per era/channel, then resubmit. Always keep `--dy-estimation`.

```bash
for ERA in 2016PreVFP 2016PostVFP 2017 2018; do
  for CH in MuMu ElEl MuEl; do
    D=CondorJobs/condorLog/$STUDY/UL$ERA/$CH
    out=CondorJobs/rerun_${STUDY}_UL${ERA}_${CH}.txt
    : > $out
    for e in $(grep -l "TNetXNGFile::Open" $D/*/*.err 2>/dev/null); do
      s=$(basename $(dirname $e)); j=$(basename $e .err); j=${j##*_}
      echo "$s $j OPENERR $BASE/$ERA/$s/${s}_${j}.list" >> $out
    done
    f=CondorJobs/bad_jobs_${STUDY}_UL${ERA}_${CH}.txt
    [ -f $f ] && cat $f >> $out
    sort -u -k1,2 $out -o $out
    echo "UL$ERA $CH : $(wc -l < $out) jobs to rerun"
    [ -s "$out" ] && python3 CondorJobs/submit_jobs_v1.py --study $STUDY --run-period UL$ERA \
        --channel $CH --bad-jobs $out --dy-estimation
  done
done
```

Repeat sections 3 and 4 until every era/channel shows `Bad / unfinished : 0`
and `grep -l TNetXNGFile::Open` finds no `.err` file.

## 5. hadd

```bash
for ERA in 2016PreVFP 2016PostVFP 2017 2018; do
  SAMPLES=$(ls $BASE/$ERA | grep -E '^(Data_|DYJetsToLL_)')
  for CH in MuMu ElEl MuEl; do
    python3 CondorJobs/run_hadd_v8.py --study $STUDY --run-period UL$ERA --channel $CH \
        --input-base $BASE/$ERA --samples $SAMPLES          # add --check-only first
  done
done
```

Samples whose chunks are all `EMPTY_INPUT` (e.g. UL2016PreVFP
`Run2016B-HIPM`, ver1) have no output directory and are skipped.

Expected number of merged files (`<sample>/<sample>.root`), 160 in total:

| Era | MuMu | ElEl | MuEl |
| --- | --- | --- | --- |
| UL2016PreVFP | 13 | 13 | 18 |
| UL2016PostVFP | 10 | 10 | 13 |
| UL2017 | 14 | 14 | 19 |
| UL2018 | 12 | 8 | 16 |

```bash
OUT=/pnfs/knu.ac.kr/data/cms/store/user/$USER/CPV_Run2/ULSummer20/$STUDY
n=0; for d in $OUT/*/*/*/; do s=$(basename $d); [ -f $d/$s.root ] && n=$((n+1)); done; echo $n
```

## 6. Copy the merged files

```bash
# cluster
cd /pnfs/knu.ac.kr/data/cms/store/user/$USER/CPV_Run2/ULSummer20
for d in $STUDY/*/*/*/; do s=$(basename $d); [ -f $d$s.root ] && echo $d$s.root; done > ~/${STUDY}_merged_list.txt
zip -q ~/${STUDY}_merged.zip -@ < ~/${STUDY}_merged_list.txt

# local
scp <user>@cms02.knu.ac.kr:~/NanoAODv15_DYEst_v1_merged.zip .
unzip -q NanoAODv15_DYEst_v1_merged.zip
```

The zip keeps the layout `<study>/<era>/<channel>/<sample>/<sample>.root`.

## 7. R_out/in calculation

Set `basedir` in `CalDYSF_Run2UL_v15()` to the unzipped `<study>` directory,
then run it from a separate output directory:

```bash
root -l -b -q CalDYSF_Run2UL_v15.C 2>&1 | tee log_v15.txt
```

- Inputs: Data and amc@NLO DY (`DYJetsToLL_M_10To50`, `DYJetsToLL_M_50`);
  the madgraph DY files are produced but not used.
- Steps: 2 (dilepton), 3 (>= 2 jets), 5 (>= 1 b-tag).
- `h_cf_DiLepMass_<step>`: all events at step 2/3, MET > 40 at step 4/5 (METCut).
  `h_cf_NoMETCut_DiLepMass_<step>`: all events at step 2/3, MET <= 40 at step 4/5.
- SF used: METCut for ee and mumu, NoMETCut for emu.
- Outputs: `DYSF_Run2UL_NanoAODv15.txt` and `NanoAODv15_<era>_<channel>_nomet.pdf`.
  The keys in the txt contain spaces (e.g. `UL2018_>= 1 b-tag_MuMu`), same format
  as the earlier `DYSF_Run2UL.txt`.

## Known limitations

- UL2018 `DY_madgraph` has no `btagEff_UParTAK4.root`; all DY samples
  currently use the `DY` b-tag efficiency (`SSBCorrections::GetProcessSubDir`).
- `KinematicReconstruction` is constructed for every event in
  `Analysis::SetUpKINObs()` and its histograms and `TRandom3` are never freed.
- `main_ssb.cpp` checks `argc < 2` but reads `argv[1]` to `argv[7]`.
