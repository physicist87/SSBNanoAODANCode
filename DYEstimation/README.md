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
| `DYEstimation/CalDYSF_Run2UL_v15.C` | R_out/in, k_ll and DY SF for all Run 2 UL eras (amc@NLO or madgraph DY) |

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

## 7. R_out/in calculation (`CalDYSF_Run2UL_v15.C`)

### Running

```bash
mkdir -p <output dir> && cd <output dir>
cp <repo>/DYEstimation/CalDYSF_Run2UL_v15.C .

# default: amc@NLO DY, basedir as written in the macro
root -l -b -q CalDYSF_Run2UL_v15.C 2>&1 | tee log_v15.txt

# madgraph DY
root -l -b -q 'CalDYSF_Run2UL_v15.C("madgraph")' 2>&1 | tee log_v15_madgraph.txt

# other location of the merged files
root -l -b -q 'CalDYSF_Run2UL_v15.C("amcatnlo", "/path/to/NanoAODv15_DYEst_v1")'
```

Arguments:

| Argument | Values | Default |
| --- | --- | --- |
| `dyGen` | `"amcatnlo"`, `"madgraph"` | `"amcatnlo"` |
| `basedir` | directory holding `<era>/<channel>/<sample>/<sample>.root` | the local unzipped study |

Before computing a period the macro checks that every Data and DY input file
exists. If one is missing it prints `[ERROR] Missing input: ...`, skips that
period and writes `# <period> skipped: ...` to the txt file, instead of
summing a partial set.

### DY MC generator option

| `dyGen` | DY samples | Outputs |
| --- | --- | --- |
| `amcatnlo` | `DYJetsToLL_M_10To50`, `DYJetsToLL_M_50` | `DYSF_Run2UL_NanoAODv15.txt`, `NanoAODv15_<era>_<ch>_nomet.pdf` |
| `madgraph` | `DYJetsToLL_M_10To50_madgraphMLM`, `DYJetsToLL_M_50_madgraphMLM` | `DYSF_Run2UL_NanoAODv15_madgraph.txt`, `NanoAODv15_madgraph_<era>_<ch>_nomet.pdf` |

- Both sets are produced by the Condor step (section 2 submits every
  `DYJetsToLL_*`), so switching needs no new jobs.
- UL2016PreVFP has no `DYJetsToLL_M_50_madgraphMLM` InputList, so with
  `madgraph` that period is skipped; the other periods run.
- The two generators give different R_out/in (shape outside the Z window)
  and N_out(MC); the difference can be used as a cross-check or a
  systematic on the DY normalization.
- Data inputs are identical for both options.

### Inputs and regions

Data: all `Data_*` samples of the channel (per era list in the macro;
UL2016PreVFP uses `Run2016B-HIPM-v2` and `C/D/E/F-HIPM`).

Steps (selection stage index in `Analysis`):

| Step | Selection |
| --- | --- |
| 2 | dilepton (no Z-veto in DY mode) |
| 3 | + >= 2 jets |
| 5 | + MET region of step 4 (METCut or NoMETCut) + >= 1 b-tag |

Histograms (m_ll, 1 GeV bins, 0-300 GeV):

| Histogram | Step 2, 3 | Step 4, 5 |
| --- | --- | --- |
| `h_cf_DiLepMass_<step>` (METCut) | all events | MET > 40 GeV (all channels incl. emu) |
| `h_cf_NoMETCut_DiLepMass_<step>` (NoMETCut) | all events | MET <= 40 GeV |

At steps 2 and 3 both histograms are identical, so METCut and NoMETCut give
the same SF there. The detailed table is printed for step 5 only.

`h_cf_DiLepMass_8` (stage 8, Top-Recon., after the b-tag step and the MET > 40
region) includes Z-window events in DY mode, for a DY SF at the kin. reco
step. Outputs produced before this change (the first `NanoAODv15_DYEst_v1`
production) have only Z-veto passing events at stage 8.

### Method

Z window: 76 <= m_ll < 106 GeV (bins `FindBin(76)` to `FindBin(105.9999)`).
For a histogram h:

```
N_in  = sum of h over the Z window
N_out = Integral(h) - N_in
```

Uncertainties are the quadrature sum of the bin errors in the same range.

For each step and each region (NoMETCut, METCut), with ll = ee or mumu:

```
k_ee  = sqrt( N_in^MC(ee)   / N_in^MC(mumu) )      # DY MC, NoMETCut, same step
k_mm  = sqrt( N_in^MC(mumu) / N_in^MC(ee)   )      # = 1 / k_ee

R_ll  = N_out^MC(ll) / N_in^MC(ll)                 # DY MC, same region

N_out^est(ll) = R_ll * ( N_in^Data(ll) - 0.5 * k_ll * N_in^Data(emu) )

SF_ll  = N_out^est(ll) / N_out^MC(ll)
SF_emu = sqrt( SF_ee * SF_mm )
```

- k_ll corrects the emu count for the ee/mumu efficiency difference; it is
  always taken from the NoMETCut DY MC of the same step.
- `0.5 * k_ll * N_in^Data(emu)` subtracts the non-DY (flavour-symmetric,
  mainly ttbar) contribution inside the Z window of the same-flavour channel.
- N_in^Data(emu) is the emu data in the same region, so it carries the
  MET > 40 cut in the METCut region.

Final DY yields:

```
N_DY(ll)  = N_out^MC(ll) * SF_ll        # SF_ll from METCut
N_DY(emu) = N^MC(emu)    * SF_emu       # SF_emu from NoMETCut, N^MC(emu) = full m_ll range
```

Error propagation (inputs treated as independent):

```
sigma(x1/x2)^2       = (y/x1)^2 sigma1^2 + (y/x2)^2 sigma2^2,          y = x1/x2
sigma(k)^2           = sigma1^2 / (4 x1 x2) + x1 sigma2^2 / (4 x2^3),  k = sqrt(x1/x2)
sigma(N_out^est)^2   = (N_in^ll - 0.5 k N_in^emu)^2 sigma(R)^2 + R^2 sigma(N_in^ll)^2
                     + (0.5 R k)^2 sigma(N_in^emu)^2 + (0.5 R N_in^emu)^2 sigma(k)^2
sigma(SF_emu)^2      = (SF_mm / SF_ee) sigma(SF_ee)^2 / 4 + (SF_ee / SF_mm) sigma(SF_mm)^2 / 4
sigma(N * SF)^2      = SF^2 sigma(N)^2 + N^2 sigma(SF)^2
```

### Outputs

- Console: per period the SF table for steps 2, 3, 5 and the detailed
  step 5 table (N_in, N_out, R_out/in, k_ll, N_in(Data), N_out(est), SF)
  for NoMETCut and METCut.
- txt: one line per period, step and channel:
  `key  SF  SF_err  N_MC  N_MC_err  N_DY  N_DY_err`.
  The key contains spaces (e.g. `UL2018_>= 1 b-tag_MuMu`), same format as
  the earlier `DYSF_Run2UL.txt`, so split from the right (last six fields
  are numbers).
- pdf: m_ll at step 5 NoMETCut, Data vs DY (M > 50 and 10-50) per period
  and channel.

## Known limitations

- UL2018 `DY_madgraph` has no `btagEff_UParTAK4.root`; all DY samples
  currently use the `DY` b-tag efficiency (`SSBCorrections::GetProcessSubDir`).
- `KinematicReconstruction` is constructed for every event in
  `Analysis::SetUpKINObs()` and its histograms and `TRandom3` are never freed.
- `main_ssb.cpp` checks `argc < 2` but reads `argv[1]` to `argv[7]`.
