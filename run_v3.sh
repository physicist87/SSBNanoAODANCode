#!/bin/bash

# === [CONFIGURATION] ===
runPeriod="UL2017"
runPeriod="UL2016PreVFP"
StudyName="Testv1"
runPeriod="UL2018"
Channels="MuMu"

echo "[INFO] runPeriod = $runPeriod"
echo "[INFO] Channels  = $Channels"

# List of input samples in the form: "Directory/FileBase"
inputlists=(
    "TTbar_Signal_1"
    #"Data_SingleMuon_Run2017C/Data_SingleMuon_Run2017C_1"
    #"Data_SingleMuon_Run2017E/Data_SingleMuon_Run2017E_1"
    #"Data_SingleMuon_Run2017F/Data_SingleMuon_Run2017F_1"
    #"Data_SingleMuon_Run2018B/Data_SingleMuon_Run2018B_1"
    #"Data_SingleMuon_Run2018C/Data_SingleMuon_Run2018C_1"
    #"Data_DoubleMuon_Run2018C/Data_DoubleMuon_Run2018C_1"
    #"Data_DoubleMuon_Run2018A/Data_DoubleMuon_Run2018A_1"
    #"Data_SingleMuon_Run2018A/Data_SingleMuon_Run2018A_1"
    #"Data_SingleMuon_Run2018A/Data_SingleMuon_Run2018A_1"
    #"TTbar_Signal/TTbar_Signal_1"
)
config="dimuon_Data_RunCtoF.config"
config="dimuon.config"

branch_list="${runPeriod}/branch_list_Run2017CtoF.txt"
branch_list="${runPeriod}/branch_list.txt"
branch_list="${runPeriod}/branch_list_v15.txt"


# === [MAIN LOOP OVER INPUT LISTS] ===
for entry in "${inputlists[@]}"; do
    dir=$(dirname "$entry")     # e.g. Data_SingleMuon_Run2017C
    base=$(basename "$entry")   # e.g. Data_SingleMuon_Run2017C_1

    # Build relative path and full output directory path
    relpath="${StudyName}/${runPeriod}/${Channels}/${dir}"
    outdir="output/${relpath}"
    mkdir -p "$outdir"

    echo "[INFO] Output directory will be: $outdir"

    # Build input list and output file path
    inputlist="${runPeriod}/${dir}/${base}.list"
    outputfile="${relpath}/${base}.root"

    #echo "[RUN] ./ssb_analysis $inputlist $outputfile $config None $runPeriod -1 $branch_list"
    echo "[RUN] ./ssb_analysis $inputlist $outputfile ULSummer20/$runPeriod/$config None $runPeriod -1 $branch_list"
    ./ssb_analysis "$inputlist" "$outputfile" "ULSummer20/$runPeriod/$config" "None" "$runPeriod" 10 "$branch_list"
done
