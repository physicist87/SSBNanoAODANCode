#!/bin/bash

set -e

# ============================================================
# Configuration
# ============================================================

runPeriod="UL2018"
StudyName="Testv1"
Channels="MuMu"

# External XRootD input-list repository
InputListBase="/u/user/sha/Develop/CPviolation/SSB/AnalysisCode/NanoAODNtuple_v1/NtupleList_v1/2018_v6-FromGuks/FileList/InputList"

config="dimuon.config"
branch_list="${runPeriod}/branch_list_v15.txt"

# Interactive test: process only a few events
maxEvents=-1


# ============================================================
# Input samples
# Format:
#   SAMPLE_DIRECTORY/FILE_BASE
# ============================================================

inputlists=(
    #"TTbar_Signal/TTbar_Signal_1"
    #"/TTbar_Signal_1"
    "Data_SingleMuon_Run2018A/Data_SingleMuon_Run2018A_1"
)


echo "============================================================"
echo "[INFO] Interactive analysis test"
echo "============================================================"
echo "[INFO] runPeriod    = ${runPeriod}"
echo "[INFO] StudyName    = ${StudyName}"
echo "[INFO] Channel      = ${Channels}"
echo "[INFO] InputListBase= ${InputListBase}"
echo "[INFO] Config       = ${config}"
echo "[INFO] Branch list  = ${branch_list}"
echo "[INFO] Max events   = ${maxEvents}"
echo "============================================================"


# ============================================================
# Main loop
# ============================================================

for entry in "${inputlists[@]}"; do

    sample=$(dirname "${entry}")
    base=$(basename "${entry}")

    # --------------------------------------------------------
    # Input list
    #
    # Example:
    # /u/user/sha/.../InputList/TTbar_Signal/TTbar_Signal_1.list
    # --------------------------------------------------------

    inputlist="${InputListBase}/${sample}/${base}.list"

    if [ ! -f "${inputlist}" ]; then
        echo "[ERROR] Input list does not exist:"
        echo "        ${inputlist}"
        exit 1
    fi


    # --------------------------------------------------------
    # Output
    # --------------------------------------------------------

    relpath="${StudyName}/${runPeriod}/${Channels}/${sample}"
    outdir="output/${relpath}"
    outputfile="${relpath}/${base}.root"

    mkdir -p "${outdir}"

    echo
    echo "============================================================"
    echo "[INFO] Sample      : ${sample}"
    echo "[INFO] Input list  : ${inputlist}"
    echo "[INFO] Output      : output/${outputfile}"
    echo "============================================================"


    # --------------------------------------------------------
    # Run
    # --------------------------------------------------------

    echo "[RUN]"
    echo "./ssb_analysis \\"
    echo "  ${inputlist} \\"
    echo "  ${outputfile} \\"
    echo "  ULSummer20/${runPeriod}/${config} \\"
    echo "  None \\"
    echo "  ${runPeriod} \\"
    echo "  ${maxEvents} \\"
    echo "  ${branch_list}"
    echo

    ./ssb_analysis \
        "${inputlist}" \
        "${outputfile}" \
        "ULSummer20/${runPeriod}/${config}" \
        "None" \
        "${runPeriod}" \
        "${maxEvents}" \
        "${branch_list}"

done
