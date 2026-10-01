#!/bin/bash

set -e

# ============================================================
# Configuration
# ============================================================

runPeriod="UL2017"
runPeriod="UL2018"
StudyName="Testv1"
Channels="MuMu"

# External XRootD input-list repository
InputListBase="/u/user/sha/Develop/CPviolation/SSB/AnalysisCode/NanoAODNtuple_v1/NtupleList_v1/NanoAODv15_v2/InputList/2018"

# Interactive test: process only a few events
maxEvents=-1


# ============================================================
# Config / branch-list selection
#
# UL2017's HLT menu changed after Run2017B (e.g.
# HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8 only exists from RunC
# onward, and RunB's config also requires HLT_IsoMu24_eta2p1 in addition
# to HLT_IsoMu27) - so a single flat config/branch-list pair is wrong for
# UL2017 data. Mirrors get_job_config() in CondorJobs/submit_jobs_v1.py.
# ============================================================

get_job_config() {
    local run_period="$1"
    local channel="$2"
    local sample="$3"

    local config_map_MuMu="dimuon.config"
    local config_map_ElEl="dielec.config"
    local config_map_MuEl="muelec.config"

    local base_prefix
    case "${channel}" in
        MuMu) base_prefix="dimuon" ;;
        ElEl) base_prefix="dielec" ;;
        MuEl) base_prefix="muelec" ;;
        *)    echo "[ERROR] Unknown channel: ${channel}" >&2; exit 1 ;;
    esac

    local config_file="${base_prefix}.config"
    local branch_list="${run_period}/branch_list_v15.txt"

    if [ "${run_period}" = "UL2017" ] && [[ "${sample}" == Data_* ]]; then
        if [[ "${sample}" == *Run2017B* ]]; then
            config_file="${base_prefix}_Data_RunB.config"
            branch_list="${run_period}/branch_list_Run2017B_v15.txt"
        else
            config_file="${base_prefix}_Data_RunCtoF.config"
            branch_list="${run_period}/branch_list_Run2017CtoF_v15.txt"
        fi
    fi

    echo "${config_file} ${branch_list}"
}


# ============================================================
# Input samples
# Format:
#   SAMPLE_DIRECTORY/FILE_BASE
# ============================================================

inputlists=(
    "TTbar_Signal/TTbar_Signal_1"
    #"/TTbar_Signal_1"
    #"Data_SingleMuon_Run2018A/Data_SingleMuon_Run2018A_1"
    #"Data_SingleMuon_Run2017B/Data_SingleMuon_Run2017B_1"
    #"TTbar_DiLepBKG/TTbar_DiLepBKG_1"
    #"ST_s-channel_4f_leptonDecays/ST_s-channel_4f_leptonDecays_1"
)


echo "============================================================"
echo "[INFO] Interactive analysis test"
echo "============================================================"
echo "[INFO] runPeriod    = ${runPeriod}"
echo "[INFO] StudyName    = ${StudyName}"
echo "[INFO] Channel      = ${Channels}"
echo "[INFO] InputListBase= ${InputListBase}"
echo "[INFO] Max events   = ${maxEvents}"
echo "============================================================"


# ============================================================
# Main loop
# ============================================================

for entry in "${inputlists[@]}"; do

    sample=$(dirname "${entry}")
    base=$(basename "${entry}")

    # --------------------------------------------------------
    # Config / branch list (per-sample, era-aware for UL2017 data)
    # --------------------------------------------------------

    read -r config branch_list <<< "$(get_job_config "${runPeriod}" "${Channels}" "${sample}")"

    echo "[INFO] Sample      : ${sample}"
    echo "[INFO] Config       = ${config}"
    echo "[INFO] Branch list  = ${branch_list}"

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
