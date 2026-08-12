#!/usr/bin/env bash

set -euo pipefail


# ============================================================
# Usage
#
# run_condor_v1.sh \
#   <runPeriod> \
#   <studyName> \
#   <channel> \
#   <sample> \
#   <configFile> \
#   <branchList> \
#   <maxEvents> \
#   <inputList>
#
# Example:
#
# run_condor_v1.sh \
#   UL2018 \
#   Testv1 \
#   MuMu \
#   TTbar_Signal \
#   dimuon.config \
#   UL2018/branch_list_v15.txt \
#   -1 \
#   TTbar_Signal_1.list
# ============================================================


if [ "$#" -ne 8 ]; then
    echo "[ERROR] Wrong number of arguments: $#"
    echo
    echo "Usage:"
    echo "  $0 <runPeriod> <studyName> <channel> <sample> <configFile> <branchList> <maxEvents> <inputList>"
    exit 1
fi


# ============================================================
# Arguments from Condor JDL
# ============================================================

RUN_PERIOD="$1"
STUDY_NAME="$2"
CHANNEL="$3"
SAMPLE="$4"
CONFIG_FILE="$5"
BRANCH_LIST="$6"
MAX_EVENTS="$7"
INPUT_LIST="$8"


# ============================================================
# Basic configuration
# ============================================================

WORK_DIR="${_CONDOR_SCRATCH_DIR:-$(pwd)}"

PACKAGE_TARBALL="SSBNanoAODANCode.tar.gz"

SCRAM_ARCH_VERSION="el9_amd64_gcc12"
CMSSW_VERSION="CMSSW_14_0_19"

# User-dependent SE path
USER_ID="$(whoami)"

SE_HOST="root://cluster142.knu.ac.kr"
SE_BASE="/store/user/${USER_ID}/CPV_Run2/ULSummer20"


echo "============================================================"
echo "=== [WORKER START] =========================================="
echo "============================================================"
echo "[INFO] Host          : $(hostname)"
echo "[INFO] Date          : $(date)"
echo "[INFO] User ID       : ${USER_ID}"
echo "[INFO] User info     : $(id)"
echo "[INFO] Work Dir      : ${WORK_DIR}"
echo
echo "[INFO] Run Period    : ${RUN_PERIOD}"
echo "[INFO] Study Name    : ${STUDY_NAME}"
echo "[INFO] Channel       : ${CHANNEL}"
echo "[INFO] Sample        : ${SAMPLE}"
echo "[INFO] Config File   : ${CONFIG_FILE}"
echo "[INFO] Branch List   : ${BRANCH_LIST}"
echo "[INFO] Max Events    : ${MAX_EVENTS}"
echo "[INFO] Input List    : ${INPUT_LIST}"
echo
echo "[INFO] SE Host       : ${SE_HOST}"
echo "[INFO] SE Base       : ${SE_BASE}"
echo "============================================================"

cd "${WORK_DIR}"


# ============================================================
# 1. XRootD I/O settings
# ============================================================

export XRD_REQUESTTIMEOUT=300
export XRD_STREAMTIMEOUT=300
export XRD_CONNECTIONRETRY=5
export XRD_WORKERTHREADS=4

echo
echo "[INFO] XRootD settings configured."


# ============================================================
# 2. Check Condor-transferred files
# ============================================================

if [ ! -f "${PACKAGE_TARBALL}" ]; then
    echo "[ERROR] Analysis package tarball not found:"
    echo "        ${WORK_DIR}/${PACKAGE_TARBALL}"
    exit 10
fi

if [ ! -f "${INPUT_LIST}" ]; then
    echo "[ERROR] Input list not found:"
    echo "        ${WORK_DIR}/${INPUT_LIST}"
    exit 11
fi

if [ ! -s "${INPUT_LIST}" ]; then
    echo "[ERROR] Input list is empty:"
    echo "        ${WORK_DIR}/${INPUT_LIST}"
    exit 12
fi

echo
echo "[INFO] Files transferred to worker:"
ls -lh


# ============================================================
# 3. CMS environment
# ============================================================

echo
echo "============================================================"
echo "=== [CMSSW SETUP] ==========================================="
echo "============================================================"

source /cvmfs/cms.cern.ch/cmsset_default.sh

export SCRAM_ARCH="${SCRAM_ARCH_VERSION}"

echo "[INFO] SCRAM_ARCH    : ${SCRAM_ARCH}"
echo "[INFO] CMSSW Version : ${CMSSW_VERSION}"


# ============================================================
# 4. Create local CMSSW project
# ============================================================

cd "${WORK_DIR}"

if [ ! -d "${CMSSW_VERSION}" ]; then
    echo "[INFO] Creating local CMSSW project..."
    scram project CMSSW "${CMSSW_VERSION}"
else
    echo "[INFO] Existing ${CMSSW_VERSION} directory found."
fi

cd "${WORK_DIR}/${CMSSW_VERSION}/src"

eval "$(scram runtime -sh)"

echo "[INFO] CMSSW_BASE    : ${CMSSW_BASE}"
echo "[INFO] ROOT version  : $(root-config --version)"

cd "${WORK_DIR}"


# ============================================================
# 5. Unpack analysis package
# ============================================================

echo
echo "============================================================"
echo "=== [UNPACK PACKAGE] ========================================"
echo "============================================================"

tar -xzf "${PACKAGE_TARBALL}"

ANALYSIS_DIR="${WORK_DIR}/SSBNanoAODANCode"

if [ ! -d "${ANALYSIS_DIR}" ]; then
    echo "[ERROR] Analysis directory not found after unpacking:"
    echo "        ${ANALYSIS_DIR}"
    echo
    echo "[INFO] Current directory contents:"
    ls -lah "${WORK_DIR}"
    exit 20
fi

cd "${ANALYSIS_DIR}"

echo "[INFO] Analysis directory:"
pwd


# ============================================================
# 6. Check required runtime files
# ============================================================

#
# IMPORTANT:
#
# ssb_analysis internally prepends:
#
#   ./configs/
#
# to the config argument.
#
# Therefore:
#
# CONFIG_ARG:
#   ULSummer20/UL2018/dimuon.config
#
# Actual local file:
#   configs/ULSummer20/UL2018/dimuon.config
#
# Same idea for branchlist.
#

CONFIG_ARG="ULSummer20/${RUN_PERIOD}/${CONFIG_FILE}"
BRANCH_ARG="${BRANCH_LIST}"

CONFIG_LOCAL="configs/${CONFIG_ARG}"
BRANCH_LOCAL="branchlist/${BRANCH_ARG}"


if [ ! -f "Makefile_ssb" ]; then
    echo "[ERROR] Makefile_ssb not found."
    exit 21
fi

if [ ! -f "${CONFIG_LOCAL}" ]; then
    echo "[ERROR] Config file not found:"
    echo "        ${ANALYSIS_DIR}/${CONFIG_LOCAL}"
    exit 22
fi

if [ ! -f "${BRANCH_LOCAL}" ]; then
    echo "[ERROR] Branch list not found:"
    echo "        ${ANALYSIS_DIR}/${BRANCH_LOCAL}"
    exit 23
fi


# ============================================================
# 7. Compile analysis on worker
# ============================================================

echo
echo "============================================================"
echo "=== [COMPILE] ==============================================="
echo "============================================================"

make -f Makefile_ssb clean
make -f Makefile_ssb

if [ ! -x "./ssb_analysis" ]; then
    echo "[ERROR] ssb_analysis was not created."
    exit 30
fi

echo
echo "[INFO] Compilation successful:"
ls -lh ./ssb_analysis


# ============================================================
# 8. Prepare input/output paths
# ============================================================

INPUT_LIST_PATH="${WORK_DIR}/${INPUT_LIST}"

INPUT_BASE="$(basename "${INPUT_LIST}" .list)"

#
# With SE directory = "None",
# ssb_analysis writes to:
#
#   output/<argv[2]>
#
OUTPUT_RELATIVE="${STUDY_NAME}/${RUN_PERIOD}/${CHANNEL}/${SAMPLE}/${INPUT_BASE}.root"

LOCAL_OUTPUT="${ANALYSIS_DIR}/output/${OUTPUT_RELATIVE}"

mkdir -p "$(dirname "${LOCAL_OUTPUT}")"


echo
echo "============================================================"
echo "=== [ANALYSIS CONFIGURATION] ================================"
echo "============================================================"
echo "[INFO] Input List    : ${INPUT_LIST_PATH}"
echo "[INFO] Local Output  : ${LOCAL_OUTPUT}"
echo "[INFO] Config Arg    : ${CONFIG_ARG}"
echo "[INFO] Config Local  : ${CONFIG_LOCAL}"
echo "[INFO] Branch Arg    : ${BRANCH_ARG}"
echo "[INFO] Branch Local  : ${BRANCH_LOCAL}"
echo "============================================================"


# ============================================================
# 9. Show XRootD input files
# ============================================================

echo
echo "[INFO] Input ROOT file(s):"
echo "------------------------------------------------------------"
cat "${INPUT_LIST_PATH}"
echo "------------------------------------------------------------"


# ============================================================
# 10. Run analysis
# ============================================================

echo
echo "============================================================"
echo "=== [RUN ANALYSIS] =========================================="
echo "============================================================"

echo "./ssb_analysis \\"
echo "  ${INPUT_LIST_PATH} \\"
echo "  ${OUTPUT_RELATIVE} \\"
echo "  ${CONFIG_ARG} \\"
echo "  None \\"
echo "  ${RUN_PERIOD} \\"
echo "  ${MAX_EVENTS} \\"
echo "  ${BRANCH_ARG}"
echo "============================================================"
echo


./ssb_analysis \
    "${INPUT_LIST_PATH}" \
    "${OUTPUT_RELATIVE}" \
    "${CONFIG_ARG}" \
    "None" \
    "${RUN_PERIOD}" \
    "${MAX_EVENTS}" \
    "${BRANCH_ARG}"


# ============================================================
# 11. Validate local output
# ============================================================

echo
echo "============================================================"
echo "=== [LOCAL OUTPUT CHECK] ===================================="
echo "============================================================"

if [ ! -f "${LOCAL_OUTPUT}" ]; then
    echo "[ERROR] Output ROOT file was not created:"
    echo "        ${LOCAL_OUTPUT}"
    exit 40
fi

if [ ! -s "${LOCAL_OUTPUT}" ]; then
    echo "[ERROR] Output ROOT file exists but is empty:"
    echo "        ${LOCAL_OUTPUT}"
    exit 41
fi

echo "[INFO] Local ROOT output created successfully:"
ls -lh "${LOCAL_OUTPUT}"


# ============================================================
# 12. Prepare SE destination
#
# NFS view:
#
# /pnfs/knu.ac.kr/data/cms/store/user/${USER_ID}/CPV_Run2/ULSummer20/
#
# XRootD view:
#
# root://cluster142.knu.ac.kr//store/user/${USER_ID}/CPV_Run2/ULSummer20/
# ============================================================

SE_DEST_DIR="${SE_BASE}/${STUDY_NAME}/${RUN_PERIOD}/${CHANNEL}/${SAMPLE}"

SE_DEST_PATH="${SE_DEST_DIR}/${INPUT_BASE}.root"

SE_DEST_URL="${SE_HOST}/${SE_DEST_PATH}"


echo
echo "============================================================"
echo "=== [STAGE-OUT PREPARATION] ================================"
echo "============================================================"
echo "[INFO] Local file    : ${LOCAL_OUTPUT}"
echo "[INFO] SE directory  : ${SE_DEST_DIR}"
echo "[INFO] SE URL        : ${SE_DEST_URL}"
echo "============================================================"


# ============================================================
# 13. Create destination directory on SE
# ============================================================

echo
echo "[INFO] Creating SE destination directory..."

xrdfs "${SE_HOST}" mkdir -p "${SE_DEST_DIR}"


# ============================================================
# 14. Stage-out with XRootD
# ============================================================

echo
echo "============================================================"
echo "=== [STAGE-OUT] ============================================="
echo "============================================================"

xrdcp -f "${LOCAL_OUTPUT}" "${SE_DEST_URL}"


# ============================================================
# 15. Verify staged-out file
# ============================================================

echo
echo "============================================================"
echo "=== [STAGE-OUT CHECK] ======================================="
echo "============================================================"

xrdfs "${SE_HOST}" stat "${SE_DEST_PATH}"


# ============================================================
# Done
# ============================================================

echo
echo "============================================================"
echo "=== [WORKER DONE] ==========================================="
echo "============================================================"
echo "[INFO] Host       : $(hostname)"
echo "[INFO] Date       : $(date)"
echo "[INFO] User       : ${USER_ID}"
echo "[INFO] Local      : ${LOCAL_OUTPUT}"
echo "[INFO] SE output  : ${SE_DEST_URL}"
echo "============================================================"

exit 0
