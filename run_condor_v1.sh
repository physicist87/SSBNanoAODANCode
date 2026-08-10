#!/usr/bin/env bash

set -euo pipefail

# ============================================================
# Arguments from HTCondor JDL
# ============================================================

if [ "$#" -ne 7 ]; then
    echo "[ERROR] Wrong number of arguments"
    echo "Usage:"
    echo "  $0 <runPeriod> <StudyName> <Channel> <configFile> <branchList> <maxEvents> <inputList>"
    exit 1
fi

RUN_PERIOD="$1"
STUDY_NAME="$2"
CHANNEL="$3"
CONFIG_FILE="$4"
BRANCH_LIST="$5"
MAX_EVENTS="$6"
INPUT_LIST="$7"

WORK_DIR="${_CONDOR_SCRATCH_DIR:-$(pwd)}"

echo "============================================================"
echo " [WORKER START]"
echo "============================================================"
echo "Host        : $(hostname)"
echo "Work Dir    : ${WORK_DIR}"
echo "Run Period  : ${RUN_PERIOD}"
echo "Study Name  : ${STUDY_NAME}"
echo "Channel     : ${CHANNEL}"
echo "Config File : ${CONFIG_FILE}"
echo "Branch List : ${BRANCH_LIST}"
echo "Max Events  : ${MAX_EVENTS}"
echo "Input List  : ${INPUT_LIST}"
echo "============================================================"


# ============================================================
# XRootD I/O stability
# ============================================================

export XRD_REQUESTTIMEOUT=300
export XRD_STREAMTIMEOUT=300
export XRD_CONNECTIONRETRY=5
export XRD_WORKERTHREADS=4


# ============================================================
# CMSSW environment
# ============================================================

source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=el9_amd64_gcc12

CMSSW_VERSION="CMSSW_14_0_19"

echo "[INFO] Setting up ${CMSSW_VERSION}"

cd "${WORK_DIR}"

scram project CMSSW "${CMSSW_VERSION}"

cd "${CMSSW_VERSION}/src"
eval "$(scramv1 runtime -sh)"

cd "${WORK_DIR}"


# ============================================================
# Unpack analysis package
# ============================================================

PACKAGE_TARBALL="SSBNanoAODANCode.tar.gz"

if [ ! -f "${PACKAGE_TARBALL}" ]; then
    echo "[ERROR] Package tarball not found: ${PACKAGE_TARBALL}"
    exit 10
fi

tar -xzf "${PACKAGE_TARBALL}"

ANALYSIS_DIR="${WORK_DIR}/SSBNanoAODANCode"

if [ ! -d "${ANALYSIS_DIR}" ]; then
    echo "[ERROR] Analysis directory not found after unpacking:"
    echo "        ${ANALYSIS_DIR}"
    exit 11
fi

cd "${ANALYSIS_DIR}"


# ============================================================
# Check runtime inputs
# ============================================================

if [ ! -f "input/${INPUT_LIST}" ]; then
    echo "[ERROR] Input list not found:"
    echo "        input/${INPUT_LIST}"
    exit 20
fi

if [ ! -f "configs/ULSummer20/${RUN_PERIOD}/${CONFIG_FILE}" ]; then
    echo "[ERROR] Config file not found:"
    echo "        configs/ULSummer20/${RUN_PERIOD}/${CONFIG_FILE}"
    exit 21
fi

if [ ! -f "branchlist/${BRANCH_LIST}" ]; then
    echo "[ERROR] Branch list not found:"
    echo "        branchlist/${BRANCH_LIST}"
    exit 22
fi


# ============================================================
# Compile analysis
# ============================================================

echo "[INFO] Compiling ssb_analysis"

make -f Makefile_ssb clean
make -f Makefile_ssb

if [ ! -x "./ssb_analysis" ]; then
    echo "[ERROR] ssb_analysis compilation failed"
    exit 30
fi


# ============================================================
# Prepare local output
# ============================================================

SAMPLE_DIR=$(dirname "${INPUT_LIST}")
INPUT_BASE=$(basename "${INPUT_LIST}" .list)

RELPATH="${STUDY_NAME}/${RUN_PERIOD}/${CHANNEL}/${SAMPLE_DIR}"
OUTPUT_FILE="${RELPATH}/${INPUT_BASE}.root"

mkdir -p "output/${RELPATH}"

echo "[INFO] Local output:"
echo "       output/${OUTPUT_FILE}"


# ============================================================
# Run analysis
# ============================================================

echo "============================================================"
echo "[RUN]"
echo "./ssb_analysis \\"
echo "  ${INPUT_LIST} \\"
echo "  ${OUTPUT_FILE} \\"
echo "  ULSummer20/${RUN_PERIOD}/${CONFIG_FILE} \\"
echo "  None \\"
echo "  ${RUN_PERIOD} \\"
echo "  ${MAX_EVENTS} \\"
echo "  ${BRANCH_LIST}"
echo "============================================================"

./ssb_analysis \
    "${INPUT_LIST}" \
    "${OUTPUT_FILE}" \
    "ULSummer20/${RUN_PERIOD}/${CONFIG_FILE}" \
    "None" \
    "${RUN_PERIOD}" \
    "${MAX_EVENTS}" \
    "${BRANCH_LIST}"


# ============================================================
# Validate local output
# ============================================================

LOCAL_OUTPUT="output/${OUTPUT_FILE}"

if [ ! -s "${LOCAL_OUTPUT}" ]; then
    echo "[ERROR] Output ROOT file was not created or is empty:"
    echo "        ${LOCAL_OUTPUT}"
    exit 40
fi

echo "[INFO] Output created successfully:"
ls -lh "${LOCAL_OUTPUT}"

echo "============================================================"
echo " [WORKER DONE]"
echo "============================================================"
