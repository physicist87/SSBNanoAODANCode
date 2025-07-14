import os
import re
import subprocess
import glob

def isDirectoryInDataList(directory, channel, runPeriod, debug=False):
    dimuonList = ["SingleMuon", "DoubleMuon"]
    dielecList = ["SingleElectron", "DoubleEG"]
    muelecList = ["SingleMuon", "SingleElectron", "MuonEG"]
    dataList = []

    if runPeriod == "UL2018":
        dielecList = ["EGamma"]
        muelecList = ["SingleMuon", "EGamma", "MuonEG"]

    if channel == "MuMu":
        dataList = dimuonList
    elif channel == "ElEl":
        dataList = dielecList
    elif channel == "MuEl":
        dataList = muelecList

    if debug:
        print(f"Data list for {channel}: {dataList}")

    for item in dataList:
        if item in directory:
            return True
    return False

def contains_str_prefix(disc, string):
    return disc in string

def has_number_suffix(filename):
    return bool(re.match(r".+_\d+$", filename))

def needs_resubmission(outputPath, sampleDir, listFileName):
    if not has_number_suffix(listFileName):
        return False
    expected_root_file = os.path.join(outputPath, sampleDir, f"{listFileName}.root")
    if not os.path.isfile(expected_root_file):
        return True
    elif os.path.getsize(expected_root_file) < 5120:
        return True
    return False

def submit_jobs(sampleType, inputListPath, runScriptPath, logDir, studyName, runPeriod, channel, outputPath, submitDir, MainPath, samples, configFile, branchList, maxEvents="-1", debug=False, resubmit=False):
    submitDir = os.path.join(submitDir, studyName, runPeriod, channel)
    if debug:
        print(f"Creating submission directory: {submitDir}")
    os.makedirs(submitDir, exist_ok=True)

    if "all" in samples:
        samples = [name for name in os.listdir(inputListPath) if os.path.isdir(os.path.join(inputListPath, name))]
        if debug:
            print(f"Detected 'all' in samples. Using all directories in {inputListPath}: {samples}")

    for sampleDir in samples:
        if debug:
            print(f"Processing sample directory: {sampleDir}")
        
        if contains_str_prefix("Data_", sampleDir):
            if debug:
                print("Detected data sample")
            if not isDirectoryInDataList(sampleDir, channel, runPeriod, debug):
                if debug:
                    print(f"Skipping data sample directory (not in data list): {sampleDir}")
                continue
        else:
            if debug:
                print("Detected MC sample")

        sampleInputPath = os.path.join(inputListPath, sampleDir)
        sampleOutputPath = os.path.join(outputPath, sampleDir)
        jobLogDir = os.path.join(logDir, studyName, runPeriod, channel, sampleDir)
        os.makedirs(jobLogDir, exist_ok=True)

        file_prefix = "resubmit_" if resubmit else ""

        submit_file_content = f"""Universe   = vanilla
Executable = {runScriptPath}
Log        = {jobLogDir}/{file_prefix}{sampleDir}_$(InputListName).log
Output     = {jobLogDir}/{file_prefix}{sampleDir}_$(InputListName).out
Error      = {jobLogDir}/{file_prefix}{sampleDir}_$(InputListName).err
RequestMemory = 1 GB
RequestCpus = 1
should_transfer_files = YES
when_to_transfer_output = ON_EXIT
transfer_input_files = {runScriptPath}
Arguments  = {MainPath} {runPeriod} {studyName} {channel} {configFile} {branchList} {maxEvents} {outputPath} {sampleDir}/{sampleDir}_$(InputListName)
Queue InputListName from (
"""

        if any(x in jobLogDir for x in ["TTbar_Signal","TTZToLLNuNu","TTWJetsToLNu"]):
            submit_file_content = f"""Universe   = vanilla
Executable = {runScriptPath}
Log        = {jobLogDir}/{file_prefix}{sampleDir}_$(InputListName).log
Output     = {jobLogDir}/{file_prefix}{sampleDir}_$(InputListName).out
Error      = {jobLogDir}/{file_prefix}{sampleDir}_$(InputListName).err
RequestMemory = 10 GB
RequestCpus = 1
# 필요한 환경 파일 전송
should_transfer_files = YES
when_to_transfer_output = ON_EXIT
transfer_input_files = {runScriptPath}
# CMSSW 환경 설정 스크립트를 통해 실행
Arguments  = {MainPath} {runPeriod} {studyName} {channel} {configFile} {branchList} {maxEvents} {outputPath} {sampleDir}/{sampleDir}_$(InputListName)
Queue InputListName from (
"""


        resubmit_list = []
        for listFile in os.listdir(sampleInputPath):
            if listFile.endswith(".list"):
                listFileName = os.path.splitext(listFile)[0]
                if listFileName == sampleDir or not listFileName.startswith(f"{sampleDir}_"):
                    if debug:
                        print(f"Skipping file that doesn't match pattern: {listFileName}")
                    continue
                number_part = listFileName[len(sampleDir)+1:]
                if not number_part.isdigit():
                    if debug:
                        print(f"Skipping file without numeric suffix: {listFileName}")
                    continue

                if resubmit and needs_resubmission(outputPath, sampleDir, listFileName):
                    resubmit_list.append(number_part)
                elif not resubmit:
                    submit_file_content += f"{number_part}\n"

        if resubmit and resubmit_list:
            for file_name in resubmit_list:
                submit_file_content += f"{file_name}\n"
            submit_file_content += ")"
            submit_file_path = os.path.join(submitDir, f"submit_{sampleDir}.sub")
            with open(submit_file_path, "w") as submit_file:
                submit_file.write(submit_file_content)
            if debug:
                print(f"Submitting Condor job for {sampleDir} with resubmission list: {resubmit_list}")
            subprocess.run(["condor_submit", "-append", "accounting_group = group_cms", submit_file_path])
        elif not resubmit:
            submit_file_content += ")"
            submit_file_path = os.path.join(submitDir, f"submit_{sampleDir}.sub")
            with open(submit_file_path, "w") as submit_file:
                submit_file.write(submit_file_content)
            if debug:
                print(f"Submitting Condor job for {sampleDir}")
            subprocess.run(["condor_submit", "-append", "accounting_group = group_cms", submit_file_path])

def main():
    # ===== Main Configuration =====
    MainPath = os.path.abspath("..")
    inputList = os.path.join(MainPath, "input/")
    logDir = os.path.join(os.getcwd(), "condorLog")
    submitDir = os.path.join(os.getcwd(), "condorSubmit")
    runScriptPath = os.path.join(MainPath, "run_cmd_v7.sh")
    outputPath = os.path.join(MainPath,"output/")
    runPeriod = "UL2018"
    studyName = "TestRun"
    channel = "MuMu"
    maxEvents = "-1"
    #maxEvents = "100"

    #make the list of sample name in input directory 
    #path = os.path.join(inputList,runPeriod)
    #dir_paths = glob.glob(os.path.join(path, '*/'))
    #samples = [os.path.basename(os.path.normpath(d)) for d in dir_paths]

    samples = "all"
    
    branchList = "UL2018/branch_list.txt"
    configFile = "dimuon.config"
    debug = True

    submit_jobs(
        "None", 
        f"{inputList}/{runPeriod}", 
        runScriptPath, 
        logDir, 
        studyName, 
        runPeriod, 
        channel, 
        outputPath, 
        submitDir, 
        MainPath, 
        samples, 
        configFile, 
        branchList, 
        maxEvents=maxEvents, 
        debug=debug, 
        resubmit=False
    )

if __name__ == "__main__":
    main()
