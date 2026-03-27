import os
import sys
import subprocess
import re
import time
import argparse
from collections import defaultdict

def isDirectoryInDataList(directory, channel, run_period):
    """
    Filter data samples based on the analysis channel and year.
    Matches the logic used in year-specific scripts.
    """
    dimuon_list = ["SingleMuon", "DoubleMuon"]
    dielec_list = ["SingleElectron", "DoubleEG"]
    muelec_list = ["SingleMuon", "SingleElectron", "MuonEG"]
    
    # Special handling for UL2018: Electron/EG samples were merged into EGamma
    if run_period == "UL2018":
        dielec_list = ["EGamma"]
        muelec_list = ["SingleMuon", "EGamma", "MuonEG"]

    if channel == "MuMu":
        data_list = dimuon_list
    elif channel == "ElEl":
        data_list = dielec_list
    elif channel == "MuEl":
        data_list = muelec_list
    else:
        return False

    # Check if any keyword matches the directory name
    return any(item in directory for item in data_list)

def submit_condor_job(submit_file_path):
    """
    Execute condor_submit command.
    """
    cmd = ["condor_submit", submit_file_path]
    print(f"[INFO] Submitting: {' '.join(cmd)}")
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] condor_submit failed: {e}")

def get_job_config(year, channel, sample_name):
    """
    Determine the correct .config and branch_list.txt based on year and sample type.
    Implemented based on logic from all 4 year-specific scripts.
    """
    base_config_map = {"MuMu": "dimuon", "ElEl": "dielec", "MuEl": "muelec"}
    base_conf = base_config_map[channel]
    
    # Default settings for MC and standard data
    config_file = f"{base_conf}.config"
    branch_list = f"{year}/branch_list.txt"

    # Mapping logic for UL2017 Data
    if year == "UL2017" and sample_name.startswith("Data_"):
        if "Run2017B" in sample_name:
            config_file = f"{base_conf}_Data_RunB.config"
            branch_list = f"{year}/branch_list_Run2017B.txt"
        else:
            config_file = f"{base_conf}_Data_RunCtoF.config"
            branch_list = f"{year}/branch_list_Run2017CtoF.txt"
            
    # Mapping logic for UL2016PostVFP Data (specifically RunH)
    elif year == "UL2016PostVFP" and sample_name.startswith("Data_"):
        if "Run2016H" in sample_name:
            config_file = f"{base_conf}_Data_RunH.config"
            branch_list = f"{year}/branch_list_RunH.txt"
            
    return config_file, branch_list

def check_job_status(log_dir, study_name, run_period, channel, sample_dir, list_num):
    """
    Verify if a job is completed by checking the end of the output log.
    Used for resubmission logic and check reporting.
    """
    job_log_dir = os.path.join(log_dir, study_name, run_period, channel, sample_dir)
    for prefix in ["resubmit_", ""]:
        out_file = os.path.join(job_log_dir, f"{prefix}{sample_dir}_{list_num}.out")
        if os.path.exists(out_file):
            try:
                with open(out_file, 'r') as f:
                    content = f.read()
                if "Analysis destructor completed." in content:
                    return 'completed'
            except:
                pass
    return 'not_completed'

def run_process(args):
    """
    Main workflow for generating JDLs and submitting HTCondor jobs.
    """
    this_script_dir = os.path.dirname(os.path.abspath(__file__))
    main_path = os.path.abspath(os.path.join(this_script_dir, ".."))
    
    input_list_path = os.path.join(main_path, f"input/{args.year}")
    log_dir = os.path.join(os.getcwd(), "condorLog")
    submit_dir = os.path.join(os.getcwd(), "condorSubmit", args.study, args.year, args.channel)
    run_script_path = os.path.join(main_path, "run_cmd_v7.sh")
    output_path = "/pnfs/knu.ac.kr/data/cms/store/user/sha/CPV_Run2/ULSummer20"
    
    os.makedirs(submit_dir, exist_ok=True)

    if not os.path.exists(input_list_path):
        print(f"[ERROR] Input directory not found: {input_list_path}")
        return

    all_samples = sorted([d for d in os.listdir(input_list_path) if os.path.isdir(os.path.join(input_list_path, d))])
    
    for sample in all_samples:
        # Filter Data samples by channel and year
        if sample.startswith("Data_"):
            if not isDirectoryInDataList(sample, args.channel, args.year):
                continue

        # Get specific config mapping
        config_file, branch_list = get_job_config(args.year, args.channel, sample)
        sample_input_path = os.path.join(input_list_path, sample)
        job_log_dir = os.path.join(log_dir, args.study, args.year, args.channel, sample)
        os.makedirs(job_log_dir, exist_ok=True)

        # Collect list file ID numbers
        list_nums = []
        for f in sorted(os.listdir(sample_input_path)):
            if f.endswith(".list"):
                match = re.search(r'_(\d+)\.list$', f)
                if match:
                    list_nums.append(match.group(1))

        # Filter jobs for resubmission
        if args.resubmit:
            to_submit = [n for n in list_nums if check_job_status(log_dir, args.study, args.year, args.channel, sample, n) != 'completed']
        else:
            to_submit = list_nums

        if not to_submit:
            continue

        # Generate JDL (.sub) file using robust syntax to avoid parsing errors
        prefix = "resubmit_" if args.resubmit else ""
        jdl_path = os.path.join(submit_dir, f"{prefix}{sample}.sub")
        
        with open(jdl_path, "w") as f:
            f.write(f"Universe              = vanilla\n")
            f.write(f"Executable            = {run_script_path}\n")
            f.write(f"Log                   = {job_log_dir}/{prefix}{sample}_$(InputListName).log\n")
            f.write(f"Output                = {job_log_dir}/{prefix}{sample}_$(InputListName).out\n")
            f.write(f"Error                 = {job_log_dir}/{prefix}{sample}_$(InputListName).err\n")
            f.write(f"RequestMemory         = 2 GB\n")
            f.write(f"RequestCpus           = 1\n")
            f.write(f"accounting_group      = group_cms\n")
            f.write(f"should_transfer_files = YES\n")
            f.write(f"when_to_transfer_output = ON_EXIT\n")
            f.write(f"transfer_input_files  = {run_script_path}\n")
            # Arguments wrapped in quotes to handle complex paths/special characters
            f.write(f'Arguments             = "{main_path} {args.year} {args.study} {args.channel} {config_file} {branch_list} -1 {output_path} {sample}/{sample}_$(InputListName)"\n')
            # Queue command at the end with standard syntax
            f.write(f"Queue InputListName from (\n")
            for n in to_submit:
                f.write(f"{n}\n")
            f.write(f")\n")
        
        if args.check:
            print(f"[CHECK] Sample: {sample} | Jobs to submit: {len(to_submit)}")
        else:
            submit_condor_job(jdl_path)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Unified HTCondor Submission Script v1')
    parser.add_argument('--year', required=True, choices=['UL2016PreVFP', 'UL2016PostVFP', 'UL2017', 'UL2018'])
    parser.add_argument('--channel', required=True, choices=['MuMu', 'ElEl', 'MuEl'])
    parser.add_argument('--study', required=True, help='Study version (e.g., AN_v6p3)')
    parser.add_argument('--resubmit', action='store_true', help='Resubmit failed jobs only')
    parser.add_argument('--check', action='store_true', help='Check status without submitting')
    
    args = parser.parse_args()
    run_process(args)
