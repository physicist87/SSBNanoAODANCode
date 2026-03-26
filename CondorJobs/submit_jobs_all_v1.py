import os
import sys
import subprocess
import re
import time
import argparse
from collections import defaultdict

# --- Common Helper Functions ---

def submit_condor_job(submit_file_path):
    """Submit a job to HTCondor with appropriate settings"""
    user_home = os.path.expanduser("~")
    is_kisti = user_home.startswith("/cms/ldap_home")

    if is_kisti:
        uid = os.getuid()
        proxy_path = f"/tmp/x509up_u{uid}"
        if not os.path.isfile(proxy_path):
            print(f"[ERROR] Proxy file not found: {proxy_path}")
            sys.exit(1)

    cmd = ["condor_submit"]
    if is_kisti:
        cmd += ["-append", "accounting_group = group_cms"]
    cmd.append(submit_file_path)

    print(f"[INFO] Submitting job: {' '.join(cmd)}")
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] condor_submit failed: {e}")
        sys.exit(1)

def get_job_config(year, channel, sample_name):
    """
    Unified logic to determine config file and branch list based on year and sample
    """
    base_config_map = {"MuMu": "dimuon", "ElEl": "dielec", "MuEl": "muelec"}
    base_conf = base_config_map[channel]
    
    # Default values
    config_file = f"{base_conf}.config"
    branch_list = f"{year}/branch_list.txt"

    # Year-specific overrides
    if year == "UL2017" and sample_name.startswith("Data_"):
        if "Run2017B" in sample_name:
            config_file = f"{base_conf}_Data_RunB.config"
            branch_list = f"{year}/branch_list_Run2017B.txt"
        else:
            config_file = f"{base_conf}_Data_RunCtoF.config"
            branch_list = f"{year}/branch_list_Run2017CtoF.txt"
            
    elif year == "UL2016PostVFP" and sample_name.startswith("Data_"):
        if "Run2016H" in sample_name:
            config_file = f"{base_conf}_Data_RunH.config"
            branch_list = f"{year}/branch_list_RunH.txt"
            
    # UL2018 and UL2016PreVFP use default logic
    return config_file, branch_list

# --- Status & Check Functions ---

def check_job_status(logDir, studyName, runPeriod, channel, sampleDir, listFileName):
    """Check the status of a specific job"""
    jobLogDir = os.path.join(logDir, studyName, runPeriod, channel, sampleDir)
    log_prefixes = ["resubmit_", ""]
    
    for prefix in log_prefixes:
        err_file = os.path.join(jobLogDir, f"{prefix}{sampleDir}_{listFileName}.err")
        out_file = os.path.join(jobLogDir, f"{prefix}{sampleDir}_{listFileName}.out")
        log_file = os.path.join(jobLogDir, f"{prefix}{sampleDir}_{listFileName}.log")
        
        if not (os.path.exists(err_file) or os.path.exists(out_file) or os.path.exists(log_file)):
            continue
        
        if os.path.exists(log_file) and not (os.path.exists(err_file) and os.path.exists(out_file)):
            return 'not_started'
            
        if os.path.exists(err_file) and os.path.getsize(err_file) > 0:
            return 'failed'
            
        if os.path.exists(out_file):
            try:
                with open(out_file, 'r') as f:
                    content = f.read()
                if "Total number of events after merging root files: 0" in content:
                    return 'zero_events'
                
                required_endings = ["fout successfully deleted.", "Analysis destructor completed."]
                if all(ending in content for ending in required_endings):
                    return 'completed'
            except:
                return 'failed'
    return 'not_submitted'

# --- Main Submission Logic ---

def run_process(args):
    """Main execution function handling submit, resubmit, or check"""
    this_script_dir = os.path.dirname(os.path.abspath(__file__))
    MainPath = os.path.abspath(os.path.join(this_script_dir, ".."))
    
    inputListPath = os.path.join(MainPath, f"input/{args.year}")
    logDir = os.path.join(os.getcwd(), "condorLog")
    submitDir = os.path.join(os.getcwd(), "condorSubmit", args.study, args.year, args.channel)
    runScriptPath = os.path.join(MainPath, "run_cmd_v7.sh")
    outputPath = f"/pnfs/knu.ac.kr/data/cms/store/user/sha/CPV_Run2/ULSummer20"
    
    os.makedirs(submitDir, exist_ok=True)

    # Simplified sample selection: list directories in input path
    all_samples = [d for d in os.listdir(inputListPath) if os.path.isdir(os.path.join(inputListPath, d))]
    
    for sample in all_samples:
        # Determine job config dynamically
        config_file, branch_list = get_job_config(args.year, args.channel, sample)
        
        sampleInputPath = os.path.join(inputListPath, sample)
        jobLogDir = os.path.join(logDir, args.study, args.year, args.channel, sample)
        os.makedirs(jobLogDir, exist_ok=True)

        # Collect list files
        list_nums = []
        for f in os.listdir(sampleInputPath):
            if f.endswith(".list"):
                match = re.search(r'_(\d+)\.list$', f)
                if match: list_nums.append(match.group(1))

        if args.check:
            # Just print status and move to next sample
            stats = defaultdict(int)
            for n in list_nums:
                status = check_job_status(logDir, args.study, args.year, args.channel, sample, n)
                stats[status] += 1
            print(f"Sample: {sample} -> {dict(stats)}")
            continue

        # Submission logic (Normal or Resubmit)
        to_submit = []
        if args.resubmit:
            to_submit = [n for n in list_nums if check_job_status(logDir, args.study, args.year, args.channel, sample, n) != 'completed']
        else:
            to_submit = list_nums

        if not to_submit: continue

        # Build Condor JDL
        prefix = "resubmit_" if args.resubmit else ""
        jdl_path = os.path.join(submitDir, f"{prefix}{sample}.sub")
        
        with open(jdl_path, "w") as f:
            f.write(f"Universe = vanilla\nExecutable = {runScriptPath}\n")
            f.write(f"Log = {jobLogDir}/{prefix}{sample}_$(In).log\n")
            f.write(f"Output = {jobLogDir}/{prefix}{sample}_$(In).out\n")
            f.write(f"Error = {jobLogDir}/{prefix}{sample}_$(In).err\n")
            f.write(f"Arguments = {MainPath} {args.year} {args.study} {args.channel} {config_file} {branch_list} -1 {outputPath} {sample}/{sample}_$(In)\n")
            f.write(f"Queue In from (\n" + "\n".join(to_submit) + "\n)\n")
        
        submit_condor_job(jdl_path)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--year', required=True, choices=['UL2016PreVFP', 'UL2016PostVFP', 'UL2017', 'UL2018'])
    parser.add_argument('--channel', required=True, choices=['MuMu', 'ElEl', 'MuEl'])
    parser.add_argument('--study', required=True, help="e.g., AN_v6p2-13")
    parser.add_argument('--resubmit', action='store_true', help="Resubmit failed jobs only")
    parser.add_argument('--check', action='store_true', help="Check job status without submitting")
    
    args = parser.parse_args()
    run_process(args)
