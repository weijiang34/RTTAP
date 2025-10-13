import os 
import shlex
import logging
logging.basicConfig(
    level=logging.DEBUG, 
    format="%(asctime)s [%(levelname)s]: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

import envs
from .system_operations import safe_run_cmd, safe_remove

def split_reads(
    out_dir, 
    fileHeader, 
    kraken_out, 
    kraken_report, 
    nohost_fq_gz, 
    reads_type, 
    conda_dir, 
    env_name,
    split_script_path
):  
    taxid_dict = {
        "bacteria": 2,
        "viruses": 10239,
        "fungi": 4751
    }
    # logging.info(f"Splitting {reads_type} reads")
    taxid = taxid_dict.get(reads_type.lower())
    fq = os.path.join(out_dir, f"{fileHeader}.{reads_type}.fq")
    tmp_script_path = os.path.join(out_dir, f"split_{reads_type}_tmp.sh")

    script_lines = [
        "#!/bin/bash",
        f"source {os.path.join(conda_dir, 'bin', 'activate')} {env_name}",
        f"python {split_script_path} -k {kraken_out} -s {nohost_fq_gz} -r {kraken_report} -t {taxid} -o {fq} --include-children --fastq-output",
        f"gzip --fast {fq} -f"
    ]
    script_lines = [line + '\n' for line in script_lines]

    # write, run, and timing the script
    with open(tmp_script_path, 'w') as f:
        f.writelines(script_lines)
    os.chmod(tmp_script_path, 0o755)
    safe_run_cmd(shlex.quote(tmp_script_path), log_file=os.path.join(out_dir, f"{fileHeader}.{reads_type}.split.log"))
    safe_remove(tmp_script_path)