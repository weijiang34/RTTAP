import os
import shlex

import logging
logging.basicConfig(
    level=logging.DEBUG, 
    format="%(asctime)s [%(levelname)s]: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

from .system_operations import safe_run_cmd, safe_remove

'''Functions for QC'''

def run_fastp(
    # input & output
    input_1, 
    input_2, 
    out_dir, 
    fileHeader, 
    # paths
    tool_path, 
    conda_dir,
    env_name,
    # parameter configs
    threads, 
    min_len=15, 
    cmd=''
):
    """Run fastp for quality control."""
    # if os.path.exists(os.path.join(out_dir, f"{fileHeader}.clean.R1.fq.gz")) or os.path.exists(os.path.join(out_dir, f"{fileHeader}.clean.fq.gz")):
    #     logging.warning(f"Output files already exist in {out_dir}. Skipping fastp step.")
    #     return
    script_lines = [
        "#!/bin/bash",
        f"source {os.path.join(conda_dir, 'bin', 'activate')} {env_name}",
    ]
    if input_2 is not None:
        script_lines.extend(
            [
                f"{tool_path} -i {shlex.quote(input_1)} -I {shlex.quote(input_2)} \\",
                f"    -o {shlex.quote(os.path.join(out_dir, f'{fileHeader}.clean.R1.fq.gz'))} -O {shlex.quote(os.path.join(out_dir, f'{fileHeader}.clean.R2.fq.gz'))} \\",
                f"    --json {shlex.quote(os.path.join(out_dir, f'{fileHeader}.json'))} --html {shlex.quote(os.path.join(out_dir, f'{fileHeader}.html'))} \\",
                f"    --thread {threads} --length_required {min_len} -D {cmd}",
            ]
        )
    else:
        script_lines.extend(
            [
                f"{tool_path} -i {shlex.quote(input_1)} \\",
                f"    -o {shlex.quote(os.path.join(out_dir, f'{fileHeader}.clean.fq.gz'))} \\",
                f"    --json {shlex.quote(os.path.join(out_dir, f'{fileHeader}.json'))} --html {shlex.quote(os.path.join(out_dir, f'{fileHeader}.html'))} \\",
                f"    --thread {threads} --length_required {min_len} -D {cmd}",
            ]
        )
    script_lines = [line + '\n' for line in script_lines]
    tmp_script_path = os.path.join(out_dir, f"run_fastp_tmp.sh")
    with open(tmp_script_path, 'w') as f:
        f.writelines(script_lines)
    os.chmod(tmp_script_path, 0o755)
    safe_run_cmd(shlex.quote(tmp_script_path), log_file=os.path.join(out_dir, f"{fileHeader}.fastp.log"))
    safe_remove(tmp_script_path)
