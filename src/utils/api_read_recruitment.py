import os
import shlex

import envs
from .system_operations import safe_run_cmd, safe_remove

def run_kraken2(
    input, out_dir, tool_path, db_path, conda_dir, env_name, fileHeader, threads, cmd=''
):
    """Run Kraken2 for taxonomic classification."""
    kraken_out = os.path.join(out_dir, f"{fileHeader}.kraken2.out")
    kraken_report = os.path.join(out_dir, f"{fileHeader}.kraken2.report")
    tmp_script_path = os.path.join(out_dir, f"run_kraken2_tmp.sh")

    ### Consider move elsewhere
    # if os.path.exists(kraken_out) or os.path.exists(kraken_report):
    #     logging.warning(f"Output files already exist in {out_dir}. Skipping Kraken2 step.")
    #     return

    script_lines = [
        "#!/bin/bash",
        f"source {os.path.join(conda_dir, 'bin', 'activate')} {env_name}",
        f"{tool_path} --db {shlex.quote(db_path)} --report {shlex.quote(kraken_report)} --output {shlex.quote(kraken_out)} \\",
        f"    --use-names --threads {threads} --gzip-compressed {shlex.quote(input)} {cmd}",
    ]
    script_lines = [line + '\n' for line in script_lines]
    # Write, run, and timing the script
    with open(tmp_script_path, 'w') as f:
        f.writelines(script_lines)
    os.chmod(tmp_script_path, 0o755)
    safe_run_cmd(shlex.quote(tmp_script_path), log_file=os.path.join(out_dir, f"{fileHeader}.kraken2.log"))
    safe_remove(tmp_script_path)