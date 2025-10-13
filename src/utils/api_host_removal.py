import os
import shlex

from .system_operations import safe_run_cmd, safe_remove

def remove_host(
    input_1, 
    out_dir, 
    fileHeader, 
    tool_path, 
    conda_dir, 
    env_name, 
    threads, 
    db_dir, 
    input_2=None, 
    cmd=''
):
        tmp_script_path = os.path.join(out_dir, f"remove_host_tmp.sh")
        script_lines = [
            "#!/bin/bash",
            f"source {os.path.join(conda_dir, 'bin', 'activate')} {env_name}",
        ]
        if input_2 is not None:
            un_conc_path = os.path.join(out_dir, f"{fileHeader}.nohost.fq")
            bt2_out_path = os.path.join(out_dir, f"{fileHeader}.nohost.bowtie2.out")
            script_lines.extend([
                f"{tool_path} -x {shlex.quote(db_dir)} -1 {shlex.quote(input_1)} -2 {shlex.quote(input_2)} \\",
                f"    --un-conc {shlex.quote(un_conc_path)} -S {shlex.quote(bt2_out_path)} \\", 
                f"    --very-sensitive-local --no-unal -I 1 \\", 
                f"    -p {threads} {cmd}",
                f"mv {os.path.join(out_dir, f'{fileHeader}.nohost.1.fq')} {os.path.join(out_dir, f'{fileHeader}.nohost.R1.fq')}",
                f"mv {os.path.join(out_dir, f'{fileHeader}.nohost.2.fq')} {os.path.join(out_dir, f'{fileHeader}.nohost.R2.fq')}",
                f"cat {os.path.join(out_dir, f'{fileHeader}.nohost.R1.fq')} {os.path.join(out_dir, f'{fileHeader}.nohost.R2.fq')} | sed 's/ /_/g' | gzip -f -c -1 > {os.path.join(out_dir, f'{fileHeader}.nohost.fq.gz')}",
                f"rm {os.path.join(out_dir, f'{fileHeader}.nohost.R1.fq')} {os.path.join(out_dir, f'{fileHeader}.nohost.R2.fq')}",
                f"rm {bt2_out_path}",
            ])
        else:
            nohost_path = os.path.join(out_dir, f"{fileHeader}.nohost.fq")
            nohost_gz_path = os.path.join(out_dir, f"{fileHeader}.nohost.fq.gz")
            bt2_out_path = os.path.join(out_dir, f"{fileHeader}.nohost.bowtie2.out")
            script_lines.extend([
                f"{tool_path} -x {shlex.quote(db_dir)} -U {shlex.quote(input_1)} \\",
                f"    --un {shlex.quote(nohost_path)} -S {shlex.quote(bt2_out_path)} \\", 
                f"    --very-sensitive-local --no-unal -I 1 \\", 
                f"    -p {threads} {cmd}",
                f"cat {nohost_path} | sed 's/ /_/g' | gzip -f -c -1 > {nohost_gz_path}",
                f"rm {nohost_path}",
                f"rm {bt2_out_path}",
            ])
        script_lines = [line + '\n' for line in script_lines]
        with open(tmp_script_path, 'w') as f:
            f.writelines(script_lines)
        os.chmod(tmp_script_path, 0o755)
        safe_run_cmd(shlex.quote(tmp_script_path), log_file=os.path.join(out_dir, f"{fileHeader}.remove_host.log"))
        safe_remove(tmp_script_path)