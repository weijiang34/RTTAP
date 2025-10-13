import os 
import shlex
import logging
logging.basicConfig(
    level=logging.DEBUG, 
    format="%(asctime)s [%(levelname)s]: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

import pandas as pd

import envs
from .system_operations import safe_run_cmd, safe_remove


def generate_overall_report(
    viruses_report_path, 
    bacteria_report_path, 
    fungi_report_path, 
    out_path,
    bac_rpm,
    vir_rpm,
    fun_rpm
):

    bt_viruses = pd.read_table(viruses_report_path, sep='\t', header=0)
    mpa_bacteria = pd.read_table(bacteria_report_path, sep='\t', header=0)
    bt_fungi = pd.read_table(fungi_report_path, sep='\t', header=0)

    # threshold filtering
    mpa_bacteria = mpa_bacteria[mpa_bacteria["RPM"]>=bac_rpm]
    bt_viruses = bt_viruses[bt_viruses["RPM"]>=vir_rpm]
    bt_fungi = bt_fungi[bt_fungi["RPM"]>=fun_rpm]

    ranks = ["SuperKingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]

    report = pd.DataFrame({})
    
    for idx, rank in enumerate(ranks):
        report = pd.concat([report, bt_viruses[bt_viruses["rank"]==rank]], axis=0)
        report = pd.concat([report, mpa_bacteria[mpa_bacteria["rank"]==rank]], axis=0)
        report = pd.concat([report, bt_fungi[bt_fungi["rank"]==rank]], axis=0)

    report.to_csv(out_path, sep='\t', index=None)

def run_rgi(
    input_fq, 
    out_dir, 
    fileHeader, 

    tool_path, 
    db_path, 
    conda_dir, 
    env_name,

    threads,
    cmd='',
):
    tmp_script_path = os.path.join(out_dir, f"ARG_indentification_tmp.sh")
    bash_commands = [
        "#!/bin/bash",
        f"exec 2>&1",
        f"source {conda_dir}/bin/activate {env_name}",
        f"cd {out_dir}",
        f"if [ -d ./localDB ];then",
        f"    rm -r ./localDB",
        f"fi",
        f"{tool_path} clean --local",
        f"{tool_path} load --card_json {db_path}/card.json --local",
        f"{tool_path} card_annotation -i {db_path}/card.json > card_annotation.log",
        f"{tool_path} load -i {db_path}/card.json --card_annotation card_database_v3.2.7.fasta --local",
        f"{tool_path} bwt -1 {input_fq} -n {threads} -o {out_dir}/{fileHeader}.rgi.out -a kma --local --clean {cmd} 2>&1",
        f"rm {os.path.join(out_dir, f'{fileHeader}.*.bam')}",
        f"rm {os.path.join(out_dir, f'{fileHeader}.*.bam.bai')}",
        f"{tool_path} clean --local",
    ]
    bash_commands = [line + '\n' for line in bash_commands]
    with open(tmp_script_path, 'w') as f:
        f.writelines(bash_commands)
    os.chmod(tmp_script_path, 0o755)
    safe_run_cmd(shlex.quote(tmp_script_path), log_file=os.path.join(out_dir, f"{fileHeader}.rgi.log"))
    safe_remove(tmp_script_path)

def run_virstrain(
    input_fq_gz, 
    out_dir, 
    fileHeader,

    tool_path, 
    db_path, 
    virus_report_path, 
    virstrain_list_path,
    conda_dir, 
    env_name,

    threads=1,
    cmd=''
):
    tmp_script_path = os.path.join(out_dir, f"virus_strain_profiling_tmp.sh")
    db_list = []

    report_virus = pd.read_table(virus_report_path, header=0, index_col=None, sep='\t')
    virstrain_list = pd.read_table(virstrain_list_path, sep='\t', header=0).set_index(["DB_name"])

    s_name_list = report_virus["taxa"].str.lower()
    # for viruses that only have name containing the db name
    for name_contains in virstrain_list["name_contains"].dropna().str.lower():
        if s_name_list[s_name_list.str.contains(name_contains)].shape[0] != 0:
            db_name = virstrain_list["name_contains"].loc[virstrain_list["name_contains"].str.lower()==name_contains].index.tolist()
            db_list += db_name
    # for viruses that have exact name match
    for species_name in virstrain_list["species_name"].dropna().str.lower():
        if species_name in s_name_list.tolist():
            db_name = virstrain_list["species_name"].loc[virstrain_list["species_name"].str.lower()==species_name].index.tolist()
            db_list += db_name
    db_list = pd.DataFrame({"DB_name":list(set(db_list))})
    db_list.to_csv(os.path.join(out_dir, f"{fileHeader}.virStrain_db_matched.txt"), sep='\t', index=False)

    # generate empty report file
    with open(os.path.join(out_dir, f"{fileHeader}.virstrain.report"), 'w') as f:
        f.writelines(["\n"])

    # run virstrain
    # if no strains in db, generate empty report and return
    if db_list.shape[0]==0:
        logging.warning("No strain found in VirStrain DB.")
        return 
    # bash header, activate env, decompress fq
    script_lines = [
        "#!/bin/bash",
        f"source {os.path.join(conda_dir, 'bin', 'activate')} {env_name}",
        f"cd {out_dir}",
        f"gzip -d -c {input_fq_gz} > {os.path.join(out_dir, f'{fileHeader}.viruses.fq')}"
    ]
    # call virstrain commands
    for db_name in db_list["DB_name"].tolist():
        script_lines.append(f"python {tool_path} -i {os.path.join(out_dir, f'{fileHeader}.viruses.fq')} -d {os.path.join(db_path, db_name)} -o {os.path.join(out_dir, db_name)} {cmd}")
        # process report if exists
        if os.path.exists(os.path.join(out_dir, fileHeader, db_name, "VirStrain_report.txt")):
            script_lines.append(f"grep '>>Most possible' -A1 {os.path.join(out_dir, fileHeader, db_name, 'VirStrain_report.txt')} | sed '1d' | sed 's/^[ \\t]*>//' | cut -f1,2,3,4,5,6,7,10 >> {os.path.join(out_dir, f'{fileHeader}.virstrain.report')}")
    script_lines.append(f"rm {os.path.join(out_dir, f'{fileHeader}.viruses.fq')}")
    script_lines = [line + '\n' for line in script_lines]
    
    # write and run the script
    with open(tmp_script_path, 'w') as f:
        f.writelines(script_lines)
    os.chmod(tmp_script_path, 0o755)
    safe_run_cmd(shlex.quote(tmp_script_path), log_file=os.path.join(out_dir, f"{fileHeader}.virstrain.log"))
    safe_remove(tmp_script_path)