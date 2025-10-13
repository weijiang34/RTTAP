import os 
import shlex
import pandas as pd
import logging
logging.basicConfig(
    level=logging.DEBUG, 
    format="%(asctime)s [%(levelname)s]: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

import envs
from .system_operations import safe_run_cmd, safe_remove


def run_metaphlan4(
    input, 
    out_dir, 
    fileHeader, 

    tool_path, 
    metaphlan4_db_path, 
    
    threads=1, 
    cmd=''
):
    """Run MetaPhlAn4 for taxonomic classification."""
    mpa_bt_out = os.path.join(out_dir, f"{fileHeader}.metaphlan4.bowtie2.out")
    mpa_out = os.path.join(out_dir, f"{fileHeader}.metaphlan4.report")
    script_path = os.path.join(out_dir, f"run_metaphlan4_tmp.sh")

    if os.path.exists(mpa_bt_out):
        safe_remove(mpa_bt_out)

    script_lines = [
        "#!/bin/bash",
        f"source {os.path.join(envs.CONDA_PATH, 'bin', 'activate')} {envs.RTTAP_ENV_NAME}",
        f"{tool_path} {shlex.quote(input)} --bowtie2db {shlex.quote(os.path.dirname(metaphlan4_db_path))} -x {shlex.quote(os.path.basename(metaphlan4_db_path))} -t rel_ab_w_read_stats -o {shlex.quote(mpa_out)} --nproc {threads} --input_type fastq --CAMI_format_output --bowtie2out {shlex.quote(mpa_bt_out)} --unclassified_estimation {cmd}"
    ]
    script_lines = [line + '\n' for line in script_lines]
    with open(script_path, 'w') as f:
        f.writelines(script_lines)
    os.chmod(script_path, 0o755)
    safe_run_cmd(shlex.quote(script_path), log_file=os.path.join(out_dir, f"{fileHeader}.metaphlan4.log"))
    safe_remove(script_path)

def generate_mpa_report(
    k2_report_path,
    bac_report_path,
    mpa_report_path, 
    RPM_threshold=1
):
    """Generate final taxonomic report."""
    # detect header line
    def find_metaphlan_header_line(filepath):
        with open(filepath, 'r') as f:
            for i, line in enumerate(f):
                if "#clade_name" in line:
                    return i
        raise ValueError("No header line with '#clade_name' found.")

    header_line = find_metaphlan_header_line(mpa_report_path)
    mpa_report = pd.read_table(mpa_report_path, header=header_line, sep='\t').rename({"#clade_name":"taxaPath", "estimated_number_of_reads_from_the_clade":"reads", "clade_taxid":"taxaIDPath"}, axis=1)
    k2_report = pd.read_table(k2_report_path, sep='\t',header=None,names=["percent","reads","reads_direct","rank","taxid","taxa"])
    total_reads_count = k2_report[k2_report["taxa"].isin(["root","unclassified"])]["reads"].sum(axis=0)
    
    # remove unclassified
    if mpa_report.shape[0]==0 or mpa_report[mpa_report['taxaPath']!='unclassified'].shape[0]==1:
        report = pd.DataFrame({
            "rank":[],
            "taxa":[],
            "taxid":[],
            "reads":[],
            "RPM":[]
        })
        logging.warning("No bacteria found in MetaPhlAn4 results.")
        report.to_csv(bac_report_path, sep='\t', index=None)
        return
    else:
        mpa_report = mpa_report.iloc[1:,:].reset_index(drop=True)
        # pre-process
        mpa_report["taxaPath"] = mpa_report["taxaPath"].str.replace("|", ";", case=True, regex=False)
        mpa_report["taxaIDPath"] = mpa_report["taxaIDPath"].str.replace("|", ";", case=True, regex=False)

        report = pd.DataFrame({})
        report_tmp = pd.DataFrame({})

        taxaPath_tmp = mpa_report["taxaPath"].str.split(";", expand=True)
        taxaIDPath_tmp = mpa_report["taxaIDPath"].str.split(";", expand=True)
        taxa = []
        taxaID = []
        for idx in range(mpa_report.shape[0]):
            taxa.append(taxaPath_tmp.iloc[idx,:][taxaPath_tmp.iloc[idx,:].T.last_valid_index()])
            taxaID.append(taxaIDPath_tmp.iloc[idx,:][taxaIDPath_tmp.iloc[idx,:].T.last_valid_index()])
        report_tmp["taxa"] = pd.Series(taxa).str.split("__", expand=True)[1]
        report_tmp["taxid"] = taxaID
        report_tmp["rank"] = pd.Series(taxa).str.split("__", expand=True)[0]
        report_tmp["reads"] = mpa_report["reads"]

        ranks = ["SuperKingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain"]
        for idx, item in enumerate(["k","p","c","o","f","g","s","t"]):
            report_tmp["rank"] = report_tmp["rank"].replace(item, ranks[idx], regex=False)
        report_tmp = report_tmp.groupby(["rank", "taxa", "taxid"]).agg({"reads":"sum"})

        for idx, rank in enumerate(ranks[:7]):
            report = pd.concat([report, report_tmp[report_tmp.index.get_level_values("rank")==rank].sort_values("reads", ascending=False)], axis=0)
        
        report = report.reset_index()
        report["taxa"] = report["taxa"].str.replace('_', ' ')
        report["RPM"] = report["reads"]*1000000/total_reads_count
        report = report[report["RPM"]>=RPM_threshold]
        report.to_csv(bac_report_path, sep='\t', index=None)