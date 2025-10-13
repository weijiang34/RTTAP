# All Chinese characters in the following code need to be translated into English
import os
import sys
import shlex
import subprocess as sp
import multiprocessing

import logging
logging.basicConfig(
    level=logging.DEBUG, 
    format="%(asctime)s [%(levelname)s]: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
import time
import glob

import pandas as pd
import numpy as np
import dask.dataframe as dd

from .. import envs


def remove_host(input_1, out_dir, tool_path, db_path, fileHeader, threads=1, input_2=None, cmd=None):
    if cmd is None:
        cmd = ''
    step_name = 'Remove host reads'
    os.makedirs(out_dir, exist_ok=True)
    tmp_script_path = os.path.join(out_dir, f"{'_'.join(step_name.lower().split(' ')+['_tmp'])}.sh")
    logging.info(f"Starting {step_name.lower()} ...")

    scripts = []
    # bash header
    scripts.append("#!/bin/bash\n")
    # activate conda environment
    scripts.append(f"source {os.path.join(envs.CONDA_PATH, 'bin', 'activate')} {envs.RTTAP_ENV_NAME}\n")
    # execute bowtie2
    if input_2 is not None:
        un_conc_path = os.path.join(out_dir, f"{fileHeader}.nohost.fq")
        bt2_out_path = os.path.join(out_dir, f"{fileHeader}.nohost.bowtie2.out")
        scripts.append(f"{tool_path} -x {db_path} -1 {input_1} -2 {input_2} --un-conc {un_conc_path} -S {bt2_out_path} --very-sensitive-local --no-unal -I 1 -p {threads} {cmd}\n")
        scripts.extend([
            f"mv {os.path.join(out_dir, f'{fileHeader}.nohost.1.fq')} {os.path.join(out_dir, f'{fileHeader}.nohost.R1.fq')}\n",
            f"mv {os.path.join(out_dir, f'{fileHeader}.nohost.2.fq')} {os.path.join(out_dir, f'{fileHeader}.nohost.R2.fq')}\n",
            f"cat {os.path.join(out_dir, f'{fileHeader}.nohost.R1.fq')} {os.path.join(out_dir, f'{fileHeader}.nohost.R2.fq')} | sed 's/ /_/g' | gzip -f -c -1 > {os.path.join(out_dir, f'{fileHeader}.nohost.fq.gz')}\n",
            f"rm {os.path.join(out_dir, f'{fileHeader}.nohost.R1.fq')} {os.path.join(out_dir, f'{fileHeader}.nohost.R2.fq')}\n",
            f"rm {bt2_out_path}\n",
        ])
    else:
        nohost_path = os.path.join(out_dir, f"{fileHeader}.nohost.fq")
        nohost_gz_path = os.path.join(out_dir, f"{fileHeader}.nohost.fq.gz")
        bt2_out_path = os.path.join(out_dir, f"{fileHeader}.nohost.bowtie2.out")
        scripts.append(f"{tool_path} -x {db_path} -U {input_1} --un {nohost_path} -S {bt2_out_path} --very-sensitive-local --no-unal -I 1 -p {threads} {cmd}\n")
        scripts.extend([
            f"cat {nohost_path} | sed 's/ /_/g' | gzip -f -c -1 > {nohost_gz_path}\n",
            f"rm {nohost_path}\n",
            f"rm {bt2_out_path}\n",
        ])

    # write, run, and timing the script
    with open(tmp_script_path, 'w') as f:
        f.writelines(scripts)
    os.system(f"chmod +x {tmp_script_path}")
    start_time = time.time()
    os.system(tmp_script_path)
    elapsed_time = time.time() - start_time
    logging.info(f"{step_name.capitalize()} finished in {elapsed_time:.2f} seconds. Output to {out_dir}/")
    os.remove(tmp_script_path)

def run_fastp(input_1, out_dir, tool_path, fileHeader, threads=1, input_2=None, min_len=15, cmd=None):
    if cmd is None:
        cmd = ''
    step_name = 'Quality control'
    os.makedirs(out_dir, exist_ok=True)
    tmp_script_path = os.path.join(out_dir, f"{'_'.join(step_name.lower().split(' ')+['_tmp'])}.sh")
    logging.info(f"Starting {step_name.lower()} ...")

    scripts = []
    # bash header 
    scripts.append("#!/bin/bash\n")
    # activate conda environment
    scripts.append(f"source {os.path.join(envs.CONDA_PATH, 'bin', 'activate')} {envs.RTTAP_ENV_NAME}\n")
    # execute fastp
    if input_2 is not None:
        scripts.append(f"{tool_path} -i {input_1} -I {input_2} -o {out_dir}/{fileHeader}.clean.R1.fq.gz -O {out_dir}/{fileHeader}.clean.R2.fq.gz --json {out_dir}/{fileHeader}.json --html {out_dir}/{fileHeader}.html --thread {threads} --length_required {min_len} -D {cmd}\n")
    else:
        scripts.append(f"{tool_path} -i {input_1} -o {out_dir}/{fileHeader}.clean.fq.gz --json {out_dir}/{fileHeader}.json --html {out_dir}/{fileHeader}.html --thread {threads} --length_required {min_len} -D\n")

    # write, run, and timing the script
    with open(tmp_script_path, 'w') as f:
        f.writelines(scripts)
    os.system(f"chmod +x {tmp_script_path}")
    start_time = time.time()
    os.system(tmp_script_path)
    elapsed_time = time.time() - start_time
    logging.info(f"{step_name.capitalize()} finished in {elapsed_time:.2f} seconds. Output to {out_dir}/")
    os.remove(tmp_script_path)

def split_reads(kraken_out, kraken_report, fq_gz, out_dir, fileHeader, reads_type, script_path="./extract_kraken_reads_nostdout.py"):
    taxid_dict = {
        "bacteria": 2,
        "viruses": 10239,
        "fungi": 4751
    }
    if reads_type.lower() not in taxid_dict.keys():
        raise ValueError(f"Invalid reads_type: {reads_type}. Must be one of {list(taxid_dict.keys())}.")
    taxid = taxid_dict.get(reads_type.lower())

    step_name = f"Split {reads_type.lower()} reads"
    os.makedirs(out_dir, exist_ok=True)
    fq = os.path.join(out_dir, f"{fileHeader}.{reads_type}.fq")
    tmp_script_path = os.path.join(out_dir, f"{'_'.join(step_name.lower().split(' ')+['_tmp'])}.sh")

    logging.info(f"Starting {step_name} ...")
    scripts = []
    # bash header
    scripts.append("#!/bin/bash\n")
    # activate conda environment
    scripts.append(f"source {os.path.join(envs.CONDA_PATH, 'bin', 'activate')} {envs.RTTAP_ENV_NAME}\n")
    # execute python script
    scripts.append(f"python {script_path} -k {kraken_out} -s {fq_gz} -r {kraken_report} -t {taxid} -o {fq} --include-children --fastq-output\n")
    scripts.append(f"gzip --fast {fq} -f\n")

    # write, run, and timing the script
    with open(tmp_script_path, 'w') as f:
        f.writelines(scripts)
    os.system(f"chmod +x {tmp_script_path}")
    start_time = time.time()
    os.system(tmp_script_path)
    elapsed_time = time.time() - start_time
    logging.info(f"{step_name.capitalize()} finished in {elapsed_time:.2f} seconds. Output to {out_dir}/")
    os.remove(tmp_script_path)
    return


def run_bowtie2(input, out_dir, bt_path, taxonkit_path, fileHeader, db_type, acc2taxid_dir, threads=1):
    db_path_dict = {
        "viruses": envs.BT_VIRAL_DB_PATH,
        "fungi": envs.BT_FUNGAL_DB_PATH
    }
    acc2taxid_dict = {
        "viruses": f"{acc2taxid_dir}/viruses.acc2taxid.txt",
        "fungi": f"{acc2taxid_dir}/EuPathDB46.seq2taxid.txt",
    }
    if db_type not in db_path_dict.keys():
        raise ValueError(f"Invalid db_type: {db_type}. Must be one of {list(db_path_dict.keys())}.")
    
    step_name = f"Bowtie2 mapping for {db_type}"
    os.makedirs(out_dir, exist_ok=True)
    db_path = db_path_dict.get(db_type)
    acc2taxid_path = acc2taxid_dict.get(db_type)

    logging.info(f"Starting {step_name.lower()} ...")
    # Step 1-2: run bowtie2 and extract uniq and multi mapping readID-accID records using awk
    scripts = []
    # bash header
    scripts.append("#!/bin/bash\n")
    # activate conda environment
    scripts.append(f"source {os.path.join(envs.CONDA_PATH, 'bin', 'activate')} {envs.RTTAP_ENV_NAME}\n")
    # execute bowtie2
    scripts.append(f"{bt_path} -x {db_path} -U {input}  --xeq -S {out_dir}/{fileHeader}.{db_type}.bowtie2.all.out -a -p {threads} --no-unal --no-head\n")
    # extract uniq and multi mapping reads
    scripts.append(
        f"awk '{{if ($0 ~ /AS:/ && $0 !~ /XS:/) print $1 \"\\t\" $3 > \"{out_dir}/{fileHeader}.{db_type}.uniq.tmp\"; "
        f"else if ($0 ~ /AS:/ && $0 ~ /XS:/) print $1 \"\\t\" $3 > \"{out_dir}/{fileHeader}.{db_type}.multi.tmp\"}}' "
        f"{out_dir}/{fileHeader}.{db_type}.bowtie2.all.out\n"
    )
    # write, run, and timing the script
    tmp_script_path = os.path.join(out_dir, f'run_Bowtie2_{db_type}_tmp.sh')
    with open(tmp_script_path, 'w') as f:
        f.writelines(scripts)
    os.system(f"chmod +x {tmp_script_path}")
    start_time = time.time()
    os.system(tmp_script_path)
    elapsed_time = time.time() - start_time
    logging.info(f"{step_name.capitalize()} finished in {elapsed_time:.2f} seconds. Output to {out_dir}/")
    os.remove(tmp_script_path)

    # Step 3: append taxid to uniq and multi mapping results
    acc2taxid_dict = {
        "viruses": f"{acc2taxid_dir}/viruses.acc2taxid.txt",
        "fungi": f"{acc2taxid_dir}/EuPathDB46.seq2taxid.txt",
    }
    acc2taxid_path = acc2taxid_dict.get(db_type)
    uniq_read2acc_path = f"{out_dir}/{fileHeader}.{db_type}.uniq.tmp"
    multi_read2acc_path = f"{out_dir}/{fileHeader}.{db_type}.multi.tmp"
    step_name = f"Appending taxid to bowtie2 {db_type} results"
    # read files
    logging.info(f"Start {step_name.lower()} ...")
    start_time = time.time()
    acc2taxid = pd.read_table(acc2taxid_path, header=None, names=["acc.ver", "taxid"]).astype({"taxid":"str"})
    uniqread2acc = pd.read_table(f"{out_dir}/{fileHeader}.{db_type}.uniq.tmp", header=None, names=["readId", "acc.ver"])
    multiread2acc = pd.read_table(f"{out_dir}/{fileHeader}.{db_type}.multi.tmp", header=None, names=["readId", "acc.ver"])
    # add taxid after acc.ver column
    uniqacc2taxid = pd.merge(uniqread2acc, acc2taxid, how="left", on="acc.ver")
    multiacc2taxid = pd.merge(multiread2acc, acc2taxid, how="left", on="acc.ver")
    # save temporary files
    uniqacc2taxid.to_csv(f"{out_dir}/{fileHeader}.{db_type}.uniq.reads2acc2taxid.tmp", header=None, index=False, sep='\t')
    multiacc2taxid.to_csv(f"{out_dir}/{fileHeader}.{db_type}.multi.reads2acc2taxid.tmp", header=None, index=False, sep='\t')
    elapsed_time = time.time() - start_time
    logging.info(f"{step_name.capitalize()} finished in {elapsed_time:.2f} seconds.")
    
    # Step 4: convert taxid to taxonomic path using taxonkit
    taxonkit_cmds = [
        f"cat {out_dir}/{fileHeader}.{db_type}.uniq.reads2acc2taxid.tmp | {taxonkit_path} reformat -I3 -t > {out_dir}/{fileHeader}.{db_type}.uniq.reads2taxid_path.tmp\n",
        f"cat {out_dir}/{fileHeader}.{db_type}.multi.reads2acc2taxid.tmp | {taxonkit_path} reformat -I3 -t > {out_dir}/{fileHeader}.{db_type}.multi.reads2taxid_path.tmp\n",
    ]
    step_name = f"Converting taxid to taxonomic path using taxonkit for {db_type}"
    logging.info(f"Starting {step_name.lower()} ...")
    start_time = time.time()
    for cmd in taxonkit_cmds:
        sp.run(cmd, shell=True, check=True)
    elapsed_time = time.time() - start_time
    logging.info(f"{step_name.capitalize()} finished in {elapsed_time:.2f} seconds.")

    # Step 5: run LCA on bowtie2 results
    def lineage_lca(lineage):
        splitted_lineage = lineage.str.split(";")
        splitted_lineage = splitted_lineage.apply(lambda x: [i for i in x if pd.notnull(i)])
        arrs = np.array(splitted_lineage.tolist())
        # 找每列都一致的部分
        mask = (arrs == arrs[0]).all(axis=0)
        lca = ";".join(arrs[0][mask])
        return lca
    step_name = f"Running LCA on bowtie2 {db_type} results"
    logging.info(f"Starting {step_name.lower()} ...")
    start_time = time.time()
    ddf = dd.read_csv(f"{out_dir}/{fileHeader}.{db_type}.multi.reads2taxid_path.tmp", sep='\t', names=["readId", "refId", "taxid", "namePath", "taxaPath"])
    result = ddf.groupby("readId").agg({"taxaPath": lambda x: lineage_lca(x)}).compute()
    result.to_csv(f"{out_dir}/{fileHeader}.{db_type}.multi_reads_path.tmp", header=None, index=None, sep='\t')
    elapsed_time = time.time() - start_time
    logging.info(f"{step_name.capitalize()} finished in {elapsed_time:.2f} seconds.")
    
    # Step 6: merge uniq and multi mapping results, clean up temporary files
    step_name = f"Merging results"
    logging.info(f"Starting {step_name.lower()} ...")
    start_time = time.time()
    uniq_results = pd.read_csv(f"{out_dir}/{fileHeader}.{db_type}.uniq.reads2taxid_path.tmp", sep='\t', header=None)
    multi_results = pd.read_csv(f"{out_dir}/{fileHeader}.{db_type}.multi.reads2taxid_path.tmp", sep='\t', header=None)
    merged_results = pd.concat([uniq_results, multi_results]).drop_duplicates()
    merged_results.to_csv(f"{out_dir}/{fileHeader}.{db_type}.LCA.out", header=None, index=None, sep='\t')
    elapsed_time = time.time() - start_time
    logging.info(f"{step_name.capitalize()} finished in {elapsed_time:.2f} seconds.")
    # clean up temporary files
    temp_files = [
        f"{out_dir}/{fileHeader}.{db_type}.uniq.tmp",
        f"{out_dir}/{fileHeader}.{db_type}.multi.tmp",
        f"{out_dir}/{fileHeader}.{db_type}.uniq.reads2acc2taxid.tmp",
        f"{out_dir}/{fileHeader}.{db_type}.multi.reads2acc2taxid.tmp",
        f"{out_dir}/{fileHeader}.{db_type}.uniq.reads2taxid_path.tmp",
        f"{out_dir}/{fileHeader}.{db_type}.multi.reads2taxid_path.tmp",
    ]
    for temp_file in temp_files:
        if os.path.exists(temp_file):
            os.remove(temp_file)

    return

    if db_type=="viruses":
        logging.info("Running Bowtie2 on viral reads ...")
        # print("INFO: Run Bowtie2 on viral reads ...")
        bash_commands_1 = [
            f"{bt_path} -x {db_path} -U {input}  --xeq -S {out_dir}/{fileHeader}.viruses.bowtie2.all.out -a -p {threads} --no-unal --no-head\n",
            f"grep \"AS:\" {out_dir}/{fileHeader}.viruses.bowtie2.all.out | grep \"XS:\" -v | cut -f1,3 > {out_dir}/{fileHeader}.viruses.uniq.tmp\n",
            f"grep \"AS:\" {out_dir}/{fileHeader}.viruses.bowtie2.all.out | grep \"XS:\" | cut -f1,3 > {out_dir}/{fileHeader}.viruses.multi.tmp\n"
        ]
        with open(f"{out_dir}/bt_1v_tmp.sh", 'w') as f:
            f.writelines("#!/bin/bash\n")
            f.writelines(bash_commands_1)
        os.system(f"chmod +x {out_dir}/bt_1v_tmp.sh")
        start_time = time.time()
        os.system(f"{out_dir}/bt_1v_tmp.sh")
        elapsed_time = time.time() - start_time
        logging.info(f"Bowtie2 on viral reads finished in {elapsed_time:.2f} seconds.")
        os.remove(f"{out_dir}/bt_1v_tmp.sh")
        
        append_taxid(bt_results_path=out_dir, acc2taxid_path=acc2taxid_path, reads_type=db_type, fileHeader=fileHeader)
        
        bash_commands_2 = [
            f"cat {out_dir}/{fileHeader}.viruses.uniq.reads2acc2taxid.tmp | {taxonkit_path} reformat -I3 -t > {out_dir}/{fileHeader}.viruses.uniq.reads2taxid_path.tmp\n",
            f"cat {out_dir}/{fileHeader}.viruses.multi.reads2acc2taxid.tmp | {taxonkit_path} reformat -I3 -t > {out_dir}/{fileHeader}.viruses.multi.reads2taxid_path.tmp\n",
        ]
        with open(f"{out_dir}/bt_2v_tmp.sh", 'w') as f:
            f.writelines("#!/bin/bash\n")
            f.writelines(bash_commands_2)
        os.system(f"chmod +x {out_dir}/bt_2v_tmp.sh")
        os.system(f"{out_dir}/bt_2v_tmp.sh")
        os.remove(f"{out_dir}/bt_2v_tmp.sh")
        
        logging.info("Running LCA on bowtie2 viruses results ...")
        start_time = time.time()
        taxa_path_LCA(resultsPath=out_dir, fileHeader=fileHeader, reads_type=db_type)
        elapsed_time = time.time() - start_time
        logging.info(f"LCA on bowtie2 viruses results finished in {elapsed_time:.2f} seconds.")

        bash_commands_3 = [
            f"cat <(cut -f1,4,5 {out_dir}/{fileHeader}.viruses.uniq.reads2taxid_path.tmp) {out_dir}/{fileHeader}.viruses.multi_reads_path.tmp > {out_dir}/{fileHeader}.viruses.LCA.out\n",
            f"rm {out_dir}/{fileHeader}.viruses.*.tmp\n",
            f"rm {out_dir}/{fileHeader}.viruses.bowtie2.all.out\n",
        ]
        with open(f"{out_dir}/bt_3v_tmp.sh", 'w') as f:
            f.writelines("#!/bin/bash\n")
            f.writelines(bash_commands_3)
        os.system(f"chmod +x {out_dir}/bt_3v_tmp.sh")
        os.system(f"{out_dir}/bt_3v_tmp.sh")
        os.remove(f"{out_dir}/bt_3v_tmp.sh")

    if db_type=="fungi":
        logging.info("Running Bowtie2 on fungal reads ...")
        # print("INFO: Run Bowtie2 on fungal reads ...")
        bash_commands_1 = [
            f"time {bt_path} -x {db_path} -U {input}  --xeq -S {out_dir}/{fileHeader}.fungi.bowtie2.all.out -a -p {threads} --no-unal --no-head\n",
            f"grep \"AS:\" {out_dir}/{fileHeader}.fungi.bowtie2.all.out | grep \"XS:\" -v | cut -f1,3 > {out_dir}/{fileHeader}.fungi.uniq.tmp\n",
            f"grep \"AS:\" {out_dir}/{fileHeader}.fungi.bowtie2.all.out | grep \"XS:\" | cut -f1,3 > {out_dir}/{fileHeader}.fungi.multi.tmp\n"
        ]
        with open(f"{out_dir}/bt_1f_tmp.sh", 'w') as f:
            f.writelines("#!/bin/bash\n")
            f.writelines(bash_commands_1)
        os.system(f"chmod +x {out_dir}/bt_1f_tmp.sh")
        start_time = time.time()
        os.system(f"{out_dir}/bt_1f_tmp.sh")
        elapsed_time = time.time() - start_time
        logging.info(f"Bowtie2 on fungal reads finished in {elapsed_time:.2f} seconds.")
        os.remove(f"{out_dir}/bt_1f_tmp.sh")
        
        append_taxid(bt_results_path=out_dir, acc2taxid_path=acc2taxid_path, reads_type=db_type, fileHeader=fileHeader)
        
        bash_commands_2 = [
            f"cat {out_dir}/{fileHeader}.fungi.uniq.reads2acc2taxid.tmp | {taxonkit_path} reformat -I3 -t > {out_dir}/{fileHeader}.fungi.uniq.reads2taxid_path.tmp\n",
            f"cat {out_dir}/{fileHeader}.fungi.multi.reads2acc2taxid.tmp | {taxonkit_path} reformat -I3 -t > {out_dir}/{fileHeader}.fungi.multi.reads2taxid_path.tmp\n",
        ]
        with open(f"{out_dir}/bt_2f_tmp.sh", 'w') as f:
            f.writelines("#!/bin/bash\n")
            f.writelines(bash_commands_2)
        os.system(f"chmod +x {out_dir}/bt_2f_tmp.sh")
        os.system(f"{out_dir}/bt_2f_tmp.sh")
        os.remove(f"{out_dir}/bt_2f_tmp.sh")
        
        start_time = time.time()
        taxa_path_LCA(resultsPath=out_dir, fileHeader=fileHeader, reads_type=db_type)
        elapsed_time = time.time() - start_time
        logging.info(f"LCA on bowtie2 fungal results finished in {elapsed_time:.2f} seconds.")

        bash_commands_3 = [
            f"cat <(cut -f1,4,5 {out_dir}/{fileHeader}.fungi.uniq.reads2taxid_path.tmp) {out_dir}/{fileHeader}.fungi.multi_reads_path.tmp > {out_dir}/{fileHeader}.fungi.LCA.out\n",
            f"rm {out_dir}/{fileHeader}.fungi.*.tmp\n",
            f"rm {out_dir}/{fileHeader}.fungi.bowtie2.all.out\n",
        ]
        with open(f"{out_dir}/bt_3f_tmp.sh", 'w') as f:
            f.writelines("#!/bin/bash\n")
            f.writelines(bash_commands_3)
        os.system(f"chmod +x {out_dir}/bt_3f_tmp.sh")
        os.system(f"{out_dir}/bt_3f_tmp.sh")
        os.remove(f"{out_dir}/bt_3f_tmp.sh")

def run_RGI(input_fq, out_dir, CARD_db_path, tool_path, fileHeader, threads=8):
    os.makedirs(out_dir, exist_ok=True)
    logging.info(" Running RGI on ...")
    # print("INFO: Run RGI on ...")
    bash_commands = [
        f"exec 2>&1\n",
        f"source {envs.CONDA_PATH}/bin/activate {envs.RTTAP_ENV_NAME}\n",
        f"cd {out_dir}\n",
        f"if [ ! -d ./localDB ];then\n",
        f"    rm -r ./localDB\n",
        f"fi\n",
        f"{tool_path} clean --local\n",
        f"{tool_path} load --card_json {CARD_db_path}/card.json --local\n",
        f"{tool_path} card_annotation -i {CARD_db_path}/card.json > card_annotation.log\n",
        f"{tool_path} load -i {CARD_db_path}/card.json --card_annotation card_database_v3.2.7.fasta --local\n",
        f"{tool_path} bwt -1 {input_fq} -n {threads} -o {out_dir}/{fileHeader}.rgi.out -a kma --local --clean 2>&1\n",
        f"rm {out_dir}/{fileHeader}.*.bam\n",
        f"rm {out_dir}/{fileHeader}.*.bam.bai\n",
        f"{tool_path} clean --local\n",
        f"echo \"INFO: RGI finished!\"\n",
    ]
    with open(f"{out_dir}/run_RGI_tmp.sh", 'w') as f:
        f.writelines("#!/bin/bash\n")
        f.writelines(bash_commands)
    os.system(f"chmod +x {out_dir}/run_RGI_tmp.sh")
    start_time = time.time()
    os.system(f"{out_dir}/run_RGI_tmp.sh")
    elapsed_time = time.time() - start_time
    logging.info(f"RGI finished in {elapsed_time:.2f} seconds. Output to {out_dir}")
    os.remove(f"{out_dir}/run_RGI_tmp.sh")

def generateBowtie2Report(fileHeader, bowtie2VirusesResultsPath, bowtie2FungiResultsPath, kraken2_report, RPM_threshold=1, type="viruses"):
    out = pd.read_table("{}/{}.{}.LCA.out".format(bowtie2VirusesResultsPath if type=="viruses" else bowtie2FungiResultsPath, fileHeader, type), 
                                header=None, 
                                names=["readID", "taxaPath","taxaIDPath"])
    k2_report = pd.read_table(kraken2_report, sep='\t',header=None,names=["percent","reads","reads_direct","rank","taxid","taxa"])
    total_reads_count = k2_report[k2_report["taxa"].isin(["root","unclassified"])]["reads"].sum(axis=0)
    tmp = pd.DataFrame({})
    report_tmp = pd.DataFrame({})
    ranks = ["SuperKingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
    if out.shape[0]==0:
        report = pd.DataFrame({
            "rank":[],
            "taxa":[],
            "taxaID":[],
            "reads":[],
            "RPM":[]
        })
    else:
        for idx, rank in enumerate(ranks):
            if idx <= out["taxaPath"].str.split(";", expand=True).shape[1]:
                tmp["taxa"] = out["taxaPath"].str.split(";", expand=True)[idx].replace("","No Name")
                tmp["taxaID"] = out["taxaIDPath"].str.split(";", expand=True)[idx].replace("","No ID")
                tmp["rank"] = ranks[idx]
                tmp["reads"] = 1
            
            report_tmp = pd.concat([report_tmp, tmp], axis=0).reset_index(drop=True)

        report_tmp = report_tmp.groupby(["rank", "taxa", "taxaID"]).agg({"reads":"sum"})

        report = pd.DataFrame({})
        for idx, rank in enumerate(ranks):
            report = pd.concat([report, report_tmp[report_tmp.index.get_level_values("rank")==rank].sort_values("reads", ascending=False)], axis=0)
        report["RPM"] = report["reads"]*1000000/total_reads_count
        report = report.reset_index()
        report = report[report["RPM"]>=RPM_threshold]
    report.to_csv("{}/{}.{}.report".format(bowtie2VirusesResultsPath if type=="viruses" else bowtie2FungiResultsPath, fileHeader, type), sep='\t', index=None)

def generateMetaphlan4Report(fileHeader, metaphlan4ResultsPath, kraken2_report, RPM_threshold=10):
    # detect header line
    metaphlan4_out_path = "{}/{}.metaphlan4.out".format(metaphlan4ResultsPath, fileHeader)
    def find_metaphlan_header_line(filepath):
        with open(filepath, 'r') as f:
            for i, line in enumerate(f):
                if "#clade_name" in line:
                    return i
        raise ValueError("No header line with '#clade_name' found.")

    header_line = find_metaphlan_header_line(metaphlan4_out_path)
    out = pd.read_table(metaphlan4_out_path, header=header_line, sep='\t').rename({"#clade_name":"taxaPath", "estimated_number_of_reads_from_the_clade":"reads", "clade_taxid":"taxaIDPath"}, axis=1)
    k2_report = pd.read_table(kraken2_report, sep='\t',header=None,names=["percent","reads","reads_direct","rank","taxid","taxa"])
    total_reads_count = k2_report[k2_report["taxa"].isin(["root","unclassified"])]["reads"].sum(axis=0)
    
    # remove unclassified
    if out.shape[0]==0 or out.shape[0]==1:
        report = pd.DataFrame({
            "rank":[],
            "taxa":[],
            "taxaID":[],
            "reads":[],
            "RPM":[]
        })
        logging.warning("No bacteria found in MetaPhlAn4 results.")
        report.to_csv("{}/{}.bacteria.report".format(metaphlan4ResultsPath, fileHeader), sep='\t', index=None)
        return
    else:
        out = out.iloc[1:,:].reset_index(drop=True)
        # pre-process
        out["taxaPath"] = out["taxaPath"].str.replace("|", ";", case=True, regex=False)
        out["taxaIDPath"] = out["taxaIDPath"].str.replace("|", ";", case=True, regex=False)

        report = pd.DataFrame({})
        report_tmp = pd.DataFrame({})
        
        taxaPath_tmp = out["taxaPath"].str.split(";", expand=True)
        taxaIDPath_tmp = out["taxaIDPath"].str.split(";", expand=True)
        taxa = []
        taxaID = []
        for idx in range(out.shape[0]):
            taxa.append(taxaPath_tmp.iloc[idx,:][taxaPath_tmp.iloc[idx,:].T.last_valid_index()])
            taxaID.append(taxaIDPath_tmp.iloc[idx,:][taxaIDPath_tmp.iloc[idx,:].T.last_valid_index()])
        report_tmp["taxa"] = pd.Series(taxa).str.split("__", expand=True)[1]
        report_tmp["taxaID"] = taxaID
        report_tmp["rank"] = pd.Series(taxa).str.split("__", expand=True)[0]
        report_tmp["reads"] = out["reads"]
            
        ranks = ["SuperKingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain"]
        for idx, item in enumerate(["k","p","c","o","f","g","s","t"]):
            report_tmp["rank"] = report_tmp["rank"].replace(item, ranks[idx], regex=False)
        report_tmp = report_tmp.groupby(["rank", "taxa", "taxaID"]).agg({"reads":"sum"})

        for idx, rank in enumerate(ranks[:7]):
            report = pd.concat([report, report_tmp[report_tmp.index.get_level_values("rank")==rank].sort_values("reads", ascending=False)], axis=0)
        
        report = report.reset_index()
        report["taxa"] = report["taxa"].str.replace('_', ' ')
        report["RPM"] = report["reads"]*1000000/total_reads_count
        report = report[report["RPM"]>=RPM_threshold]
        report.to_csv("{}/{}.bacteria.report".format(metaphlan4ResultsPath, fileHeader), sep='\t', index=None)
        return

def generateOveralReport(fileHeader, viral_report_path, bacterial_report_path, fungal_report_path, out_path):

    # bt_viruses = pd.read_table("{}/{}.{}.report".format(bowtie2VirusesResultsPath, fileHeader, "viruses"), sep='\t', header=0)
    # mpa_bacteria = pd.read_table("{}/{}.metaphlan4.report".format(metaphlan4ResultsPath, fileHeader), sep='\t', header=0)
    # bt_fungi = pd.read_table("{}/{}.{}.report".format(bowtie2FungiResultsPath, fileHeader, "fungi"), sep='\t', header=0)
    
    bt_viruses = pd.read_table(viral_report_path, sep='\t', header=0)
    mpa_bacteria = pd.read_table(bacterial_report_path, sep='\t', header=0)
    bt_fungi = pd.read_table(fungal_report_path, sep='\t', header=0)

    ranks = ["SuperKingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]

    report = pd.DataFrame({})
    
    for idx, rank in enumerate(ranks):
        report = pd.concat([report, bt_viruses[bt_viruses["rank"]==rank]], axis=0)
        report = pd.concat([report, mpa_bacteria[mpa_bacteria["rank"]==rank]], axis=0)
        report = pd.concat([report, bt_fungi[bt_fungi["rank"]==rank]], axis=0)

    # report.to_csv(os.path.join(out_dir, f"{fileHeader}.RTTAP.report"), sep='\t', index=None)
    report.to_csv(out_path, sep='\t', index=None)

def run_VirStrain(viral_fq_gz, viral_report, virstrain_list, out_dir, tool_path, db_path, fileHeader):
    os.makedirs(out_dir, exist_ok=True)
    logging.info(f"Run VirStrain on {fileHeader} ...")
    # print(f"INFO: Run VirStrain on {fileHeader} ...")
    report = pd.read_table(viral_report, header=0, index_col=None, sep='\t')
    virstrain_list = pd.read_table(virstrain_list, sep='\t', header=0).set_index(["DB_name"])

    db_list = []
    s_name_list = report["taxa"].str.lower()
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

    # run virstrain
    if db_list.shape[0]==0:
        logging.warning("No strain found in VirStrain DB.")
        # print("INFO: No strain found in VirStrain DB.")
        with open(os.path.join(out_dir, f"{fileHeader}.virstrain.report"), 'w') as f:
            f.writelines(["\n"])
    else:
        bash_commands = []
        bash_commands.append(f"gzip -d -c {viral_fq_gz} > {os.path.join(out_dir, f'{fileHeader}.viruses.fq')}\n")
        for db_name in db_list["DB_name"].tolist():
            bash_commands.append(f"python {tool_path} -i {os.path.join(out_dir, f'{fileHeader}.viruses.fq')} -d {os.path.join(db_path, db_name)} -o {os.path.join(out_dir, db_name)}\n")
        bash_commands.append(f"rm {os.path.join(out_dir, f'{fileHeader}.viruses.fq')}\n")
        with open(f"{out_dir}/run_virstrain_tmp.sh", 'w') as f:
            f.writelines("#!/bin/bash\n")
            f.writelines(bash_commands)
        os.system(f"chmod +x {out_dir}/run_virstrain_tmp.sh")
        start_time = time.time()
        os.system(f"{out_dir}/run_virstrain_tmp.sh")
        elapsed_time = time.time() - start_time
        logging.info(f"VirStrain finished in {elapsed_time:.2f} seconds. Output to {out_dir}")
        os.remove(f"{out_dir}/run_virstrain_tmp.sh")
        
        if os.path.exists(os.path.join(out_dir, f"{fileHeader}.virstrain.report")):
            os.remove(os.path.join(out_dir, f"{fileHeader}.virstrain.report"))
        bash_commands_1 = []
        for db_name in db_list:
            if os.path.exists(os.path.join(out_dir, fileHeader, db_name, "VirStrain_report.txt")):
                bash_commands_1.append(f"grep '>>Most possible' -A1 {os.path.join(out_dir, fileHeader, db_name, 'VirStrain_report.txt')} | sed '1d' | sed 's/^[ \\t]*>//' | cut -f1,2,3,4,5,6,7,10 >> {os.path.join(out_dir, f'{fileHeader}.virstrain.report')}\n")
        with open(f"{out_dir}/gen_virstrain_report_tmp.sh", 'w') as f:
            f.writelines("#!/bin/bash\n")
            f.writelines(bash_commands_1)
        os.system(f"chmod +x {out_dir}/gen_virstrain_report_tmp.sh")
        os.system(f"{out_dir}/gen_virstrain_report_tmp.sh")
        os.remove(f"{out_dir}/gen_virstrain_report_tmp.sh")
        
        if not os.path.exists(os.path.join(out_dir, f'{fileHeader}.virstrain.report')):
            with open(os.path.join(out_dir, f"{fileHeader}.virstrain.report"), 'w') as f:
                f.writelines(["\n"])

def bacterial_stage(out_dir, fileHeader, split_script, mpa_path, rgi_path, CARD_db, metaphlan4_db_path, kraken2_report, RPM_threshold=10, threads=8, skip_rgi=False):
    if not os.path.isfile(os.path.join(out_dir,"splited_reads",f"{fileHeader}.bacteria.fq.gz")):
        split_reads(
            kraken_out=os.path.join(out_dir,"Kraken2_results", f"{fileHeader}.nohost.kraken2ntmicrodb.out"), 
            kraken_report=os.path.join(out_dir,"Kraken2_results", f"{fileHeader}.nohost.kraken2ntmicrodb.report_official"), 
            fq_gz=os.path.join(out_dir,"no_host", f"{fileHeader}.nohost.fq.gz"), 
            out_dir=os.path.join(out_dir,"splited_reads"), 
            fileHeader=fileHeader, reads_type="bacteria", 
            script_path=split_script
        )
    p_list = []
    if not os.path.isfile(os.path.join(out_dir, "Bacteria_results", f"{fileHeader}.metaphlan4.out")):
        p_list.append(
            multiprocessing.Process(
                target=run_Metaphlan4, 
                kwargs={
                    "input": os.path.join(out_dir, "splited_reads", f"{fileHeader}.bacteria.fq.gz"),
                    "out_dir": os.path.join(out_dir, "Bacteria_results"),
                    "tool_path": mpa_path,
                    "fileHeader": fileHeader,
                    "metaphlan4_db_path": metaphlan4_db_path,
                    "threads": threads
                }
            )
        )
    # rgi process
    if not skip_rgi and not os.path.isfile(os.path.join(out_dir, "ARG",f"{fileHeader}.rgi.out.gene_mapping_data.txt")):
        p_list.append(
            multiprocessing.Process(
                target=run_RGI, 
                kwargs={
                    "input_fq": os.path.join(out_dir, "splited_reads", f"{fileHeader}.bacteria.fq.gz"),
                    "out_dir": os.path.join(out_dir, "ARG"),
                    "CARD_db_path": CARD_db,
                    "tool_path": rgi_path,
                    "fileHeader": fileHeader,
                    "threads": threads
                }
            )
        )
    if len(p_list)>0:
        for p in p_list:
            p.start()
        for p in p_list:
            p.join()
    if not os.path.isfile(os.path.join(out_dir, "Bacteria_results", f"{fileHeader}.bacteria.report")):
        generateMetaphlan4Report(
            fileHeader=fileHeader,
            metaphlan4ResultsPath=os.path.join(out_dir, "Bacteria_results"),
            kraken2_report=kraken2_report, RPM_threshold=RPM_threshold
        )

def viral_stage(out_dir, fileHeader, split_script, taxonkit_path, bt_path, #bt_viral_db_path, 
                virstrain_path, virstrain_db_path, virstrain_db_list, acc2taxid_dir,
                kraken2_report, RPM_threshold=1, threads=8, skip_virstrain=False):
    if not os.path.isfile(os.path.join(out_dir,"splited_reads",f"{fileHeader}.viruses.fq.gz")):
        split_reads(
            kraken_out=os.path.join(out_dir,"Kraken2_results", f"{fileHeader}.nohost.kraken2ntmicrodb.out"), 
            kraken_report=os.path.join(out_dir,"Kraken2_results", f"{fileHeader}.nohost.kraken2ntmicrodb.report_official"), 
            fq_gz=os.path.join(out_dir,"no_host", f"{fileHeader}.nohost.fq.gz"), 
            out_dir=os.path.join(out_dir,"splited_reads"), 
            fileHeader=fileHeader, reads_type="viruses", 
            script_path=split_script
        )
    
    if not os.path.isfile(os.path.join(out_dir,"Viruses_results",f"{fileHeader}.viruses.LCA.out")):
        run_bowtie2(
            input=os.path.join(out_dir,"splited_reads", f"{fileHeader}.viruses.fq.gz"),
            out_dir=os.path.join(out_dir,"Viruses_results"),
            bt_path=bt_path,
            taxonkit_path=taxonkit_path,
            # db_path=bt_viral_db_path,
            fileHeader=fileHeader,
            db_type="viruses",
            acc2taxid_dir=acc2taxid_dir,
            threads=threads
        )
    if not os.path.isfile(os.path.join(out_dir,"Viruses_results", f"{fileHeader}.viruses.report")):
        generateBowtie2Report(
            fileHeader=fileHeader, 
            bowtie2FungiResultsPath=os.path.join(out_dir,"Fungi_results"),
            bowtie2VirusesResultsPath=os.path.join(out_dir,"Viruses_results"),
            kraken2_report=kraken2_report, RPM_threshold=RPM_threshold,
            type="viruses"
        )
    # virstrain
    if not skip_virstrain and not os.path.isfile(os.path.join(out_dir,"Viral_strains", f"{fileHeader}.virstrain.report")):
        run_VirStrain(
            viral_fq_gz=os.path.join(out_dir,"splited_reads", f"{fileHeader}.viruses.fq.gz"),
            viral_report=os.path.join(out_dir, "Viruses_results", f"{fileHeader}.viruses.report"),
            virstrain_list=virstrain_db_list,
            out_dir=os.path.join(out_dir,"Viral_strains"),
            tool_path=virstrain_path,
            db_path=virstrain_db_path,
            fileHeader=fileHeader
        )
        
    return

def fungal_stage(out_dir, fileHeader, split_script, taxonkit_path, bt_path, #bt_fungal_db_path, 
                 acc2taxid_dir, kraken2_report, RPM_threshold=10, threads=8):
    if not os.path.isfile(os.path.join(out_dir,"splited_reads",f"{fileHeader}.fungi.fq.gz")):
        split_reads(
            kraken_out=os.path.join(out_dir,"Kraken2_results", f"{fileHeader}.nohost.kraken2ntmicrodb.out"), 
            kraken_report=os.path.join(out_dir,"Kraken2_results", f"{fileHeader}.nohost.kraken2ntmicrodb.report_official"), 
            fq_gz=os.path.join(out_dir,"no_host", f"{fileHeader}.nohost.fq.gz"), 
            out_dir=os.path.join(out_dir,"splited_reads"), 
            fileHeader=fileHeader, reads_type="fungi", 
            script_path=split_script
        )
    if not os.path.isfile(os.path.join(out_dir,"Fungi_results",f"{fileHeader}.fungi.LCA.out")):
        run_bowtie2(
            input=os.path.join(out_dir,"splited_reads", f"{fileHeader}.fungi.fq.gz"),
            out_dir=os.path.join(out_dir,"Fungi_results"),
            bt_path=bt_path,
            taxonkit_path=taxonkit_path,
            # db_path=bt_fungal_db_path,
            fileHeader=fileHeader,
            db_type="fungi",
            acc2taxid_dir=acc2taxid_dir,
            threads=threads
        )
    if not os.path.isfile(os.path.join(out_dir,"Fungi_results",f"{fileHeader}.fungi.report")):
        generateBowtie2Report(
            fileHeader=fileHeader, 
            bowtie2FungiResultsPath=os.path.join(out_dir,"Fungi_results"),
            bowtie2VirusesResultsPath=os.path.join(out_dir,"Viruses_results"),
            kraken2_report=kraken2_report, RPM_threshold=RPM_threshold,
            type="fungi"
        )

    return

def run_pipeline(input_1, out_dir, RPM_b, RPM_v, RPM_f, input_2=None, env_variables={}, threads=8):
    
    fileHeader = os.path.basename(out_dir)
    logging.info(f"RTTAP on analysing {fileHeader} started.")
    # print(f"INFO: RTTAP on analysing {fileHeader} started.")

    fastp_path = env_variables.get("fastp_path")
    if not fastp_path:
        raise ValueError("fastp_path is not set in env_variables.")
    bowtie2_path = env_variables.get("bowtie2_path")
    if not bowtie2_path:
        raise ValueError("bowtie2_path is not set in env_variables.")
    kraken2_path = env_variables.get("kraken2_path")
    if not kraken2_path:
        raise ValueError("kraken2_path is not set in env_variables.")
    mpa_path = env_variables.get("mpa_path")
    if not mpa_path:
        raise ValueError("mpa_path is not set in env_variables.")
    rgi_path = env_variables.get("rgi_path")
    if not rgi_path:
        raise ValueError("rgi_path is not set in env_variables.")
    virstrain_path = env_variables.get("virstrain_path")
    if not virstrain_path:
        raise ValueError("virstrain_path is not set in env_variables.")
    taxonkit_path = env_variables.get("taxonkit_path")
    if not taxonkit_path:
        raise ValueError("taxonkit_path is not set in env_variables.")
    split_script = env_variables.get("split_script")
    if not split_script:
        raise ValueError("split_script is not set in env_variables.")
    metaphlan4_db_path = env_variables.get("metaphlan4_db_path")
    if not metaphlan4_db_path:
        raise ValueError("metaphlan4_db_path is not set in env_variables.")
    CARD_db = env_variables.get("CARD_db")
    if not CARD_db:
        raise ValueError("CARD_db is not set in env_variables.")
    bt_viral_db_path = env_variables.get("bt_viral_db_path")
    if not bt_viral_db_path:
        raise ValueError("bt_viral_db_path is not set in env_variables.")
    bt_fungal_db_path = env_variables.get("bt_fungal_db_path")
    if not bt_fungal_db_path:
        raise ValueError("bt_fungal_db_path is not set in env_variables.")
    virstrain_db_path = env_variables.get("virstrain_db_path")
    if not virstrain_db_path:
        raise ValueError("virstrain_db_path is not set in env_variables.")
    virstrain_db_list = env_variables.get("virstrain_db_list")
    if not virstrain_db_list:
        raise ValueError("virstrain_db_list is not set in env_variables.")
    nt_microbial = env_variables.get("nt_microbial")
    if not nt_microbial:
        raise ValueError("nt_microbial is not set in env_variables.")
    ref_human_db = env_variables.get("ref_human_db")
    if not ref_human_db:
        raise ValueError("ref_human_db is not set in env_variables.")
    acc2taxid_dir = env_variables.get("acc2taxid_dir")
    if not acc2taxid_dir:
        raise ValueError("acc2taxid_dir is not set in env_variables.")

    if input_2 is not None:
        run_fastp(
            input_1=input_1, 
            input_2=input_2, 
            out_dir=os.path.join(out_dir,"QC"), 
            tool_path=fastp_path, 
            fileHeader=fileHeader, 
            threads=threads,
            min_len=15
        )
        
        remove_host(
            input_1=os.path.join(out_dir,"QC",f"{fileHeader}_1.clean.fq.gz"), 
            input_2=os.path.join(out_dir,"QC",f"{fileHeader}_2.clean.fq.gz"), 
            out_dir=os.path.join(out_dir,"no_host"), 
            tool_path=bowtie2_path, 
            db_path=ref_human_db,
            fileHeader=fileHeader, 
            threads=threads, 
            min_len=15,
        )
    else:
        run_fastp(
            input_1=input_1, 
            out_dir=os.path.join(out_dir,"QC"), 
            tool_path=fastp_path, 
            fileHeader=fileHeader, 
            threads=threads,
            min_len=15
        )
        remove_host(
            input_1=os.path.join(out_dir,"QC",f"{fileHeader}.clean.fq.gz"), 
            out_dir=os.path.join(out_dir,"no_host"), 
            tool_path=bowtie2_path, db_path=ref_human_db,
            fileHeader=fileHeader, threads=threads, min_len=15,
        )

    run_Kraken2(
        input=os.path.join(out_dir,"no_host", f"{fileHeader}.nohost.fq.gz"), 
        out_dir=os.path.join(out_dir,"Kraken2_results"), 
        fileHeader=fileHeader, tool_path=kraken2_path,
        db_path=nt_microbial, threads=threads
    )

    p1 = multiprocessing.Process(target=bacterial_stage, 
        kwargs={
            "out_dir": out_dir, 
            "fileHeader": fileHeader, 
            "split_script": split_script, 
            "mpa_path": mpa_path, 
            "rgi_path": rgi_path, 
            "metaphlan4_db_path": metaphlan4_db_path, 
            "CARD_db": CARD_db, 
            "kraken2_report": os.path.join(out_dir,"Kraken2_results", f"{fileHeader}.nohost.kraken2ntmicrodb.report_official"),
            "RPM_threshold": RPM_b,
            "threads": threads
        }
    )

    p2 = multiprocessing.Process(target=viral_stage, 
        kwargs={
            "out_dir": out_dir, 
            "fileHeader": fileHeader, 
            "split_script": split_script, 
            "taxonkit_path": taxonkit_path, 
            "bt_path": bowtie2_path, 
            "bt_viral_db_path": bt_viral_db_path, 
            "virstrain_path": virstrain_path, 
            "virstrain_db_path": virstrain_db_path, 
            "virstrain_db_list": virstrain_db_list, 
            "acc2taxid_path": acc2taxid_dir,
            "kraken2_report": os.path.join(out_dir,"Kraken2_results", f"{fileHeader}.nohost.kraken2ntmicrodb.report_official"),
            "RPM_threshold": RPM_v,
            "threads": threads
        }
    )

    p3 = multiprocessing.Process(target=fungal_stage, 
        kwargs={
            "out_dir": out_dir, 
            "fileHeader": fileHeader, 
            "split_script": split_script, 
            "taxonkit_path": taxonkit_path, 
            "bt_path": bowtie2_path, 
            "bt_fungal_db_path": bt_fungal_db_path, 
            "acc2taxid_path": acc2taxid_dir,
            "kraken2_report": os.path.join(out_dir,"Kraken2_results", f"{fileHeader}.nohost.kraken2ntmicrodb.report_official"),
            "RPM_threshold": RPM_f,
            "threads": threads
        }
    )
    
    p1.start()
    p2.start()
    p3.start()
    p1.join()
    p2.join()
    p3.join()
    logging.info(f"RTTAP on analysing {fileHeader} finished.")
    # print(f"INFO: RTTAP on analysing {fileHeader} finished.")


if __name__=="__main__":
    pass