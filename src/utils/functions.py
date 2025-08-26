import os
import pandas as pd
import subprocess as sp
import multiprocessing
import logging
logging.basicConfig(
    level=logging.DEBUG, 
    format="%(asctime)s [%(levelname)s]: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
import time
import envs

def run_fastp(input_1, out_dir, tool_path, fileHeader, threads=8, input_2=None, min_len=15, cmd=None):
    if cmd is None:
        cmd = ''
    os.makedirs(out_dir, exist_ok=True)
    tmp_script = os.path.join(out_dir, 'run_fastp_tmp.sh')
    logging.info("Running Fastp ...")
    # print("INFO: Running Fastp ...")
    if input_2 is not None:
        bash_commands = [
            f"{tool_path} -i {input_1} -I {input_2} -o {out_dir}/{fileHeader}.clean.R1.fq.gz -O {out_dir}/{fileHeader}.clean.R2.fq.gz --json {out_dir}/{fileHeader}.json --html {out_dir}/{fileHeader}.html --thread {threads} --length_required {min_len} -D {cmd}\n",
        ]
    else:
        bash_commands = [
            f"{tool_path} -i {input_1} -o {out_dir}/{fileHeader}.clean.fq.gz --json {out_dir}/{fileHeader}.json --html {out_dir}/{fileHeader}.html --thread {threads} --length_required {min_len} -D\n",
        ]

    with open(tmp_script, 'w') as f:
        f.writelines("#!/bin/bash\n")
        f.writelines(bash_commands)
    os.system(f"chmod +x {tmp_script}")
    start_time = time.time()
    os.system(f"{tmp_script}")
    elapsed_time = time.time() - start_time
    logging.info(f"Fastp finished in {elapsed_time:.2f} seconds. Output to {out_dir}")
    os.remove(tmp_script)

def remove_rRNA(input_1, out_dir, tool_path, db_path, fileHeader, threads=8, input_2=None, cmd=None):
    if cmd is None:
        cmd = ''
    os.makedirs(out_dir, exist_ok=True)
    logging.info("Removing human rRNA reads ...")
    # print("INFO: Removing human rRNA reads ...")
    if input_2 is not None:
        logging.info(f"Input:\t{input_1}\n\t\t{input_2}")
        un_conc_path = os.path.join(out_dir, f"{fileHeader}.norRNA.fq")
        bt2_out_path = os.path.join(out_dir, f"{fileHeader}.norRNA.bowtie2.out")
        bash_commands = [
            f"{tool_path} -x {db_path} -1 {input_1} -2 {input_2} --un-conc {un_conc_path} -S {bt2_out_path} --very-sensitive-local --no-unal -I 1 -p {threads} {cmd}\n",
            f"mv {os.path.join(out_dir, f'{fileHeader}.norRNA.1.fq')} {os.path.join(out_dir, f'{fileHeader}.norRNA.R1.fq')}\n",
            f"mv {os.path.join(out_dir, f'{fileHeader}.norRNA.2.fq')} {os.path.join(out_dir, f'{fileHeader}.norRNA.R2.fq')}\n",
            f"cat {os.path.join(out_dir, f'{fileHeader}.norRNA.R1.fq')} {os.path.join(out_dir, f'{fileHeader}.norRNA.R2.fq')} | sed 's/ /_/g' | gzip -f -c -1 > {os.path.join(out_dir, f'{fileHeader}.norRNA.fq.gz')}\n",
            f"rm {os.path.join(out_dir, f'{fileHeader}.norRNA.R1.fq')} {os.path.join(out_dir, f'{fileHeader}.norRNA.R2.fq')}\n",
            f"rm {bt2_out_path}\n",
        ]
    else:
        logging.info(f"Input:\t{input_1}")
        bash_commands = [
            f"{tool_path} -x {db_path} -U {input_1} --un {os.path.join(out_dir, f'{fileHeader}.norRNA.fq')} -S {os.path.join(out_dir, f'{fileHeader}.norRNA.bowtie2.out')} --very-sensitive-local --no-unal -I 1 -p {threads} {'' if cmd is None else cmd}\n",
            f"cat {os.path.join(out_dir, f'{fileHeader}.norRNA.fq')} | sed 's/ /_/g' | gzip -f -c -1 > {os.path.join(out_dir, f'{fileHeader}.norRNA.fq.gz')}\n",
            f"rm {os.path.join(out_dir, f'{fileHeader}.norRNA.fq')}\n",
            f"rm {os.path.join(out_dir, f'{fileHeader}.norRNA.bowtie2.out')}\n",
        ]
    with open(f"{out_dir}/remove_rRNA_tmp.sh", 'w') as f:
        f.writelines("#!/bin/bash\n")
        f.writelines(bash_commands)
    os.system(f"chmod +x {out_dir}/remove_rRNA_tmp.sh")
    start_time = time.time()
    os.system(f"{out_dir}/remove_rRNA_tmp.sh")
    elapsed_time = time.time() - start_time
    logging.info(f"Removing human rRNA reads finished in {elapsed_time:.2f} seconds. Output to {out_dir}")
    os.remove(f"{out_dir}/remove_rRNA_tmp.sh")

def run_Kraken2(input, out_dir, fileHeader, tool_path, db_path, threads=8, cmd=None):
    if cmd is None:
        cmd = ''
    os.makedirs(out_dir, exist_ok=True)
    logging.info("Running Kraken2 ...")
    kraken_out = os.path.join(out_dir, f"{fileHeader}.norRNA.kraken2ntmicrodb.out")
    kraken_report = os.path.join(out_dir, f"{fileHeader}.norRNA.kraken2ntmicrodb.report_official")
    bash_commands = [
        f"time {tool_path} --db {db_path} --threads {threads} --output {kraken_out} {input} --report {kraken_report} --gzip-compressed {cmd}\n",
        f"echo 'INFO: Kraken2 finished.'\n"
    ]
    with open(f"{out_dir}/run_Kraken2_tmp.sh", 'w') as f:
        f.writelines("#!/bin/bash\n")
        f.writelines(bash_commands)
    os.system(f"chmod +x {out_dir}/run_Kraken2_tmp.sh")
    start_time = time.time()
    os.system(f"{out_dir}/run_Kraken2_tmp.sh")
    elapsed_time = time.time() - start_time
    logging.info(f"Kraken2 finished in {elapsed_time:.2f} seconds. Output to {out_dir}")
    os.remove(f"{out_dir}/run_Kraken2_tmp.sh")

def split_reads(kraken_out, kraken_report, fq_gz, out_dir, fileHeader, reads_type, script_path="./extract_kraken_reads_nostdout.py"):
    os.makedirs(out_dir, exist_ok=True)
    logging.info(f"Splitting {reads_type} reads...")
    fq_bac = os.path.join(out_dir, f"{fileHeader}.bacteria.fq")
    fq_vir = os.path.join(out_dir, f"{fileHeader}.viruses.fq")
    fq_fun = os.path.join(out_dir, f"{fileHeader}.fungi.fq")
    if reads_type=="bacteria":
        bash_commands = [
            f"python {script_path} -k {kraken_out} -s {fq_gz} -r {kraken_report} -t 2 -o {fq_bac} --include-children --fastq-output\n",
            f"gzip --fast {fq_bac} -f\n",
            # f"echo 'INFO: Splitting {reads_type} reads finished.'\n"
        ]
    if reads_type=="viruses":
        bash_commands = [
            f"python {script_path} -k {kraken_out} -s {fq_gz} -r {kraken_report} -t 10239 -o {fq_vir} --include-children --fastq-output\n",
            f"gzip --fast {fq_vir} -f\n",
            # f"echo 'INFO: Splitting {reads_type} reads finished.'\n"
        ]
    if reads_type=="fungi":
        bash_commands = [
            f"python {script_path} -k {kraken_out} -s {fq_gz} -r {kraken_report} -t 4751 -o {fq_fun} --include-children --fastq-output\n",
            f"gzip --fast {fq_fun} -f\n",
            # f"echo 'INFO: Splitting {reads_type} reads finished.'\n"
        ]

    tmp_script = os.path.join(out_dir, f"split_reads_{reads_type}_tmp.sh")
    with open(tmp_script, 'w') as f:
        f.writelines("#!/bin/bash\n")
        f.writelines(bash_commands)
    os.system(f"chmod +x {tmp_script}")
    start_time = time.time()
    os.system(tmp_script)
    elapsed_time = time.time() - start_time
    logging.info(f"Splitting {reads_type} reads finished in {elapsed_time:.2f} seconds. Output to {out_dir}")
    os.remove(tmp_script)

def run_Metaphlan4(input, out_dir, tool_path, fileHeader, metaphlan4_db_path, threads=8):
    os.makedirs(out_dir, exist_ok=True)
    mpa_bt_out = os.path.join(out_dir, f"{fileHeader}.metaphlan4.bowtie2.out")
    mpa_out = os.path.join(out_dir, f"{fileHeader}.metaphlan4.out")
    tmp_script = os.path.join(out_dir, 'run_Metaphlan4_tmp.sh')
    logging.info("Running MetaPhlAn4 ...")
    
    if os.path.exists(mpa_bt_out):
        os.remove(mpa_bt_out)
    bash_commands = [
        f"source {os.path.join(envs.CONDA_PATH, 'bin', 'activate')} {envs.RTTAP_ENV_NAME}\n",
        f"{tool_path} {input} --bowtie2db {os.path.dirname(metaphlan4_db_path)} -x {os.path.basename(metaphlan4_db_path)} -t rel_ab_w_read_stats -o {mpa_out} --nproc {threads} --input_type fastq --CAMI_format_output --bowtie2out {mpa_bt_out} --unclassified_estimation\n",
    ]
    with open(tmp_script, 'w') as f:
        f.writelines("#!/bin/bash\n")
        f.writelines(bash_commands)
    os.system(f"chmod +x {tmp_script}")
    start_time = time.time()
    os.system(tmp_script)
    elapsed_time = time.time() - start_time
    logging.info(f"MetaPhlAn4 finished in {elapsed_time:.2f} seconds. Output to {out_dir}")
    os.remove(tmp_script)

def append_taxid(bt_results_path, acc2taxid_path, reads_type, fileHeader):
    if reads_type in ["viruses", "fungi", "archaea"]:
        resultsPath = bt_results_path
        if reads_type == "viruses":
            acc2taxid = pd.read_table(f"{acc2taxid_path}/viruses.acc2taxid.txt", header=None, names=["acc.ver", "taxid"]).astype({"taxid":"str"})
        if reads_type == "fungi":
            acc2taxid = pd.read_table(f"{acc2taxid_path}/EuPathDB46.seq2taxid.txt", header=None, names=["acc.ver", "taxid"]).astype({"taxid":"str"})
    
        uniqread2acc = pd.read_table(f"{resultsPath}/{fileHeader}.{reads_type}.uniq.tmp", header=None, names=["readId", "acc.ver"])
        multiread2acc = pd.read_table(f"{resultsPath}/{fileHeader}.{reads_type}.multi.tmp", header=None, names=["readId", "acc.ver"])
        
        # add taxid after acc.ver column
        uniqacc2taxid = pd.merge(uniqread2acc, acc2taxid, how="left", on="acc.ver")
        multiacc2taxid = pd.merge(multiread2acc, acc2taxid, how="left", on="acc.ver")
        
        # save temporary file
        uniqacc2taxid.to_csv(f"{resultsPath}/{fileHeader}.{reads_type}.uniq.reads2acc2taxid.tmp", header=None, index=False, sep='\t')
        multiacc2taxid.to_csv(f"{resultsPath}/{fileHeader}.{reads_type}.multi.reads2acc2taxid.tmp", header=None, index=False, sep='\t')

def taxa_path_LCA(resultsPath, fileHeader, reads_type):
    def readLCA(path):
        splittedPath = path.str.split(";", expand=True)
        values = splittedPath.iloc[0, :][splittedPath.apply(pd.Series.nunique) == 1].tolist()
        # change into string and filter out NaN
        values = [str(v) for v in values if pd.notnull(v)]
        return ";".join(values)

    if reads_type in ["viruses", "fungi", "archaea"]:
        read2taxaPath = pd.read_csv(f"{resultsPath}/{fileHeader}.{reads_type}.multi.reads2taxid_path.tmp", header=None, names=["readId", "refId", "taxid", "namePath", "taxaPath"], sep='\t')
        reads_new = read2taxaPath.groupby("readId").agg({"namePath": readLCA, "taxaPath": readLCA}).reset_index()
        reads_new.to_csv(f"{resultsPath}/{fileHeader}.{reads_type}.multi_reads_path.tmp", header=None, index=None, sep='\t')

def run_Bowtie2(input, out_dir, bt_path, taxonkit_path, db_path, fileHeader, db_type, acc2taxid_path, threads=8):
    os.makedirs(out_dir, exist_ok=True)
    
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
            kraken_out=os.path.join(out_dir,"Kraken2_results", f"{fileHeader}.norRNA.kraken2ntmicrodb.out"), 
            kraken_report=os.path.join(out_dir,"Kraken2_results", f"{fileHeader}.norRNA.kraken2ntmicrodb.report_official"), 
            fq_gz=os.path.join(out_dir,"no_rRNA", f"{fileHeader}.norRNA.fq.gz"), 
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

def viral_stage(out_dir, fileHeader, split_script, taxonkit_path, bt_path, bt_viral_db_path, 
                virstrain_path, virstrain_db_path, virstrain_db_list, acc2taxid_path,
                kraken2_report, RPM_threshold=1, threads=8, skip_virstrain=False):
    if not os.path.isfile(os.path.join(out_dir,"splited_reads",f"{fileHeader}.viruses.fq.gz")):
        split_reads(
            kraken_out=os.path.join(out_dir,"Kraken2_results", f"{fileHeader}.norRNA.kraken2ntmicrodb.out"), 
            kraken_report=os.path.join(out_dir,"Kraken2_results", f"{fileHeader}.norRNA.kraken2ntmicrodb.report_official"), 
            fq_gz=os.path.join(out_dir,"no_rRNA", f"{fileHeader}.norRNA.fq.gz"), 
            out_dir=os.path.join(out_dir,"splited_reads"), 
            fileHeader=fileHeader, reads_type="viruses", 
            script_path=split_script
        )
    
    if not os.path.isfile(os.path.join(out_dir,"Viruses_results",f"{fileHeader}.viruses.LCA.out")):
        run_Bowtie2(
            input=os.path.join(out_dir,"splited_reads", f"{fileHeader}.viruses.fq.gz"),
            out_dir=os.path.join(out_dir,"Viruses_results"),
            bt_path=bt_path,
            taxonkit_path=taxonkit_path,
            db_path=bt_viral_db_path,
            fileHeader=fileHeader,
            db_type="viruses",
            acc2taxid_path=acc2taxid_path,
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

def fungal_stage(out_dir, fileHeader, split_script, taxonkit_path, bt_path, bt_fungal_db_path, 
                 acc2taxid_path, kraken2_report, RPM_threshold=10, threads=8):
    if not os.path.isfile(os.path.join(out_dir,"splited_reads",f"{fileHeader}.fungi.fq.gz")):
        split_reads(
            kraken_out=os.path.join(out_dir,"Kraken2_results", f"{fileHeader}.norRNA.kraken2ntmicrodb.out"), 
            kraken_report=os.path.join(out_dir,"Kraken2_results", f"{fileHeader}.norRNA.kraken2ntmicrodb.report_official"), 
            fq_gz=os.path.join(out_dir,"no_rRNA", f"{fileHeader}.norRNA.fq.gz"), 
            out_dir=os.path.join(out_dir,"splited_reads"), 
            fileHeader=fileHeader, reads_type="fungi", 
            script_path=split_script
        )
    if not os.path.isfile(os.path.join(out_dir,"Fungi_results",f"{fileHeader}.fungi.LCA.out")):
        run_Bowtie2(
            input=os.path.join(out_dir,"splited_reads", f"{fileHeader}.fungi.fq.gz"),
            out_dir=os.path.join(out_dir,"Fungi_results"),
            bt_path=bt_path,
            taxonkit_path=taxonkit_path,
            db_path=bt_fungal_db_path,
            fileHeader=fileHeader,
            db_type="fungi",
            acc2taxid_path=acc2taxid_path,
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
    
    fastp_path = env_variables["fastp_path"]
    bowtie2_path = env_variables["bowtie2_path"]
    kraken2_path = env_variables["kraken2_path"]
    mpa_path = env_variables["mpa_path"]
    rgi_path = env_variables["rgi_path"]
    virstrain_path = env_variables["virstrain_path"]
    taxonkit_path = env_variables["taxonkit_path"]
    
    split_script = env_variables["split_script"]
    
    metaphlan4_db_path = env_variables["metaphlan4_db_path"]
    CARD_db = env_variables["CARD_db"]
    bt_viral_db_path = env_variables["bt_viral_db_path"]
    bt_fungal_db_path = env_variables["bt_fungal_db_path"]
    virstrain_db_path = env_variables["virstrain_db_path"]

    virstrain_db_list = env_variables["virstrain_db_list"]

    nt_microbial = env_variables["nt_microbial"]
    ref_human_db = env_variables["ref_human_db"]

    acc2taxid_path = env_variables["acc2taxid_path"]

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
        
        remove_rRNA(
            input_1=os.path.join(out_dir,"QC",f"{fileHeader}_1.clean.fq.gz"), 
            input_2=os.path.join(out_dir,"QC",f"{fileHeader}_2.clean.fq.gz"), 
            out_dir=os.path.join(out_dir,"no_rRNA"), 
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
        remove_rRNA(
            input_1=os.path.join(out_dir,"QC",f"{fileHeader}.clean.fq.gz"), 
            out_dir=os.path.join(out_dir,"no_rRNA"), 
            tool_path=bowtie2_path, db_path=ref_human_db,
            fileHeader=fileHeader, threads=threads, min_len=15,
        )

    run_Kraken2(
        input=os.path.join(out_dir,"no_rRNA", f"{fileHeader}.norRNA.fq.gz"), 
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
            "kraken2_report": os.path.join(out_dir,"Kraken2_results", f"{fileHeader}.norRNA.kraken2ntmicrodb.report_official"),
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
            "acc2taxid_path": acc2taxid_path,
            "kraken2_report": os.path.join(out_dir,"Kraken2_results", f"{fileHeader}.norRNA.kraken2ntmicrodb.report_official"),
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
            "acc2taxid_path": acc2taxid_path,
            "kraken2_report": os.path.join(out_dir,"Kraken2_results", f"{fileHeader}.norRNA.kraken2ntmicrodb.report_official"),
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