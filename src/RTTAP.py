import argparse
import os
import sys
from utils import functions
import envs
import logging
logging.basicConfig(
    level=logging.DEBUG, 
    format="%(asctime)s [%(levelname)s]: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

# Requirements:
# 1. fastp
# 2. bowtie2
# 3. Kraken2
# 4. Metaphlan4
# 5. RGI
# 6. VirStrain
# 7. taxonkit

def check_file_exists(file_path, description):
    """Check if file exists. If not, report error and exit."""
    if not os.path.exists(file_path):
        logging.error(f"{description} file not found: {file_path}")
        sys.exit(1)
    return True

def check_step_completed(step_file, description):
    """Check if steps are finished."""
    if os.path.exists(step_file):
        logging.info(f"{description} already completed, skipping...")
        return True
    return False

def run_qc(args):
    file_header = os.path.basename(args.output_dir)
    output_html = os.path.join(args.output_dir, "QC", f"{file_header}.html")
    
    # 检查是否已存在QC输出文件
    if check_step_completed(output_html, "Quality control"):
        return
    
    # 确保输入文件存在
    check_file_exists(args.input_1, "Input file 1")
    if args.input_2 is not None:
        check_file_exists(args.input_2, "Input file 2")
    
    print("INFO: Running quality control...")
    functions.run_fastp(
        input_1=args.input_1, 
        input_2=args.input_2, 
        out_dir=os.path.join(args.output_dir, "QC"), 
        tool_path=envs.FASTP_PATH, 
        fileHeader=file_header, 
        threads=args.threads,
        min_len=15,
        cmd=args.cmd
    )
    
    # 验证QC输出
    check_file_exists(output_html, "QC output")

def run_rrna_removal(args):
    file_header = os.path.basename(args.output_dir)
    output_fq = os.path.join(args.output_dir, "no_rRNA", f"{file_header}.norRNA.fq.gz")
    
    # 检查是否已存在rRNA去除输出文件
    # if True skip
    if check_step_completed(output_fq, "rRNA removal"):
        return
    
    # 检查QC步骤是否完成
    qc_output_html = os.path.join(args.output_dir, "QC", f"{file_header}.html")
    check_file_exists(qc_output_html, "QC output")
    
    logging.info("Removing rRNA reads...")
    # print("INFO: Removing rRNA reads...")
    functions.remove_rRNA(
        input_1=os.path.join(args.output_dir, "QC", f"{file_header}.clean.R1.fq.gz") if args.input_2 
                else os.path.join(args.output_dir, "QC", f"{file_header}.clean.fq.gz"),
        input_2=os.path.join(args.output_dir, "QC", f"{file_header}.clean.R2.fq.gz") if args.input_2 else None,
        out_dir=os.path.join(args.output_dir, "no_rRNA"), 
        tool_path=envs.BOWTIE2_PATH, 
        db_path=envs.REF_HUMAN_DB,
        fileHeader=file_header, 
        threads=args.threads, 
        cmd=args.cmd
    )
    
    # 验证rRNA去除输出
    check_file_exists(output_fq, "rRNA removal output")

def run_kraken2(args):
    file_header = os.path.basename(args.output_dir)
    kraken_out = os.path.join(args.output_dir, "Kraken2_results", f"{file_header}.norRNA.kraken2ntmicrodb.report_official")
    
    # 检查是否已存在Kraken2输出文件
    if check_step_completed(kraken_out, "Kraken2 classification"):
        return
    
    # 检查rRNA去除步骤是否完成
    rrna_output = os.path.join(args.output_dir, "no_rRNA", f"{file_header}.norRNA.fq.gz")
    check_file_exists(rrna_output, "rRNA removal output")
    
    logging.info("Running Kraken2 classification...")
    # print("INFO: Running Kraken2 classification...")
    functions.run_Kraken2(
        input=rrna_output, 
        out_dir=os.path.join(args.output_dir, "Kraken2_results"), 
        fileHeader=file_header, 
        tool_path=envs.KRAKEN2_PATH,
        db_path=envs.NT_MICROBIAL, 
        threads=args.threads,
        cmd=args.cmd,
    )
    
    # 验证Kraken2输出
    check_file_exists(kraken_out, "Kraken2 output")

def run_bacterial_analysis(args):
    file_header = os.path.basename(args.output_dir)
    mpa_report = os.path.join(args.output_dir, "Bacteria_results", f"{file_header}.bacteria.report")
    
    # 检查是否已存在细菌分析结果
    if check_step_completed(mpa_report, "Bacterial analysis"):
        return
    
    # 检查Kraken2步骤是否完成
    kraken_out = os.path.join(args.output_dir, "Kraken2_results", f"{file_header}.norRNA.kraken2ntmicrodb.report_official")
    check_file_exists(kraken_out, "Kraken2 output")
    
    logging.info("Running bacterial analysis...")
    # print("INFO: Running bacterial analysis...")
    functions.bacterial_stage(
        out_dir=args.output_dir, 
        fileHeader=file_header, 
        split_script=envs.SPLIT_SCRIPT, 
        mpa_path=envs.METAPHLAN4_PATH, 
        rgi_path=envs.RGI_PATH, 
        metaphlan4_db_path=envs.METAPHLAN4_DB_PATH, 
        CARD_db=envs.CARD_DB_PATH, 
        kraken2_report=os.path.join(args.output_dir, "Kraken2_results", f"{file_header}.norRNA.kraken2ntmicrodb.report_official"),
        RPM_threshold=args.rpm_b,
        threads=args.threads,
        skip_rgi=args.skip_rgi
    )
    
    # 验证细菌分析输出
    check_file_exists(mpa_report, "Metaphlan4 output")

def run_viral_analysis(args):
    file_header = os.path.basename(args.output_dir)
    viral_report = os.path.join(args.output_dir, "Viruses_results", f"{file_header}.viruses.report")
    LCA_out = os.path.join(args.output_dir, "Viruses_results", f"{file_header}.viruses.LCA.out")
    
    # 检查是否已存在病毒分析结果
    if check_step_completed(viral_report, "Viral analysis"):
        return
    
    # 检查Kraken2步骤是否完成
    kraken_out = os.path.join(args.output_dir, "Kraken2_results", f"{file_header}.norRNA.kraken2ntmicrodb.report_official")
    check_file_exists(kraken_out, "Kraken2 output")
    
    logging.info("Running viral analysis...")
    # print("INFO: Running viral analysis...")
    functions.viral_stage(
        out_dir=args.output_dir, 
        fileHeader=file_header, 
        split_script=envs.SPLIT_SCRIPT, 
        taxonkit_path=envs.TAXONKIT_PATH, 
        bt_path=envs.BOWTIE2_PATH, 
        bt_viral_db_path=envs.BT_VIRAL_DB_PATH, 
        virstrain_path=envs.VIRSTRAIN_PATH, 
        virstrain_db_path=envs.VIRSTRAIN_DB_PATH, 
        virstrain_db_list=envs.VIRSTRAIN_DB_LIST, 
        acc2taxid_path=envs.ACC2TAXID_PATH,
        kraken2_report=os.path.join(args.output_dir, "Kraken2_results", f"{file_header}.norRNA.kraken2ntmicrodb.report_official"),
        RPM_threshold=args.rpm_v,
        threads=args.threads,
        skip_virstrain=args.skip_virstrain
    )
    
    # 验证病毒分析输出
    check_file_exists(viral_report, "Viral taxonomic analysis output")

def run_fungal_analysis(args):
    file_header = os.path.basename(args.output_dir)
    fungal_report = os.path.join(args.output_dir, "Fungi_results", f"{file_header}.fungi.report")
    
    # 检查是否已存在真菌分析结果
    if check_step_completed(fungal_report, "Fungal analysis"):
        return
    
    # 检查Kraken2步骤是否完成
    kraken_out = os.path.join(args.output_dir, "Kraken2_results", f"{file_header}.norRNA.kraken2ntmicrodb.report_official")
    check_file_exists(kraken_out, "Kraken2 output")
    
    logging.info("Running fungal analysis...")
    # print("INFO: Running fungal analysis...")
    functions.fungal_stage(
        out_dir=args.output_dir, 
        fileHeader=file_header, 
        split_script=envs.SPLIT_SCRIPT, 
        taxonkit_path=envs.TAXONKIT_PATH, 
        bt_path=envs.BOWTIE2_PATH, 
        bt_fungal_db_path=envs.BT_FUNGAL_DB_PATH, 
        acc2taxid_path=envs.ACC2TAXID_PATH,
        kraken2_report=os.path.join(args.output_dir, "Kraken2_results", f"{file_header}.norRNA.kraken2ntmicrodb.report_official"),
        RPM_threshold=args.rpm_f,
        threads=args.threads
    )
    
    # 验证真菌分析输出
    check_file_exists(fungal_report, "Fungal analysis output")

def run_generate_report(args, force=True):
    file_header = os.path.basename(args.output_dir)
    kraken2_report=os.path.join(args.output_dir, "Kraken2_results", f"{file_header}.norRNA.kraken2ntmicrodb.report_official")
    viral_report_path = os.path.join(args.output_dir, "Viruses_results", f"{file_header}.viruses.report")
    bacterial_report_path= os.path.join(args.output_dir, "Bacteria_results", f"{file_header}.bacteria.report")
    fungal_report_path = os.path.join(args.output_dir, "Fungi_results", f"{file_header}.fungi.report")
    out_dir = os.path.join(args.output_dir, f"{file_header}.RTTAP.report")
    
    paths = {
        "Viruses report": viral_report_path,
        "Bacteria report": bacterial_report_path,
        "Fungi report": fungal_report_path,
    }
    
    functions.generateBowtie2Report(
        fileHeader=file_header, 
        bowtie2FungiResultsPath=os.path.join(args.output_dir,"Fungi_results"),
        bowtie2VirusesResultsPath=os.path.join(args.output_dir,"Viruses_results"),
        kraken2_report=kraken2_report, RPM_threshold=args.rpm_v,
        type="viruses"
    )
    functions.generateMetaphlan4Report(
        fileHeader=file_header,
        metaphlan4ResultsPath=os.path.join(args.output_dir, "Bacteria_results"),
        kraken2_report=kraken2_report, RPM_threshold=args.rpm_b
    )
    functions.generateBowtie2Report(
        fileHeader=file_header, 
        bowtie2FungiResultsPath=os.path.join(args.output_dir,"Fungi_results"),
        bowtie2VirusesResultsPath=os.path.join(args.output_dir,"Viruses_results"),
        kraken2_report=kraken2_report, RPM_threshold=args.rpm_f,
        type="fungi"
    )
    flag = True
    for key in paths.keys():
        if os.path.isfile(paths[key]):
            continue
        else:
            logging.error(f"{key} of {file_header} does not exist, failed to generate overall report.")
            flag = False
    if flag==True:
        functions.generateOveralReport(file_header, viral_report_path, bacterial_report_path, fungal_report_path, out_dir)
        logging.info(f"Generated overall report to {out_dir}")
        
    return
    
def run_end_to_end(args):
    file_header = os.path.basename(args.output_dir)
    viral_report_path = os.path.join(args.output_dir, "Viruses_results", f"{file_header}.viruses.report")
    bacterial_report_path= os.path.join(args.output_dir, "Bacteria_results", f"{file_header}.bacteria.report")
    fungal_report_path = os.path.join(args.output_dir, "Fungi_results", f"{file_header}.fungi.report")
    out_dir = os.path.join(args.output_dir, f"{file_header}.RTTAP.report")
    logging.info(f"RTTAP end-to-end analysis on {file_header} started.")
    paths = {
        "Viruses report": viral_report_path,
        "Bacteria report": bacterial_report_path,
        "Fungi report": fungal_report_path,
    }
    # print(f"INFO: RTTAP end-to-end analysis on {file_header} started.")
    
    # 按顺序运行所有步骤，每个步骤会自动检查是否已完成
    run_qc(args)
    run_rrna_removal(args)
    run_kraken2(args)
    
    # 并行运行分析步骤
    import multiprocessing
    processes = []
    
    if args.rpm_b >= 0:
        p1 = multiprocessing.Process(target=run_bacterial_analysis, args=(args,))
        processes.append(p1)
    
    if args.rpm_v >= 0:
        p2 = multiprocessing.Process(target=run_viral_analysis, args=(args,))
        processes.append(p2)
    
    if args.rpm_f >= 0:
        p3 = multiprocessing.Process(target=run_fungal_analysis, args=(args,))
        processes.append(p3)
    
    for p in processes:
        p.start()
    
    for p in processes:
        p.join()
    
    # run_generate_report(args)
    flag = True
    for key in paths.keys():
        if os.path.isfile(paths[key]):
            continue
        else:
            logging.error(f"{key} of {file_header} does not exist, failed to generate overall report.")
            flag = False
    if flag==True:
        functions.generateOveralReport(file_header, viral_report_path, bacterial_report_path, fungal_report_path, out_dir)
        logging.info(f"Generated overall report to {out_dir}")
    
    logging.info(f"RTTAP end-to-end analysis on {file_header} finished.")
    # print(f"INFO: RTTAP end-to-end analysis on {file_header} finished.")


def main():
    parser = argparse.ArgumentParser(
        prog="RTTAP",
        description="A Read-based Total-infecTome Analysis Pipeline."
    )
    
    # Create subparsers for different modules
    subparsers  = parser.add_subparsers(
        title="Modules",
        dest="modules",
        description="Modules that proceed different functions.", 
        help="Please specify one and only one of these options.", 
        required=True
    )
    # Common arguments for all commands
    base_parser = argparse.ArgumentParser(add_help=False)
    base_parser.add_argument(
        "-i","--input_1",
        required=True,
        help="The first input file (\".fq\", \".fq.gz\")."
    )
    base_parser.add_argument(
        "-I","--input_2",
        required=False, default=None,
        help="The second input file (\".fq\", \".fq.gz\"; used for paired-end files)."
    )
    base_parser.add_argument(
        "-o","--output_dir",
        required=True,
        help="Output directory for results."
    )
    base_parser.add_argument(
        "-t", "--threads",
        required=False, default=8, type=int,
        help="Threads used for analysis. Default: 8"
    )
    base_parser.add_argument(
        "--rpm_b", default=10, type=int,
        help="The cut-off threshold for bacterial categories. Default: 10"
    )
    base_parser.add_argument(
        "--rpm_v", default=1, type=int,
        help="The cut-off threshold for viral categories. Default: 1"
    )
    base_parser.add_argument(
        "--rpm_f", default=10, type=int,
        help="The cut-off threshold for fungal categories. Default: 10"
    )
    base_parser.add_argument(
        "--skip_rgi", default=False, action="store_true",
        help="Skip the ARG analysis. Default: False"
    )
    base_parser.add_argument(
        "--skip_virstrain", default=False, action="store_true",
        help="Skip the viruses strain analysis. Default: False"
    )
    base_parser.add_argument(
        "--cmd", type=str, default=None,
        help="The string of commands to be added to runnning tools. Default: None. For example, when runnning kraken2, '--cmd \"--confidence 0.5\"' will set the minimum confidence to 0.5 for Kraken2 to classify a sequence. This parameter should be used very carefully as it allows great flexibility while introducing a lot of uncertainty for smoothly running the whole pipeline."
    )
    
    # QC module
    qc_parser = subparsers.add_parser('qc', parents=[base_parser], help='Run quality control only')
    qc_parser.set_defaults(func=run_qc)
    
    # rRNA removal module
    rrna_parser = subparsers.add_parser('rrna_removal', parents=[base_parser], help='Run rRNA removal only')
    rrna_parser.set_defaults(func=run_rrna_removal)
    
    # Kraken2 module
    kraken_parser = subparsers.add_parser('kraken2', parents=[base_parser], help='Run Kraken2 classification only')
    kraken_parser.set_defaults(func=run_kraken2)
    
    # Bacterial analysis module
    bacteria_parser = subparsers.add_parser('bacteria', parents=[base_parser], help='Run bacterial analysis only')
    bacteria_parser.set_defaults(func=run_bacterial_analysis)
    
    # Viral analysis module
    viral_parser = subparsers.add_parser('virus', parents=[base_parser], help='Run viral analysis only')
    viral_parser.set_defaults(func=run_viral_analysis)
    
    # Fungal analysis module
    fungal_parser = subparsers.add_parser('fungi', parents=[base_parser], help='Run fungal analysis only')
    fungal_parser.set_defaults(func=run_fungal_analysis)
    
    # Overall report module
    report_parser = subparsers.add_parser('report', parents=[base_parser], help='Generate overall report for viruses, bacteria, and fungi')
    report_parser.set_defaults(func=run_generate_report)
    
    # End-to-end analysis
    endtoend_parser = subparsers.add_parser('end_to_end', parents=[base_parser], help='Run complete end-to-end analysis')
    endtoend_parser.set_defaults(func=run_end_to_end)
    
    args = parser.parse_args()
    if not hasattr(args, 'func'):
        parser.print_help()
        return
    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)
    # Run the selected function
    args.func(args)

if __name__=="__main__":
    main()