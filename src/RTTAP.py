import os
import argparse
import multiprocessing as mp

import logging
logging.basicConfig(
    level=logging.DEBUG, 
    format="%(asctime)s [%(levelname)s]: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

import envs
from utils import workflows 
from utils.api_downstream import generate_overall_report

# Steps:
# 1. fastp
# 2. bowtie2
# 3. Kraken2
# 4. (can be parallel)
#    - Metaphlan4
#    - Virus bowtie2
#    - Fungi bowtie2
# 5. RGI
# 6. VirStrain

def run_qc(args, configs):
    wf = workflows.QualityControl(
        input_1=args.input_1,
        input_2=args.input_2,
        out_dir=configs.qc_dir,
        fileHeader=configs.fileHeader,

        tool_path=envs.FASTP_PATH,
        conda_dir=envs.CONDA_PATH,
        env_name=envs.RTTAP_ENV_NAME,

        threads=args.threads,
        min_len=args.min_len,
        cmd=args.cmd,

        force=args.force
    )
    wf.run()

def run_remove_host(args, configs):
    wf = workflows.HostRemoval(
        input_1=configs.clean_fq_1,
        input_2=configs.clean_fq_2,
        out_dir=configs.nohost_dir,
        fileHeader=configs.fileHeader,

        tool_path=envs.BOWTIE2_PATH,
        db_dir=envs.REF_HOST_DB,
        conda_dir=envs.CONDA_PATH,
        env_name=envs.RTTAP_ENV_NAME,

        threads=args.threads,
        cmd=args.cmd,

        force=args.force
    )
    wf.run()

def run_recruit_reads(args, configs):
    wf = workflows.ReadRecruitment(
        input=configs.nohost_fq,
        out_dir=configs.k2_dir,
        tool_path=envs.KRAKEN2_PATH,
        db_path=envs.NT_MICROBIAL,
        conda_dir=envs.CONDA_PATH,
        env_name=envs.RTTAP_ENV_NAME,
        fileHeader=configs.fileHeader,

        threads=args.threads,
        cmd=args.cmd,

        force=args.force
    )
    wf.run()

def run_bacteria(args, configs):
    wf1 = workflows.ReadExtractor(
        kraken_out=configs.k2_out,
        kraken_report=configs.k2_report,
        nohost_fq_gz=configs.nohost_fq,
        out_dir=configs.bac_dir,
        fileHeader=configs.fileHeader,

        conda_dir=envs.CONDA_PATH,
        env_name=envs.RTTAP_ENV_NAME,
        split_script_path=envs.SPLIT_SCRIPT_PATH,

        threads=args.threads,
        reads_type="bacteria",

        force=args.force
    )
    wf1.run()

    wf2 = workflows.TaxonomyBacteriaWorkflow(
        input=configs.bac_fq,
        k2_report_path=configs.k2_report,
        mpa_report_path=configs.mpa_report,
        bac_report_path=configs.bac_report,

        out_dir=configs.bac_dir,
        fileHeader=configs.fileHeader,

        tool_path=envs.METAPHLAN4_PATH,
        metaphlan4_db_path=envs.METAPHLAN4_DB_PATH,
        conda_dir=envs.CONDA_PATH,
        env_name=envs.RTTAP_ENV_NAME,

        threads=args.threads,
        RPM_threshold=args.rpm_b,
        cmd=args.cmd,

        force=args.force
    )
    wf2.run()

    # if not args.skip_rgi:
    #     wf = workflows.ARGIdentification(
    #         input_fq=configs.bac_fq,
    #         out_dir=configs.ARG_dir,
    #         fileHeader=configs.fileHeader,

    #         tool_path=envs.RGI_PATH,
    #         db_path=envs.CARD_DB_PATH,
    #         conda_dir=envs.CONDA_PATH,
    #         env_name=envs.RTTAP_ENV_NAME,

    #         threads=args.threads,
    #         cmd=args.cmd,

    #         force=args.force
    #     )
    #     wf.run()

def run_viruses(args, configs):
    wf1 = workflows.ReadExtractor(
        kraken_out=configs.k2_out,
        kraken_report=configs.k2_report,
        nohost_fq_gz=configs.nohost_fq,
        out_dir=configs.vir_dir,
        fileHeader=configs.fileHeader,

        conda_dir=envs.CONDA_PATH,
        env_name=envs.RTTAP_ENV_NAME,
        split_script_path=envs.SPLIT_SCRIPT_PATH,

        threads=args.threads,
        reads_type="viruses",

        force=args.force
    )
    wf1.run()

    wf2 = workflows.TaxonomyBowtie2Workflow(
        input=configs.vir_fq,
        k2_report_path=configs.k2_report,
        lca_report_path=configs.vir_lca_out,
        bt_report_path=configs.vir_report,
        out_dir=configs.vir_dir,
        fileHeader=configs.fileHeader,

        bt_path=envs.BOWTIE2_PATH,
        taxonkit_path=envs.TAXONKIT_PATH,
        db_type="viruses",
        acc2taxid_dir=envs.ACC2TAXID_DIR,
        conda_dir=envs.CONDA_PATH,
        env_name=envs.RTTAP_ENV_NAME,

        threads=args.threads,
        RPM_threshold=args.rpm_v,
        cmd=args.cmd,

        force=args.force
    )
    wf2.run()

    # if not args.skip_virstrain:
    #     wf3 = workflows.VirusStrainProfiling(
    #         input_fq_gz=configs.vir_fq,
    #         out_dir=configs.virus_strain_dir,
    #         fileHeader=configs.fileHeader,

    #         virus_report_path=configs.vir_report,
    #         tool_path=envs.VIRSTRAIN_PATH,
    #         db_path=envs.VIRSTRAIN_DB_PATH,
    #         virstrain_list_path=envs.VIRSTRAIN_DB_LIST,
    #         conda_dir=envs.CONDA_PATH,
    #         env_name=envs.RTTAP_ENV_NAME,

    #         threads=args.threads,
    #         cmd=args.cmd,

    #         force=args.force
    #     )
    #     wf3.run()

def run_fungi(args, configs):
    wf1 = workflows.ReadExtractor(
        kraken_out=configs.k2_out,
        kraken_report=configs.k2_report,
        nohost_fq_gz=configs.nohost_fq,
        out_dir=configs.fun_dir,
        fileHeader=configs.fileHeader,

        conda_dir=envs.CONDA_PATH,
        env_name=envs.RTTAP_ENV_NAME,
        split_script_path=envs.SPLIT_SCRIPT_PATH,

        threads=args.threads,
        reads_type="fungi",

        force=args.force
    )
    wf1.run()

    wf2 = workflows.TaxonomyBowtie2Workflow(
        input=configs.fun_fq,
        k2_report_path=configs.k2_report,
        lca_report_path=configs.fun_lca_out,
        bt_report_path=configs.fun_report,
        out_dir=configs.fun_dir,
        fileHeader=configs.fileHeader,

        bt_path=envs.BOWTIE2_PATH,
        taxonkit_path=envs.TAXONKIT_PATH,
        db_type="fungi",
        acc2taxid_dir=envs.ACC2TAXID_DIR,
        conda_dir=envs.CONDA_PATH,
        env_name=envs.RTTAP_ENV_NAME,

        threads=args.threads,
        RPM_threshold=args.rpm_f,
        cmd=args.cmd,

        force=args.force
    )
    wf2.run()

def run_arg(args, configs):
    wf = workflows.ARGIdentification(
        input_fq=configs.bac_fq,
        out_dir=configs.ARG_dir,
        fileHeader=configs.fileHeader,

        tool_path=envs.RGI_PATH,
        db_path=envs.CARD_DB_PATH,
        conda_dir=envs.CONDA_PATH,
        env_name=envs.RGI_ENV_NAME,

        threads=args.threads,
        cmd=args.cmd,

        force=args.force
    )
    wf.run()

def run_virstrain(args, configs):
    wf = workflows.VirusStrainProfiling(
        input_fq_gz=configs.vir_fq,
        out_dir=configs.virus_strain_dir,
        fileHeader=configs.fileHeader,

        virus_report_path=configs.vir_report,
        tool_path=envs.VIRSTRAIN_PATH,
        db_path=envs.VIRSTRAIN_DB_PATH,
        virstrain_list_path=envs.VIRSTRAIN_DB_LIST,
        conda_dir=envs.CONDA_PATH,
        env_name=envs.RTTAP_ENV_NAME,

        threads=args.threads,
        cmd=args.cmd,

        force=args.force
    )
    wf.run()

def run_report(args, configs):
    wf = workflows.GenerateReport(
        out_dir=configs.sample_dir,
        
        k2_report_path=configs.k2_report,
        mpa_report_path=configs.mpa_report,
        vir_lca_out_path=configs.vir_lca_out,
        fun_lca_out_path=configs.fun_lca_out,
        bac_report_path=configs.bac_report,
        fun_report_path=configs.fun_report,
        vir_report_path=configs.vir_report,
        overall_report_path=configs.overall_report,
        bac_rpm=args.rpm_b,
        vir_rpm=args.rpm_v,
        fun_rpm=args.rpm_f,
        bac=args.bac,
        vir=args.vir,
        fun=args.fun,
        overall=args.overall,
        all=args.all,

        force=args.force,
        threads=args.threads,
    )
    wf.run()

def run_end_to_end(args, configs):
    run_qc(args, configs)
    run_remove_host(args, configs)
    run_recruit_reads(args, configs)

    tasks = [run_bacteria, run_viruses, run_fungi]
    processes = [mp.Process(target=task, args=(args, configs)) for task in tasks]
    for p in processes:
        p.start()
    for p in processes:
        p.join()

    run_report(args, configs)

    if not args.skip_rgi:
        run_arg(args, configs)
    if not args.skip_virstrain:
        run_virstrain(args, configs)

class Configs:
    def __init__(self, output_dir, input_2=None):
        self.sample_dir = output_dir
        self.fileHeader = os.path.basename(output_dir)

        # qc_dir
        self.qc_dir = os.path.join(output_dir, '0_quality_control')
        if input_2 is not None:
            self.clean_fq_1 = os.path.join(self.qc_dir, f'{self.fileHeader}.clean.R1.fq.gz')
            self.clean_fq_2 = os.path.join(self.qc_dir, f'{self.fileHeader}.clean.R2.fq.gz')
        else:
            self.clean_fq_1 = os.path.join(self.qc_dir, f'{self.fileHeader}.clean.fq.gz')
            self.clean_fq_2 = None
        self.qc_html = os.path.join(self.qc_dir, f'{self.fileHeader}.html')
        self.qc_json = os.path.join(self.qc_dir, f'{self.fileHeader}.json')

        self.nohost_dir = os.path.join(output_dir, '1_no_host')
        self.nohost_fq = os.path.join(self.nohost_dir, f'{self.fileHeader}.nohost.fq.gz')

        self.k2_dir = os.path.join(output_dir, '2_kraken2_results')
        self.k2_out = os.path.join(self.k2_dir, f'{self.fileHeader}.kraken2.out')
        self.k2_report = os.path.join(self.k2_dir, f'{self.fileHeader}.kraken2.report')

        self.bac_dir = os.path.join(output_dir, '3_1_bacteria_results')
        self.bac_fq = os.path.join(self.bac_dir, f'{self.fileHeader}.bacteria.fq.gz')
        self.mpa_bt2_out = os.path.join(self.bac_dir, f'{self.fileHeader}.metaphlan4.bowtie2.out')
        self.mpa_report = os.path.join(self.bac_dir, f'{self.fileHeader}.metaphlan4.report')
        self.bac_report = os.path.join(self.bac_dir, f'{self.fileHeader}.bacteria.report')

        self.vir_dir = os.path.join(output_dir, '3_2_viruses_results')
        self.vir_fq = os.path.join(self.vir_dir, f'{self.fileHeader}.viruses.fq.gz')
        self.vir_bt2_out = os.path.join(self.vir_dir, f'{self.fileHeader}.viruses.bowtie2.all.out')
        self.vir_lca_out = os.path.join(self.vir_dir, f'{self.fileHeader}.viruses.LCA.out')
        self.vir_report = os.path.join(self.vir_dir, f'{self.fileHeader}.viruses.report')

        self.fun_dir = os.path.join(output_dir, '3_3_fungi_results')
        self.fun_fq = os.path.join(self.fun_dir, f'{self.fileHeader}.fungi.fq.gz')
        self.fun_bt2_out = os.path.join(self.fun_dir, f'{self.fileHeader}.fungi.bowtie2.all.out')
        self.fun_lca_out = os.path.join(self.fun_dir, f'{self.fileHeader}.fungi.LCA.out')
        self.fun_report = os.path.join(self.fun_dir, f'{self.fileHeader}.fungi.report')

        self.overall_report = os.path.join(output_dir, f'{self.fileHeader}.overall.report')

        self.ARG_dir = os.path.join(output_dir, '4_1_ARG_results')
        self.virus_strain_dir = os.path.join(output_dir, '4_2_virus_strain_results')


def main():
    parser = argparse.ArgumentParser(
        prog="RTTAP",
        description="A Read-based Total-infecTome Analysis Pipeline."
    )
    
    # Create subparsers for different modules
    subparsers = parser.add_subparsers(
        title="Modules",
        dest="modules",
        description="Modules that proceed different functions.", 
        help="Please specify one and only one of these options.", 
        required=True
    )
    # Common arguments for all commands
    # io_parser: -i -I -o -t --cmd -F
    io_parser = argparse.ArgumentParser(add_help=False)
    io_parser.add_argument(
        "-i","--input_1",
        required=True,
        help="The first input file (\".fq\", \".fq.gz\")."
    )
    io_parser.add_argument(
        "-I","--input_2",
        required=False, default=None,
        help="The second input file (\".fq\", \".fq.gz\"; used for paired-end files)."
    )
    io_parser.add_argument(
        "-o","--output_dir",
        required=True,
        help="Output directory for results."
    )
    io_parser.add_argument(
        "-t", "--threads",
        required=False, default=1, type=int,
        help="Threads used for analysis. Default: 1"
    )
    io_parser.add_argument(
        "-F", "--force", default=False, action="store_true",
        help="Force to rerun the selected module even if the output files already exist. Default: False"
    )
    io_parser.add_argument(
        "--cmd", type=str, default='',
        help="The string of commands to be added to runnning tools. Default: ''. For example, when runnning kraken2 in remove_host module, '--cmd \"--confidence 0.5\"' will set the minimum confidence to 0.5 for Kraken2 to classify a sequence. This parameter should be used very carefully as it allows great flexibility while introducing a lot of uncertainty for smoothly running the whole pipeline."
    )
    # Common parameters for different modules
    # param_parser_qc: --min_len
    param_parser_qc = argparse.ArgumentParser(add_help=False)
    param_parser_qc.add_argument("--min_len", default=15, type=int,help="Minimum length to keep reads after QC. Default: 15")
    # param_parser_bac: --rpm_b --skip_rgi
    param_parser_bac = argparse.ArgumentParser(add_help=False)
    param_parser_bac.add_argument(
        "--rpm_b", default=10, type=int,
        help="The cut-off threshold for bacterial categories. Default: 10"
    )
    param_parser_bac.add_argument(
        "--skip_rgi", default=False, action="store_true",
        help="Skip the ARG analysis. Default: False"
    )
    # param_parser_vir: --rpm_v --skip_virstrain
    param_parser_vir = argparse.ArgumentParser(add_help=False)
    param_parser_vir.add_argument(
        "--rpm_v", default=1, type=int,
        help="The cut-off threshold for viral categories. Default: 1"
    )
    param_parser_vir.add_argument(
        "--skip_virstrain", default=False, action="store_true",
        help="Skip the viruses strain analysis. Default: False"
    )
    # param_parser_fun: --rpm_f
    param_parser_fun = argparse.ArgumentParser(add_help=False)
    param_parser_fun.add_argument(
        "--rpm_f", default=10, type=int,
        help="The cut-off threshold for fungal categories. Default: 10"
    )
    # param_parser_report: --bac --vir --fun --overall --all
    param_parser_report = argparse.ArgumentParser(add_help=False)
    param_parser_report.add_argument('-b', '--bac', default=False, action="store_true", help='Include bacterial report. Default: False')
    param_parser_report.add_argument('-v', '--vir', default=False, action="store_true", help='Include viral report. Default: False')
    param_parser_report.add_argument('-f', '--fun', default=False, action="store_true", help='Include fungal report. Default: False')
    param_parser_report.add_argument('--overall', default=True, action="store_true", help='Include overall report. Default: True')
    param_parser_report.add_argument('-a', '--all', default=False, action="store_true", help='Include all reports. Default: False')

    # QC module
    qc_subparser = subparsers.add_parser('qc', parents=[io_parser, param_parser_qc], help='Run quality control only')
    # host removal module
    host_subparser = subparsers.add_parser('remove_host', parents=[io_parser], help='Run host removal only')
    # Kraken2 module
    kraken_subparser = subparsers.add_parser('recruit_reads', parents=[io_parser], help='Run Kraken2 classification only')

    # Bacterial analysis module
    bacteria_subparser = subparsers.add_parser('bacteria', parents=[io_parser, param_parser_bac], help='Run bacterial analysis only')
    # Viral analysis module
    viral_subparser = subparsers.add_parser('viruses', parents=[io_parser, param_parser_vir], help='Run viral analysis only')
    # Fungal analysis module
    fungal_subparser = subparsers.add_parser('fungi', parents=[io_parser, param_parser_fun], help='Run fungal analysis only')

    # ARG analysis module
    arg_subparser = subparsers.add_parser('ARG', parents=[io_parser], help='Run ARG analysis only')
    # Virus strain analysis module
    virstrain_subparser = subparsers.add_parser('virstrain', parents=[io_parser], help='Run virus strain analysis only')

    # Overall report module
    report_subparser = subparsers.add_parser(
        'report', parents=[io_parser, param_parser_bac, param_parser_vir, param_parser_fun, param_parser_report], 
        help='Generate report(s) for viruses, bacteria, fungi, overall, or all (specify with -v, -b, -f, --overall, or -a, options respectively). Default: overall.'
    )
    # report_subparser.add_argument('-b', '--bac', default=False, action="store_true", help='Include bacterial report. Default: False')
    # report_subparser.add_argument('-v', '--vir', default=False, action="store_true", help='Include viral report. Default: False')
    # report_subparser.add_argument('-f', '--fun', default=False, action="store_true", help='Include fungal report. Default: False')
    # report_subparser.add_argument('--overall', default=True, action="store_true", help='Include overall report. Default: True')
    # report_subparser.add_argument('-a', '--all', default=False, action="store_true", help='Include all reports. Default: False')
    # End-to-end analysis
    endtoend_subparser = subparsers.add_parser('end_to_end', parents=[io_parser, param_parser_qc, param_parser_bac, param_parser_vir, param_parser_fun, param_parser_report], help='Run complete end-to-end analysis')
    
    args = parser.parse_args()
    
    # relative path to absolute path
    args.input_1 = os.path.abspath(args.input_1)    # -i
    if args.input_2 is not None:
        args.input_2 = os.path.abspath(args.input_2)    # -I
    args.output_dir = os.path.abspath(args.output_dir)  # -o

    configs = Configs(args.output_dir, args.input_2)

    # Run the selected module
    os.makedirs(args.output_dir, exist_ok=True)
    if args.modules == 'qc':
        run_qc(args, configs)
    if args.modules == 'remove_host':
        run_remove_host(args, configs)
    if args.modules == 'recruit_reads':
        run_recruit_reads(args, configs)

    if args.modules == 'bacteria':
        run_bacteria(args, configs)
    if args.modules == 'viruses':
        run_viruses(args, configs)
    if args.modules == 'fungi':
        run_fungi(args, configs)

    if args.modules == 'ARG':
        run_arg(args, configs)
    if args.modules == 'virstrain':
        run_virstrain(args, configs)

    if args.modules == 'report':
        run_report(args, configs)

    if args.modules == 'end_to_end':
        run_end_to_end(args, configs)

if __name__=="__main__":
    main()