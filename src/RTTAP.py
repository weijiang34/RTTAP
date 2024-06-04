import argparse
import os
from utils import functions
import envs

# Requirements:
# 1. fastp
# 2. bowtie2
# 3. Kraken2
# 4. Metaphlan4
# 5. RGI
# 6. VirStrain
# 7. taxonkit


def main():
    parser = argparse.ArgumentParser(
        prog="RTTAP",
        description="A Read-based Total-infecTome Analysis Pipeline."
    )
    parser.add_argument(
        "-i","--input_1",
        required=True, default=None,
        help="The first input file (\".fq\", \".fq.gz\")."
    )
    parser.add_argument(
        "-I","--input_2",
        required=False, default=None,
        help="The second input file (\".fq\", \".fq.gz\"; used for paired-end files)."
    )
    parser.add_argument(
        "-o","--output_dir",required=False, type=str,
        help="By default, RTTAP will output to the current folder \"./\", and there will be a folder \"./out/\" contaning all the output/intermediate files."
    )
    parser.add_argument(
        "--rpm_b",required=False, default=10,
        help="The cut-off threshold for bacterial categories. Default: 10"
    )
    parser.add_argument(
        "--rpm_v",required=False, default=1,
        help="The cut-off threshold for viral categories. Default: 1"
    )
    parser.add_argument(
        "--rpm_f",required=False, default=10,
        help="The cut-off threshold for fungal categories. Default: 10"
    )
    parser.add_argument(
        "-t", "--threads",required=False, default=8,
        help="Threads used for analysis. Default: 8"
    )

    args = parser.parse_args()

    # print(args.input_1, args.input_2)
    
    functions.run_pipeline(
        args.input_1, out_dir=args.output_dir, input_2=args.input_2, threads=1, 
        RPM_b=args.rpm_b, RPM_v=args.rpm_v, RPM_f=args.rpm_f,
        env_variables={
            "fastp_path": envs.FASTP_PATH,
            "bowtie2_path": envs.BOWTIE2_PATH,
            "kraken2_path": envs.KRAKEN2_PATH,
            "mpa_path": envs.METAPHLAN4_PATH,
            "rgi_path": envs.RGI_PATH,
            "virstrain_path": envs.VIRSTRAIN_PATH,
            "taxonkit_path": envs.TAXONKIT_PATH,
            "split_script": envs.SPLIT_SCRIPT,
            "CARD_db": envs.CARD_DB_PATH,
            "metaphlan4_db_path": envs.METAPHLAN4_DB_PATH,
            "bt_viral_db_path": envs.BT_VIRAL_DB_PATH,
            "bt_fungal_db_path": envs.BT_FUNGAL_DB_PATH,
            "virstrain_db_path": envs.VIRSTRAIN_DB_PATH,
            "virstrain_db_list": envs.VIRSTRAIN_DB_LIST,
            "nt_microbial": envs.NT_MICROBIAL,
            "ref_human_db": envs.REF_HUMAN_DB,
            "acc2taxid_path": envs.ACC2TAXID_PATH
        }
    )

    # os.system("chmod +x ./runPipeline.v9.sh")
    # os.system("./runPipeline.v9.sh")
    

if __name__=="__main__":
    main()