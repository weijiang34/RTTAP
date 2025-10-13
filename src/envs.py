import os
### For users modification START
CONDA_PATH = "/computenodes/node351/wjiang34/miniconda3"
RTTAP_ENV_NAME = "RTTAP"
RGI_ENV_NAME = "rgi"
RTTAP_PATH = "/computenodes/node35/team3/wjiang34/RTTAP"
RTTAP_DB_PATH = "/computenodes/node35/team3/wjiang34/dataRef/RTTAP_DBs"

# If you installed RTTAP following exactly the steps in documentation, the fllowing will not need to be modified.

FASTP_PATH = "/computenodes/node351/wjiang34/miniconda3/envs/RTTAP/bin/fastp"
TAXONKIT_PATH = "/computenodes/node351/wjiang34/miniconda3/envs/RTTAP/bin/taxonkit"

KRAKEN2_PATH = "/computenodes/node35/team3/wjiang34/tools/Kraken2/kraken2"

METAPHLAN4_PATH = "/computenodes/node351/wjiang34/miniconda3/envs/RTTAP/bin/metaphlan"
METAPHLAN4_DB_PATH = "/computenodes/node351/wjiang34/miniconda3/envs/mpa/lib/python3.10/site-packages/metaphlan/metaphlan_databases/mpa_vJan21_CHOCOPhlAnSGB_202103"
BOWTIE2_PATH = "/computenodes/node351/wjiang34/miniconda3/envs/RTTAP/bin/bowtie2"

RGI_PATH = "/computenodes/node351/wjiang34/miniconda3/envs/rgi/bin/rgi"
VIRSTRAIN_PATH = "/computenodes/node35/team3/wjiang34/tools/VirStrain/VirStrain.py"
### For users modification END

SPLIT_SCRIPT_PATH = os.path.join(RTTAP_PATH, 'src', 'utils', 'extract_kraken_reads_nostdout.py')
NT_MICROBIAL = os.path.join(RTTAP_DB_PATH, 'kraken2ntmicrodb')
REF_HOST_DB = os.path.join(RTTAP_DB_PATH, 'humanrRNA_btdb', 'btHumanrRNAdb')
BT_VIRAL_DB_PATH = os.path.join(RTTAP_DB_PATH, 'viral_btdb', 'viral_btdb')
BT_FUNGAL_DB_PATH = os.path.join(RTTAP_DB_PATH, 'EuPathDB46_btdb', 'EuPathDB46_Clean')
CARD_DB_PATH = os.path.join(RTTAP_DB_PATH, 'CARD')
VIRSTRAIN_DB_PATH = os.path.join(RTTAP_DB_PATH, 'virstrain_db')
VIRSTRAIN_DB_LIST = os.path.join(RTTAP_PATH, 'src', 'VirStrain_namelist.txt')
ACC2TAXID_DIR = os.path.join(RTTAP_DB_PATH, 'seq2taxid')