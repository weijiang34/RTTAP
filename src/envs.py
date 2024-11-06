import os
### For users modification START
FASTP_PATH = "/computenodes/node35/team3/wjiang34/tools/fastp"
BOWTIE2_PATH = "/home/wjiang34/anaconda3/envs/mpa/bin/bowtie2"
KRAKEN2_PATH = "/computenodes/node35/team3/wjiang34/tools/Kraken2/kraken2"
METAPHLAN4_PATH = "/home/wjiang34/anaconda3/envs/mpa/bin/metaphlan"
METAPHLAN4_DB_PATH = "/home/wjiang34/anaconda3/envs/mpa/lib/python3.10/site-packages/metaphlan/metaphlan_databases/mpa_vJan21_CHOCOPhlAnSGB_202103"
RGI_PATH = "/computenodes/node35/team3/wjiang34/tools/rgi/rgi"
VIRSTRAIN_PATH = "/computenodes/node35/team3/wjiang34/tools/VirStrain/VirStrain.py"
TAXONKIT_PATH = "/home/wjiang34/anaconda3/envs/mpa/bin/taxonkit"
RTTAP_PATH = "/computenodes/node35/team3/wjiang34/RTTAP"
### For users modification END

SPLIT_SCRIPT = os.path.join(RTTAP_PATH, 'src', 'utils', 'extract_kraken_reads_nostdout.py')
NT_MICROBIAL = os.path.join(RTTAP_PATH, 'src', 'databases', 'kraken2ntmicrodb')
REF_HUMAN_DB = os.path.join(RTTAP_PATH, 'src', 'databases', 'humanrRNA_btdb', 'btHumanrRNAdb')
BT_VIRAL_DB_PATH = os.path.join(RTTAP_PATH, 'src', 'databases', 'viral_btdb', 'viral_btdb')
BT_FUNGAL_DB_PATH = os.path.join(RTTAP_PATH, 'src', 'databases', 'EuPathDB46_btdb', 'EuPathDB46_Clean')
CARD_DB_PATH = os.path.join(RTTAP_PATH, 'src', 'databases', 'CARD')
VIRSTRAIN_DB_PATH = os.path.join(RTTAP_PATH, 'src', 'databases', 'virstrain_db')
VIRSTRAIN_DB_LIST = os.path.join(RTTAP_PATH, 'src', 'VirStrain_namelist.txt')
ACC2TAXID_PATH = os.path.join(RTTAP_PATH, 'src', 'dependencies')