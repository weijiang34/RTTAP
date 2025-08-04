# RTTAP
The Read-based Total-infectome Taxonomic Analysis Pipeline.  
RTTAP (Read-based Total-infectome Taxonomic Analysis Pipeline) is a fast, accurate, and sensitive pipeline focusing analyses of the totoal-infectome of clinical metatranscriptomic data. It includes multiple useful functions to process and analyze raw sequencing reads: quality control, taxonomy profiling, ARG profiling, and virus strain-level profiling thus providing users comprehensive insights about the microbial composition in clinical samples. It is user friendly and easy-to-use, involving minimum human intervention: all its steps could be finished in a single run. 
# Installation
1. Clone this respirotry to your local computer:
    ```
    git clone https://github.com/weijiang34/RTTAP.git
    ```
2. Install depoendencies:  
- create conda environment:  
    ```
    cd ./RTTAP/
    conda create -f environment.yml
    ```
- Download databases:  
    Download requeired databases from:
    ```
    # Baidu Netdisk
    https://pan.baidu.com/s/1HYhzCu9KFAacUk5FqwhK6A?pwd=k4tg
    # or Google Drive:
    https://drive.google.com/uc?id=1SMx0wD8_z0fn44brOmlDkInM1V6JQxbw (databases)
    https://drive.google.com/uc?id=1r-PQ0x-uLNhgHKxA1TaXXmL_2xM1DAqU (dependencies)
    ```
    Unzip them:
    ```
    tar -zxvf RTTAP_DBs.tar.gz
    ```
    and put them under folder: 
    ```
    mv where/you/downloaded path/to/RTTAP/src/databases/
    ```
    The database should looks like this:
    ```
    --./RTTAP/src/databases/
        |--./CARD/
        |--./EuPathDB46_btdb/
        |--./humanrRNA_btdb/
        |--./kraken2ntmicrodb/
        |--./viral_btdb/
        |--./virstrain_db/
    ```
- Taxonkit, Fastp, MetaPhlAn4, Bowtie2  
    ```
    conda install -c bioconda fastp 
    conda install -c bioconda taxonkit 
    conda install -c bioconda bowtie2 
    conda install -c bioconda metaphlan
    ```
    Detailed instructions for installing them are listed below:  
    - Taxonkit: https://github.com/shenwei356/taxonkit  
    - Fastp: https://github.com/OpenGene/fastp  
    - MetaPhlAn4: https://github.com/biobakery/MetaPhlAn  
    - Bowtie2: https://github.com/BenLangmead/bowtie2  
- Kraken2  
    Download Kraken2 using the following command:
    ```
    git clone https://github.com/DerrickWood/kraken2.git
    ```
    And please refer to the instructions in *Installation* part: https://github.com/DerrickWood/kraken2/wiki/Manual#installation to install Kraken2
- RGI  
    RGI:
    ```
    conda install -c conda-forge -c bioconda -c defaults rgi
    ```
    GitHub: https://github.com/arpcard/rgi
- VirStrain
    ```
    conda install -c bioconda virstrain
    ```
    GiHub: https://github.com/liaoherui/VirStrain
3. Config your environments
    Open ```./RTTAP/src/envs.py``` with text editor and set all the required environmental variables according to you computer's environment. Plese only modify the texts between ```### For users modification START``` and ```### For users modification END``` accordingly.
# Test sample
Once successfully installed all the dependencies, you can run the following command to have a quick test:
- Download test sample from https://portland-my.sharepoint.com/:f:/g/personal/wjiang34-c_my_cityu_edu_hk/El6uED6jVfxNtEG21Jc27XoBNRiDw6Z08uNxxPRHJtXhBA?e=qRWWGE  
Move this test sample to your linux system and run th efollowing command to unzip.
```
tar -zxvf test_sample.tar.gz
cd test_sample/

# run RTTAP
python path/to/RTTAP/src/RTTAP.py -i test.fq.gz -o test-output/ -t 8
```
The test sample is a toy dataset which includes simulated reads from several commonly-seen pathegens' genomes. It is a non-paired fastq file, so you will only need to specify one of the input by: ```-i```.

# User guides
- Overview  
    At any time of using RTTAP, you can use the following command to check each of its parameters and corresponding meaning.
    ```
    python path/to/RTTAP/src/RTTAP.py --help

    # usage: RTTAP [-h] -i INPUT_1 [-I INPUT_2] [-o OUTPUT_DIR] [--rpm_b RPM_B] [--rpm_v RPM_V] [--rpm_f RPM_F] [-t THREADS]
    # 
    # A Read-based Total-infecTome Analysis Pipeline.
    # 
    # options:
    #   -h, --help            show this help message and exit
    #   -i INPUT_1, --input_1 INPUT_1
    #                         The first input file (".fq", ".fq.gz").
    #   -I INPUT_2, --input_2 INPUT_2
    #                         The second input file (".fq", ".fq.gz"; used for paired-end files).
    #   -o OUTPUT_DIR, --output_dir OUTPUT_DIR
    #                         By default, RTTAP will output to the current folder "./", and there will be a folder "./out/" contaning all the output/intermediate files.
    #   --rpm_b RPM_B         The cut-off threshold for bacterial categories. Default: 10
    #   --rpm_v RPM_V         The cut-off threshold for viral categories. Default: 1
    #   --rpm_f RPM_F         The cut-off threshold for fungal categories. Default: 10
    #   -t THREADS, --threads THREADS
    #                         Threads used for analysis. Default: 8
    ```
- Input  
    RTTAP supports both paried-end fastq files and single-end fastq files. 
    If you have paried-end data, you can use the ```-i``` and ```-I``` to specify paired inputs:
    ```
    python path/to/RTTAP/src/RTTAP.py -i input_1.fq.gz -I input_2.fq.gz -o output_path/
    ```
    If your data is single-end, you can use ```-i``` to specify your input and ignore the ```-I```:
    ```
    python path/to/RTTAP/src/RTTAP.py -i input_1.fq.gz -o output_path/
    ```
- Output  
    RTTAP supports user-customized output directories by using ```-o```:
    ```
    python path/to/RTTAP/src/RTTAP.py -i input_1.fq.gz -I input_2.fq.gz -o path/to/your/output_dir
    ```
- Filtering thresholds  
    RTTAP uses several thresholds to filter out extremely low abundance taxonmies. Users can freely adjust these thresholds according to their wishes. 
    ```
    python path/to/RTTAP/src/RTTAP.py -i input_1.fq.gz -I input_2.fq.gz -o output_path/ --rpm_b 10 --rpm_v 1 --rpm_f 10

    # --rpm_b RPM_B         The cut-off threshold for bacterial categories. Default: 10
    # --rpm_v RPM_V         The cut-off threshold for viral categories. Default: 1
    # --rpm_f RPM_F         The cut-off threshold for fungal categories. Default: 10
    ```
- Threads  
    RTTAP supports users to adjust working threads by using ```-t```, depending on their available computing resources.
    ```
    python path/to/RTTAP/src/RTTAP.py -i input_1.fq.gz -I input_2.fq.gz -o output_path/ -t num_of_threads
    ```
    The default valeus is set to ```-t 8```.

# Output interpretation
RTTAP generates several output files which contains useful information for you to explore and analyze. 
- Taxonomic profiles:  
    Viral: ```./test_output/Bacteria_results/test_output.bacteria.report```  
    Bactrial: ```./test_output/Viruses_results/test_output.viruses.report```  
    Fungal: ```./test_output/Fungi_results/test_output.fungi.report```  
    They all share the same ```.tsv``` format with each column names shown in the first row:
    ```
    rank	taxa	taxaID	reads	RPM
    SuperKingdom	Viruses	10239	339	2559.70763457342
    Phylum	Pisuviricota	2732408	265	2000.951395757983
    Phylum	Negarnaviricota	2497569	47	354.88571924764227
    Phylum	Uroviricota	2731618	13	98.15987979190105
    Phylum	Preplasmiviricota	2732008	10	75.50759983992388
    Phylum	Cossaviricota	2732415	4	30.203039935969556
    Class	Pisoniviricetes	2732506	265	2000.951395757983
    Class	Monjiviricetes	2497574	37	279.3781194077184
    ```
- Antibiotic Resistance Genes (ARG) profiles:
    The ARG profiels are under the folder ```./test_output/ARG/```.  
    You can find detected ARGs in ```./test_output/ARG/test_output.rgi.out.gene_mapping_data.txt```
    This files is a ```.tsv``` file with column names presented in the first row:
    ```
    ARO Term	ARO Accession	Reference Model Type	Reference DB	Alleles with Mapped Reads	Reference Allele(s) Identity to CARD Reference Protein (%)	Resistomes & Variants: Observed in Genome(s)	Resistomes & Variants: Observed in Plasmid(s)	Resistomes & Variants: Observed Pathogen(s)	Completely Mapped Reads	Mapped Reads with Flanking Sequence	All Mapped Reads	Average Percent Coverage	Average Length Coverage (bp)	Average MAPQ (Completely Mapped Reads)	Number of Mapped Baits	Number of Mapped Baits with Reads	Average Number of reads per Bait	Number of reads per Bait Coefficient of Variation (%)	Number of reads mapping to baits and mapping to complete gene	Number of reads mapping to baits and mapping to complete gene (%)	Mate Pair Linkage (# reads)	Reference Length	AMR Gene Family	Drug Class	Resistance Mechanism
    macB	3000535	protein homolog model	CARD	1	100.0	no data	no data	Neisseria gonorrhoeae	4.00	0.00	4.00	29.04	562.00	194.00	0	0	0	0	N/A	N/A		1935	ATP-binding cassette (ABC) antibiotic efflux pump	macrolide antibiotic	antibiotic efflux
    TEM-79	3000946	protein homolog model	CARD	1	100.0	no data	no data	Escherichia coli	1.00	0.00	1.00	17.31	149.00	189.00	0	0	0	0	N/A	N/A		861	TEM beta-lactamase	monobactam; cephalosporin; penam; penem	antibiotic inactivation
    mtrC	3000810	protein homolog model	CARD	1	100.0	no data	no data	Neisseria meningitidis	2.00	0.00	2.00	18.32	227.00	168.50	0	0	0	0	N/A	N/A		1239	resistance-nodulation-cell division (RND) antibiotic efflux pump	macrolide antibiotic; penam	antibiotic efflux
    mefE	3000614	protein homolog model	CARD	1	100.0	no data	no data	Streptococcus pneumoniae	3.00	0.00	3.00	29.87	362.00	117.33	0	0	0	0	N/A	N/A		1212	major facilitator superfamily (MFS) antibiotic efflux pump	macrolide antibiotic	antibiotic efflux
    ```
    For more detailed column explanation, please refer to *RGI bwt read mapping results at allele level* in: https://github.com/arpcard/rgi/blob/master/docs/rgi_bwt.rst.
- Viral strain-level profiles  
    Viral strain-level profiels are under ```./test_output/Viral_strains/```. A summary is at ```./test_output/Viral_strains/test_output.virstrain.report```
    ```
    Strain_ID	Cls_info	SubCls_info	Map_Score	Valid_Map_Rate	Total_Map_Rate	Strain_Depth	Strain_info	Unique_SNP
    >MK773085	Cluster2660_1	NA	0.9975520121037461	1014/1033	1014/1033	4602.055555555556			{'153-G': 186, '162-A': 429, ...}
    ```
    For detailed interpretation, please refer to https://github.com/liaoherui/VirStrain.

# Contact  
- Welcome to email to wjiang34-c@my.cityu.edu.hk or leave your thoughts under *Issues*, if you have any questions to RTTAP, or have any troubles during usage, or you have any suggestions to us.
