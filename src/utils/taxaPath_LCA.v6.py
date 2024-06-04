import pandas as pd
import numpy as np
import sys

fileHeader = sys.argv[1]
readsType = sys.argv[2]

def readLCA(path):
    splittedPath = path.str.split(";", expand=True)
    return ";".join(splittedPath.iloc[0,:][splittedPath.apply(pd.Series.nunique)==1].tolist())
    
if readsType == "viruses":
    resultsPath = "/computenodes/node35/team3/wjiang34/pipeline/downstream/Bowtie2Results_Viruses"
if readsType == "fungi":
    resultsPath = "/computenodes/node35/team3/wjiang34/pipeline/downstream/Bowtie2Results_Fungi"
if readsType == "archaea":
    resultsPath = "/computenodes/node35/team3/wjiang34/pipeline/downstream/Bowtie2Results_Archaea"

if readsType in ["viruses", "fungi", "archaea"]:
    read2taxaPath = pd.read_csv("{}/{}.{}.multi.reads2taxid_path.tmp".format(resultsPath, fileHeader, readsType), header=None, names=["readId", "refId", "taxid", "namePath", "taxaPath"], sep='\t')



if readsType in ["viruses", "fungi", "archaea"]:
    reads_new = read2taxaPath.groupby("readId").agg({"namePath": readLCA, "taxaPath": readLCA}).reset_index()
    reads_new.to_csv("{}/{}.{}.multi_reads_path.tmp".format(resultsPath, fileHeader, readsType), header=None, index=None, sep='\t')