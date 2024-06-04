import pandas as pd
import numpy as np
import sys

fileHeader = sys.argv[1]
readsType = sys.argv[2]


if readsType in ["viruses", "fungi", "archaea"]:
    if readsType == "viruses":
        resultsPath = "/computenodes/node35/team3/wjiang34/pipeline/downstream/Bowtie2Results_Viruses"
        acc2taxid = pd.read_table("{}.acc2taxid.txt".format(readsType), header=None, names=["acc.ver", "taxid"]).astype({"taxid":"str"})
    if readsType == "fungi":
        resultsPath = "/computenodes/node35/team3/wjiang34/pipeline/downstream/Bowtie2Results_Fungi"
        acc2taxid = pd.read_table("EuPathDB46.seq2taxid.txt", header=None, names=["acc.ver", "taxid"]).astype({"taxid":"str"})
    if readsType == "archaea":
        resultsPath = "/computenodes/node35/team3/wjiang34/pipeline/downstream/Bowtie2Results_Archaea"
        acc2taxid = pd.read_table("{}.acc2taxid.txt".format(readsType), header=None, names=["acc.ver", "taxid"]).astype({"taxid":"str"})
    
    uniqread2acc = pd.read_table("{}/{}.{}.uniq.tmp".format(resultsPath, fileHeader, readsType), header=None, names=["readId", "acc.ver"])
    multiread2acc = pd.read_table("{}/{}.{}.multi.tmp".format(resultsPath, fileHeader, readsType), header=None, names=["readId", "acc.ver"])
    # add taxid after acc.ver column
    uniqacc2taxid = pd.merge(uniqread2acc, acc2taxid, how="left", on="acc.ver")
    multiacc2taxid = pd.merge(multiread2acc, acc2taxid, how="left", on="acc.ver")

# save temporary file
if readsType in ["viruses", "fungi", "archaea"]:

    uniqacc2taxid.to_csv("{}/{}.{}.uniq.reads2acc2taxid.tmp".format(resultsPath, fileHeader, readsType), header=None, index=False, sep='\t')
    multiacc2taxid.to_csv("{}/{}.{}.multi.reads2acc2taxid.tmp".format(resultsPath, fileHeader, readsType), header=None, index=False, sep='\t')
