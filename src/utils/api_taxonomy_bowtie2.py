import os
import shlex

import logging
logging.basicConfig(
    level=logging.DEBUG, 
    format="%(asctime)s [%(levelname)s]: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

import pandas as pd
import numpy as np
import dask.dataframe as dd

import envs
from .system_operations import safe_run_cmd, safe_remove

def run_mapping_and_extract(
    input, 
    out_dir, 
    fileHeader, 

    bt_path, 
    db_type, 
    db_path, 
    conda_dir, 
    env_name,
    
    threads,
    cmd=''
):
    """Run Bowtie2 and extract uniq/multi mapping results."""
    script_lines = [
        "#!/bin/bash",
        f"source {os.path.join(conda_dir, 'bin', 'activate')} {env_name}",
        f"{shlex.quote(bt_path)} -x {shlex.quote(db_path)} -U {shlex.quote(input)} --xeq -S {shlex.quote(os.path.join(out_dir, f'{fileHeader}.{db_type}.bowtie2.all.out'))} -a -p {threads} --no-unal --no-head {cmd}",
        f"awk '{{if ($0 ~ /AS:/ && $0 !~ /XS:/) print $1 \"\\t\" $3 > \"{os.path.join(out_dir, f'{fileHeader}.{db_type}.uniq.tmp')}\"; else if ($0 ~ /AS:/ && $0 ~ /XS:/) print $1 \"\\t\" $3 > \"{os.path.join(out_dir, f'{fileHeader}.{db_type}.multi.tmp')}\"}}' {os.path.join(out_dir, f'{fileHeader}.{db_type}.bowtie2.all.out')}"
    ]
    script_lines = [line + '\n' for line in script_lines]
    script_path = os.path.join(out_dir, f"run_Bowtie2_{db_type}_tmp.sh")
    with open(script_path, 'w') as f:
        f.writelines(script_lines)
    os.chmod(script_path, 0o755)
    safe_run_cmd(shlex.quote(script_path), log_file=os.path.join(out_dir, f"{fileHeader}.bowtie2.{db_type}.log"))
    safe_remove(script_path)

def append_taxid(
    out_dir, 
    db_type, 
    acc2taxid_path, 
    fileHeader
):
    """Append taxid to uniq/multi mapping results."""
    acc2taxid = pd.read_table(acc2taxid_path, header=None, names=["acc.ver", "taxid"]).astype({"taxid": "str"})
    for tag in ["uniq", "multi"]:
        tmp_path = os.path.join(out_dir, f"{fileHeader}.{db_type}.{tag}.tmp")
        assert os.path.exists(tmp_path), f"Missing mapping file: {tmp_path}"
        read2acc = pd.read_table(tmp_path, header=None, names=["readId", "acc.ver"])
        merged = pd.merge(read2acc, acc2taxid, how="left", on="acc.ver")
        out_path = os.path.join(out_dir, f"{fileHeader}.{db_type}.{tag}.reads2acc2taxid.tmp")
        merged.to_csv(out_path, header=None, index=False, sep='\t')

def convert_taxid_to_path(
    out_dir, 
    db_type, 
    taxonkit_path, 
    fileHeader
):
    """Use TaxonKit to convert taxid to phylogenetic path."""
    for tag in ["uniq", "multi"]:
        in_path = os.path.join(out_dir, f"{fileHeader}.{db_type}.{tag}.reads2acc2taxid.tmp")
        out_path = os.path.join(out_dir, f"{fileHeader}.{db_type}.{tag}.lineage.tmp")
        assert os.path.exists(in_path), f"Missing input for taxonkit: {in_path}"
        cmd = f"cat {shlex.quote(in_path)} | {shlex.quote(taxonkit_path)} reformat -I3 -t > {shlex.quote(out_path)}"
        safe_run_cmd(cmd, log_file=os.path.join(out_dir, f"{fileHeader}.{db_type}.{tag}.taxonkit.log"))

def run_lca(
    out_dir, 
    db_type, 
    fileHeader
):
    """Calculate LCA for multi-mapping results."""
    def lineage_lca(lineage):
        splitted = lineage.str.split(";")
        splitted = splitted.apply(lambda x: [i for i in x if pd.notnull(i)])
        arrs = np.array(splitted.tolist())
        mask = (arrs == arrs[0]).all(axis=0)
        return ";".join(arrs[0][mask])
    in_path = os.path.join(out_dir, f"{fileHeader}.{db_type}.multi.lineage.tmp")
    ddf = dd.read_csv(in_path, sep='\t', names=["readId", "refId", "taxid", "namePath", "taxidPath"])

    def lca_group(df):
        return pd.Series({
            "readId": df["readId"].iloc[0],
            # "refId": df["refId"].iloc[0],
            "namePath": lineage_lca(df["namePath"]),
            "taxidPath": lineage_lca(df["taxidPath"])
        })
    meta = {
        "readId": "object",
        # "refId": "object",
        "namePath": "object",
        "taxidPath": "object"
    }
    result = ddf.groupby("readId").apply(lca_group, meta=meta, include_groups=True).compute()

    out_path = os.path.join(out_dir, f"{fileHeader}.{db_type}.multi.LCA.tmp")

    result.to_csv(out_path, header=None, index=None, sep='\t')

def merge_and_cleanup(
    out_dir, db_type, fileHeader
):
    """Merge uniq/multi results and safely clean up temporary files."""
    uniq_path = os.path.join(out_dir, f"{fileHeader}.{db_type}.uniq.lineage.tmp")
    multi_path = os.path.join(out_dir, f"{fileHeader}.{db_type}.multi.LCA.tmp")

    uniq = pd.read_csv(uniq_path, sep='\t', header=None)
    uniq = uniq.iloc[:,[0,3,4]] # ,1

    multi = pd.read_csv(multi_path, sep='\t', header=None)

    merged = pd.concat([uniq, multi]).drop_duplicates()

    out_path = os.path.join(out_dir, f"{fileHeader}.{db_type}.LCA.out")
    merged.to_csv(out_path, header=None, index=None, sep='\t')

    temp_files = [
        os.path.join(out_dir, f"{fileHeader}.{db_type}.uniq.tmp"),
        os.path.join(out_dir, f"{fileHeader}.{db_type}.multi.tmp"),
        os.path.join(out_dir, f"{fileHeader}.{db_type}.uniq.reads2acc2taxid.tmp"),
        os.path.join(out_dir, f"{fileHeader}.{db_type}.multi.reads2acc2taxid.tmp"),
        os.path.join(out_dir, f"{fileHeader}.{db_type}.uniq.lineage.tmp"),
        os.path.join(out_dir, f"{fileHeader}.{db_type}.multi.lineage.tmp"),
        os.path.join(out_dir, f"{fileHeader}.{db_type}.multi.LCA.tmp"),
    ]
    for temp_file in temp_files:
        safe_remove(temp_file)

def generate_bt2_report(
    # out_dir: str, 
    # fileHeader: str, 
    # db_type: str, 
    lca_out_path: str, 
    k2_report_path: str, 
    bt_report_path: str,
    RPM_threshold: float = 1.0
) -> None:
    """
    Generate final taxonomic report from LCA and Kraken2 results.
    
    Args:
        out_dir: Output directory path
        fileHeader: File header prefix for input/output files
        db_type: Database type identifier
        k2_report_path: Path to Kraken2 report file
        RPM_threshold: Minimum RPM threshold for filtering (default: 1.0)
    
    Returns:
        None: Writes the report to a file
    """
    try:
        # Read LCA results
        # lca_out_path = os.path.join(out_dir, f"{fileHeader}.{db_type}.LCA.out")
        lca_out = pd.read_table(
            lca_out_path,
            header=None,
            names=["readID", "namePath", "taxidPath"]   # , "refId"
        )
        # Read Kraken2 report to get total reads count
        k2_report = pd.read_table(
            k2_report_path, 
            sep='\t',
            header=None,
            names=["percent", "reads", "reads_direct", "rank", "taxid", "taxa"]
        )
        # Calculate total reads from root and unclassified entries
        total_reads_mask = k2_report["taxa"].isin(["root", "unclassified"])
        total_reads_count = k2_report.loc[total_reads_mask, "reads"].sum()
        # Handle case with no LCA results
        if lca_out.empty:
            report = pd.DataFrame({
                "rank": [],
                "taxa": [], 
                "taxid": [],
                "reads": [],
                "RPM": []
            })
        else:
            # Process each taxonomic rank
            ranks = ["SuperKingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
            report_data = []
            
            # Split name and taxid paths
            name_paths = lca_out["namePath"].str.split(";", expand=True)
            name_paths = name_paths.fillna('')
            name_paths = name_paths[~(name_paths == '').all(axis=1)]
            taxid_paths = lca_out["taxidPath"].str.split(";", expand=True)
            taxid_paths = taxid_paths.fillna('')
            taxid_paths = taxid_paths[~(taxid_paths == '').all(axis=1)]

            max_columns = min(name_paths.shape[1], len(ranks))
            
            # Process each available rank level
            for idx in range(max_columns):
                rank_name = ranks[idx]
                
                # Extract taxa names and taxids for this rank
                taxa_names = name_paths[idx][name_paths[idx] != '']
                taxids = taxid_paths[idx][taxid_paths[idx] != '']
                
                # Create temporary DataFrame for this rank
                rank_df = pd.DataFrame({
                    "rank": rank_name,
                    "taxa": taxa_names,
                    "taxid": taxids,
                    "reads": 1  # Each row represents one read
                })
                
                report_data.append(rank_df)
            
            # Combine all rank data
            if report_data:
                combined_data = pd.concat(report_data, ignore_index=True)
                
                # Aggregate reads by rank, taxa, and taxid
                aggregated = (combined_data
                            .groupby(["rank", "taxa", "taxid"], as_index=False)
                            .agg(reads=("reads", "sum")))
                
                # Calculate RPM and filter by threshold
                aggregated["RPM"] = (aggregated["reads"] * 1_000_000) / total_reads_count
                report = aggregated[aggregated["RPM"] >= RPM_threshold].copy()
                
                # Sort by rank order and reads count
                rank_order = {rank: i for i, rank in enumerate(ranks)}
                report["rank_order"] = report["rank"].map(rank_order)
                report = report.sort_values(["rank_order", "reads"], ascending=[True, False])
                report = report.drop("rank_order", axis=1)
            else:
                report = pd.DataFrame({
                    "rank": [],
                    "taxa": [], 
                    "taxid": [],
                    "reads": [],
                    "RPM": []
                })
        # Write final report
        # output_file = os.path.join(out_dir, f"{fileHeader}.{db_type}.report")

        report.to_csv(bt_report_path, sep='\t', index=False)

    except FileNotFoundError as e:
        logging.error(f"Input file not found - {e}")
        raise
    except pd.errors.EmptyDataError:
        logging.error("Error: Input file is empty")
        raise
    except Exception as e:
        logging.error(f"Unexpected error: {e}")
        raise