########################## file information #############################
# Author: Bochao Zhang
# Date: 2025-12-5
# Description:
# Determine MHC gene LOH status
#########################################################################

import os,argparse
import pandas as pd
import numpy as np

def parse_args():
    parser = argparse.ArgumentParser(description = "A Python script for determination of MHC region LOH")
    parser.add_argument("-o", help="Output path (REQUIRED)", dest="OUTPUT_PATH")
    parser.add_argument("-b", help="Path to bumps file (REQUIRED)", dest="BUMP_PATH")
    parser.add_argument("-c", help="Path to copy number file (REQUIRED)", dest="CN_PATH")
    parser.add_argument("-e", help="Path to ensembl file (REQUIRED)", dest="ENSEMBL_PATH")
    parser.add_argument("-t", help="Hetereozygous SNP proportion needed to call LOH, default 0.4", dest="LOH_THRESHOLD")
    return parser.parse_args()

def get_mhc_genes(ensembl_path):
    """
    Read MHC gene positions from ensembl
    """
    df = pd.read_csv(ensembl_path)
    gene_pos = {}
    for i in range(len(df)):
        gene_pos[df["GeneName"].iloc[i]] = pd.Interval(df["GeneStart"].iloc[i],df["GeneEnd"].iloc[i],closed="both")
    return gene_pos

def get_loh_bumps(bump_path,cn_path):
    """
    Get bumps with LOH status
    """
    bumps = pd.read_csv(bump_path)
    cn = pd.read_csv(cn_path)
    loh_bumps = []
    for i in range(len(bumps)):
        L = sum(cn["pos"].isin(range(bumps["start"].iloc[i],bumps["end"].iloc[i]+1)))
        cn1 = cn["tumor1"][cn["pos"].isin(range(bumps["start"].iloc[i],bumps["end"].iloc[i]+1))].median()
        cn2 = cn["tumor2"][cn["pos"].isin(range(bumps["start"].iloc[i],bumps["end"].iloc[i]+1))].median()
        if min(cn1,cn2) < 0.5 and max(cn1,cn2) > 0.5 and L >= 3:
                loh_bumps.append(pd.Interval(bumps["start"].iloc[i],bumps["end"].iloc[i],closed="both"))
    return loh_bumps

def get_gene_loh(gene_pos,loh_bumps,cn_path,loh_threshold):
    """
    Determine gene level LOH
    """
    cn = pd.read_csv(cn_path)
    gene_loh = []
    for key in gene_pos:
        het_n = sum(cn["pos"].isin(range(gene_pos[key].left,gene_pos[key].right+1)))
        if het_n > 3:
            gene_loh_het = sum([get_loh_het(gene_pos[key],bump,cn) for bump in loh_bumps])
            if gene_loh_het > het_n * loh_threshold:
                gene_loh.append(True)
            else:
                gene_loh.append(False)
        else:
            gene_loh.append(np.nan)
    loh_status = pd.DataFrame({"gene":gene_pos.keys(), "LOH":gene_loh})
    return loh_status
        
def get_loh_het(gene,bump,cn):
    """
    Return number of heterozygous SNPs in LOH bumps in gene
    """
    if not gene.overlaps(bump):
        return 0
    else:
        return sum(cn["hg38_pos"].isin(range(max(gene.left,bump.left),min(gene.right,bump.right)+1)))

def main():
    args = parse_args()
    if not args.OUTPUT_PATH:
        exit("Specify output directory using -o")
    else:
        output_path = args.OUTPUT_PATH

    if not args.BUMP_PATH:
        exit("Specify path to bump file using -b")
    else:
        bump_path = args.BUMP_PATH

    if not args.CN_PATH:
        exit("Specify path to copynumber file using -c")
    else:
        cn_path = args.CN_PATH

     if not args.ENSEMBL_PATH:
        exit("Specify path to ensemble file using -e")
    else:
        ensembl_path = args.ENSEMBL_PATH

    if not args.LOH_THRESHOLD:
        loh_threshold = 0.4
    else:
        loh_threshold = args.LOH_THRESHOLD

    gene_pos = get_mhc_genes(ensembl_path)
    loh_bumps = get_loh_bumps(bump_path,cn_path)
    gene_loh = get_gene_loh(gene_pos,loh_bumps,cn_path,loh_threshold)
    gene_loh.to_csv(output_path,index=False)

if __name__ == "__main__":
    main()

