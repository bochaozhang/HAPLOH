########################## file information #############################
# Author: Bochao Zhang
# Date: 2024-3-11
# Description:
# Estimated allele copy number at each heterozygous site
#########################################################################

import argparse
import os
import pysam
import pandas as pd
import numpy as np

def parse_args():
    parser = argparse.ArgumentParser(description = "A Python script for determination of MHC region LOH")
    parser.add_argument("-o", help="Output path (REQUIRED)", dest="OUTPUT_PATH")
    parser.add_argument("-n1", help="Filtered hap1 normal bam (REQUIRED)", dest="NORMAL_1")
    parser.add_argument("-n2", help="Filtered hap2 normal bam (REQUIRED)", dest="NORMAL_2")
    parser.add_argument("-t1", help="Filtered hap1 tumor bam (REQUIRED)", dest="TUMOR_1")
    parser.add_argument("-t2", help="Filtered hap2 tumor bam (REQUIRED)", dest="TUMOR_2")
    parser.add_argument("-p", help="Path to polymorphic site tsv file (REQUIRED)", dest="POLY_PATH")
    parser.add_argument("-c", help="Path to combined copy number file (REQUIRED)", dest="CN_PATH")
    parser.add_argument("-u", help="Tumor purity (REQUIRED)", dest="PURITY")
    return parser.parse_args()

def get_heterozygous_site(poly_path):
    """
    Read polymorphic site
    """
    df = pd.read_csv(poly_path,sep="\t")
    data = {"hg38_pos":[],"hap1_contig":[],"hap1_pos":[],"hap1_nt":[],"hap2_contig":[],"hap2_pos":[],"hap2_nt":[]}
    for i in range(len(df)):
        data["hg38_pos"].append(df["hg38_pos"][i])
        data["hap1_contig"].append(df["hp1_pos"][i].split("-")[0])
        data["hap1_pos"].append(int(df["hp1_pos"][i].split("-")[1]))
        data["hap1_nt"].append(df["hp1.nt"][i])
        data["hap2_contig"].append(df["hp2_pos"][i].split("-")[0])
        data["hap2_pos"].append(int(df["hp2_pos"][i].split("-")[1]))
        data["hap2_nt"].append(df["hp2.nt"][i])
    het = pd.DataFrame(data)
    het = het[(het["hap1_contig"]!=None) & (het["hap2_contig"]!=None)]
    het.reset_index(drop=True,inplace=True)
    return het

def get_coverage(bam):
    """
    Read bam file coverage
    """
    samfile = pysam.AlignmentFile(bam,"rb")
    contigs = samfile.references
    lengths = samfile.lengths
    coverage = {}
    for i in range(len(contigs)):
        c = samfile.count_coverage(contigs[i],0,lengths[i])
        coverage[contigs[i]] = [c[0][i]+c[1][i]+c[2][i]+c[3][i] for i in range(len(c[0]))]
    samfile.close()
    return coverage

def calculate_baf(tumor_coverage1,tumor_coverage2,het):
    """
    Calculate tumor BAF
    """
    BAF = []
    for i in range(len(het)):
        if het["hap1_pos"][i] <= len(tumor_coverage1[str(het["hap1_contig"][i])]) and het["hap2_pos"][i] <= len(tumor_coverage2[str(het["hap2_contig"][i])]):
            coverage_1 = tumor_coverage1[str(het["hap1_contig"][i])][int(het["hap1_pos"][i])-1]
            coverage_2 = tumor_coverage2[str(het["hap2_contig"][i])][int(het["hap2_pos"][i])-1]
            if coverage_1 + coverage_2 > 0:
                f = coverage_1 / (coverage_1 + coverage_2)
                BAF.append(f)
            else:
                BAF.append(None)
        else:
            BAF.append(None)
    return BAF

def get_combined_copynumber(het,cn_path):
    """
    Read combined copynumber
    """
    df = pd.read_csv(cn_path,sep="\t")
    combined_cn = []
    for i in range(len(het)):
        try:
            combined_cn.append(float(df["ploidy"][df["pos"]==het["hg38_pos"][i]]))
        except TypeError:
            combined_cn.append(None)
    return combined_cn

def calculate_tumor_haplotype_copynumber(BAF,purity,combined_cn):
    """
    Calculate tumor coverage for each pseudo-haplotype
    """
    hap1 = []
    hap2 = []
    for i in range(len(combined_cn)):
        if combined_cn[i] != None:
            if combined_cn[i] < -2:
                allele1.append(None)
                allele2.append(None)
            else:
                if BAF[i] != None:
                    if BAF[i] != 0 and BAF[i] != 1:
                        a1 = ( purity - 1 + BAF[i] * (2*(1-purity) + purity*combined_cn[i]) ) / purity
                        a2 = ( purity - 1 + (1-BAF[i]) * (2*(1-purity) + purity*combined_cn[i]) ) / purity
                    elif BAF[i] == 0:
                        a1,a2 = 0,combined_cn[i]
                    elif BAF[i] == 1:
                        a1,a2 = combined_cn[i],0
                    hap1.append(a1)
                    hap2.append(a2)
                else:
                    hap1.append(None)
                    hap2.append(None)
        else:
            hap1.append(None)
            hap2.append(None)
    return hap1,hap2

def calculate_normal_haplotype_copynumber(het,normal_coverage1,normal_coverage2):
    """
    Calculate normal coverage for each pseudo-haplotype
    """
    hap1 = []
    hap2 = []
    for i in range(len(het)):
        bin1 = range(max(0,int(het["hap1_pos"][i])-75),min(len(normal_coverage1[str(het["hap1_contig"][i])]),int(het["hap1_pos"][i])+75))
        bin2 = range(max(0,int(het["hap2_pos"][i])-75),min(len(normal_coverage2[str(het["hap2_contig"][i])]),int(het["hap2_pos"][i])+75))
        normal1 = sum([normal_coverage1[str(het["hap1_contig"][i])][j] for j in bin1])
        normal2 = sum([normal_coverage2[str(het["hap2_contig"][i])][k] for k in bin2])
        if normal1 != 0 and normal2 != 0:
            cn1 = 2 * normal1/(normal1+normal2)
            cn2 = 2 - cn1
            hap1.append(cn1)
            hap2.append(cn2)
        elif normal1 == 0 and normal2 != 0:
            hap1.append(0)
            hap2.append(2)
        elif normal1 != 0 and normal2 == 0:
            hap1.append(2)
            hap2.append(0)
        elif normal1 == 0 and normal2 == 0:
            hap1.append(None)
            hap2.append(None)
    return hap1,hap2

def main():
    args = parse_args()
    if not args.OUTPUT_PATH:
        exit("Specify output directory using -o")
    else:
        output_path = args.OUTPUT_PATH

    if not args.NORMAL_1:
        exit("Specify fitst normal bam using -n1")
    else:
        normal1 = args.NORMAL_1

    if not args.NORMAL_2:
        exit("Specify second normal bam using -n2")
    else:
        normal2 = args.NORMAL_2

     if not args.TUMOR_1:
        exit("Specify fitst tumor bam using -t1")
    else:
        tumor1 = args.TUMOR_1

    if not args.TUMOR_2:
        exit("Specify second tumor bam using -n2")
    else:
        tumor2 = args.TUMOR_2

    if not args.POLY_PATH:
        exit("Specify path to polymorphic site tsv file using -p")
    else:
        poly_path = args.POLY_PATH

    if not args.CN_PATH:
        exit("Specify Path to combined copy number tsv file using -c")
    else:
        cn_path = args.CN_PATH

    if not args.PURITY:
        exit("Specify tumor purity using -u")
    else:
        purity = args.PURITY

    het = get_polymorphic_site(poly_path)

    normal_coverage1 = get_coverage(normal1)
    normal_coverage2 = get_coverage(normal2)
    normal_cn1,normal_cn2 = calculate_normal_haplotype_copynumber(het,normal_coverage1,normal_coverage2)

    tumor_coverage1 = get_coverage(tumor1)
    tumor_coverage2 = get_coverage(tumor2)
    BAF = calculate_baf(tumor_coverage1,tumor_coverage2,het)
    combined_cn = get_combined_copynumber(het,cn_path)
    tumor_cn1,tumor_cn2 = calculate_tumor_haplotype_copynumber(BAF,purity,combined_cn)

    df = pd.DataFrame({"pos":het["hg38_pos"].tolist(),"tumor1":tumor_cn1,"tumor2":tumor_cn2,"normal1":normal_cn1,"normal2":normal_cn2})
    df = df.dropna()
    df.to_csv(output_path,index=False)

if __name__ == "__main__":
    main()

