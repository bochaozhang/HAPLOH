########################## file information #############################
# Author: Bochao Zhang
# Date: 2023-2-23
# Description: 
# Select and filter 10x haplotype assembly based on N's per 10k bases and scafold length 
#########################################################################

import argparse
import os
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser(description = "A Python script for filtering and selecting of assemblies")
    parser.add_argument("-i", help="Input directory (REQUIRED)", dest="INPUT_DIR")
    parser.add_argument("-o", help="Output directory", dest="OUTPUT_DIR")
    parser.add_argument("-q", help="Quast results directory", dest="QUAST_DIR")
    parser.add_argument("-p", help="Quast script path", dest="QUAST_PATH")
    parser.add_argument("-c", help="Minimum cutoff of scafold length (default 50k)", dest="CUTOFF")
    return parser.parse_args()

def run_quast(input_dir,output_dir,quast_path):
    """
    Run Quast to evaluate assemblies
    """
    script_path = os.path.abspath(__file__)
    cmd = "bash {script_path}/01_quast.sh -p {quast_path} -r {resource_dir} -i {input_dir} -o {output_dir}".format( \
        script_path = os.path.abspath(__file__), \
        quast_path = quast_path, \
        resource_dir = os.path.join(os.path.abspath(__file__),"data"), \
        input_dir = input_dir, \
        output_dir = output_dir)
    os.system(cmd)

def select_assembly(quast_dir,haploid):
    """
    Select pseudo haplotype assembly with the longest effective length
    """
    infile = os.path.join(quast_dir,"quast_hp%d_results" % haploid,"report.tsv")
    df = pd.read_csv(infile,sep="\t")
    tag1 = "HP%d.1" % haploid
    tag2 = "HP%d.2" % haploid
    length1 = int(df[tag1][df["Assembly"]=="Total length (>= 50000 bp)"])
    length2 = int(df[tag2][df["Assembly"]=="Total length (>= 50000 bp)"])
    num_N1 = float(df[tag1][df["Assembly"]=="# N's per 100 kbp"])
    num_N2 = float(df[tag2][df["Assembly"]=="# N's per 100 kbp"])
    effective_len1 = length1 - length1*(num_N1/100000)
    effective_len2 = length2 - length2*(num_N2/100000)
    if effective_len1 >= effective_len2:
        selected_fasta = tag1 + ".fasta"
    else:
        selected_fasta = tag2 + ".fasta"
    return selected_fasta

def filter_assembly(quast_dir,output_dir,cutoff,haploid):
    '''
    Filter scafolds based on length
    '''
    selected_fasta = select_assembly(quast_dir,haploid)
    infile = os.path.join(quast_dir,selected_fasta)
    lines = open(infile).readlines()
    p = [i for i in range(len(lines)) if lines[i].startswith(">")]
    p += [len(lines)]
    outfile = "HP%d.%d.filtered.fasta" % (haploid,cutoff)
    o = open(os.path.join(output_dir,outfile),"w")
    for i in range(len(p)-1):
        seq = ""
        for j in range(p[i]+1,p[i+1]):
            seq += lines[j].strip()
        if len(seq) >= cutoff:
            o.write("%s%s\n" % (lines[p[i]],seq))
    o.close()

def main():
    args = parse_args()
    if not args.INPUT_DIR:
        exit("Specify input directory using -i")
    else:
        input_dir = args.INPUT_DIR

    if not args.OUTPUT_DIR:
        output_dir = args.INPUT_DIR
    else:
        output_dir = args.OUTPUT_DIR

    if not args.CUTOFF:
        cutoff = 50000
    else:
        cutoff = int(args.CUTOFF)

    if not args.QUAST_DIR:
        if not args.QUAST_PATH:
            exit("Specify Quast results folder using -q or Quast script path using -p")
        else:
            quast_path = args.QUAST_PATH
            quast_dir = output_dir
            run_quast(input_dir,quast_dir,quast_path)        
    else:
        quast_dir = args.QUAST_DIR

    filter_assembly(quast_dir,output_dir,cutoff,1)
    filter_assembly(quast_dir,output_dir,cutoff,2)

if __name__ == "__main__":
    main()
                                      
