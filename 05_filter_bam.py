########################## file information #############################
# Author: Bochao Zhang
# Date: 2024-2-11
# Description:
# Filter assembly referencing bam, removing reads: 1) mapped to multiple locations, or 2) mapped to multiple scaffolds, or 3) have more than one mismatch
#########################################################################

import os,pysam
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser(description = "A Python script for filtering assembly referencing bam")
    parser.add_argument("-b", help="Input bam (REQUIRED)", dest="INPUT_BAM")
    parser.add_argument("-o", help="Output directory", dest="OUTPUT_DIR")
    parser.add_argument("-p", help="Path to polymorphic site tsv file (REQUIRED)", dest="POLY_PATH")
    return parser.parse_args()

def get_polymorphic_site(poly_path):
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

def filterBam(het,input_bam,output_dir):
    """
    Filter bam, removing reads with more than one mismatch, one mismatch but on polymorphic site, or not properly paired
    """
    samfile = pysam.AlignmentFile(input_bam,"rb")
    reads = list(samfile.fetch())
    read_ref = {}
    for read in reads:
        ADD = False
        MD = read.get_tag("MD")
        if MD == "150":
            ADD = True
        else:
            mismatch = [c for c in MD if c in ["A","C","G","T","N"]]
            if len(mismatch) == 1:
                mismatch = mismatch[0]
                mismatch_pos = read.reference_start + int(MD.split(mismatch)[0])
                if hap == 1:
                    if sum((het["hap1_contig"]==read.reference_name) & (het["hap1_pos"]==mismatch_pos)) == 0:
                        ADD = True
                elif hap == 2:
                    if sum((het["hap2_contig"]==read.reference_name) & (het["hap2_pos"]==mismatch_pos)) == 0:
                        ADD = True
        if ADD:
            try:
                read_ref[read.query_name].append(read.reference_name)
            except KeyError:
                read_ref[read.query_name] = [read.reference_name]
    qualified_reads = set([key for key in read_ref.keys() if len(set(read_ref[key])) == 1 and len(read_ref[key]) == 2])
    output_bam = os.path.basename(input_bam).replace(".bam","_filtered.bam")
    filtered_sam = pysam.AlignmentFile(os.path.join(output_dir,output_bam),"wb",template=samfile)
    for read in reads:
        if read.query_name in qualified_reads:
            filtered_sam.write(read)
    samfile.close()
    filtered_sam.close()

def main():
    args = parse_args()
    if not args.INPUT_BAM:
        exit("Specify input bam using -b")
    else:
        input_bam = args.INPUT_BAM

    if not args.OUTPUT_DIR:
        output_dir = os.path.dirname(input_dir)
    else:
        output_dir = args.OUTPUT_DIR
    
    het = get_polymorphic_site(poly_path)
    filterBam(het,input_bam,output_dir)
       
if __name__ == "__main__":
    main()
