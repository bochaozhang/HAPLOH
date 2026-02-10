if [ "$#" -ne 17 ]; then
    echo ""
    echo "Usage:"
    echo "  $0 <workDir> <sample> <normal_HP1> <normal_HP2> <tumor_HP1> <tumor_HP2> \\"
    echo "     <HP1_coords> <HP2_coords> <window_size> <start> <end> <average_size> \\"
    echo "     <regions_to_analyze> <purity_csv> <min_length> <min_valid_points> <smooth_size>"
    echo ""
    echo "Arguments:"
    echo "  workDir"
    echo "      Working/output directory."
    echo "      All intermediate and final results will be generated under this directory."
    echo ""
    echo "  sample"
    echo "      Sample ID used as the prefix for all output files."
    echo ""
    echo "  normal_HP1"
    echo "      BAM file of the normal sample aligned to haplotype 1 assembly."
    echo ""
    echo "  normal_HP2"
    echo "      BAM file of the normal sample aligned to haplotype 2 assembly."
    echo ""
    echo "  tumor_HP1"
    echo "      BAM file of the tumor sample aligned to haplotype 1 assembly."
    echo ""
    echo "  tumor_HP2"
    echo "      BAM file of the tumor sample aligned to haplotype 2 assembly."
    echo ""
    echo "  HP1_coords"
    echo "      Coordinate mapping file of positions for haplotype 1"
    echo "      (haplotype assembly → hg38)."
    echo ""
    echo "  HP2_coords"
    echo "      Coordinate mapping file of positions for haplotype 2"
    echo "      (haplotype assembly → hg38)."
    echo ""
    echo "  window_size"
    echo "      Window size (bp) used for estimated ploidy calculation."
    echo ""
    echo "  start"
    echo "      Start coordinate of the analyzed region on hg38."
    echo ""
    echo "  end"
    echo "      End coordinate of the analyzed region on hg38."
    echo ""
    echo "  average_size"
    echo "      Window size (bp) for smoothing coverage in coverage plots."
    echo ""
    echo "  regions_to_analyze"
    echo "      Genomic sub-regions to analyze, provided as a Python-style list."
    echo "      Example: \"[(28903952,30915858),(30915859,31177951),(31177952,33268517)]\""
    echo ""
    echo "  purity_csv"
    echo "      CSV file containing tumor purity estimates for samples."
    echo ""
    echo "  min_length"
    echo "      Minimum deletion length (bp) to be considered."
    echo ""
    echo "  min_valid_points"
    echo "      Minimum number of valid data points required per window."
    echo ""
    echo "  smooth_size"
    echo "      Window size (bp) for smoothing the ploidy in coverage plots."
    echo ""
    exit 1
fi



# 解析输入参数
workDir=$1
sample=$2
normal_HP1=$3
normal_HP2=$4
tumor_HP1=$5
tumor_HP2=$6
HP1_coords=$7
HP2_coords=$8
window_size=$9
start=$10
end=$11
average_size=$12
regions_to_analyze=$13
purity_csv=$14
min_length=$15
min_valid_points=${16}
smooth_size=${17}

cd "$workDir"
samtools view ${normal_HP1} | cut -f1 | sort | uniq > ${sample}.normal.hap1_perfect_reads.list
samtools view ${normal_HP2} | cut -f1 | sort | uniq > ${sample}.normal.hap2_perfect_reads.list
# 找出共有的reads
comm -12 <(sort ${sample}.normal.hap1_perfect_reads.list) <(sort ${sample}.normal.hap2_perfect_reads.list) > common_normal_reads.list
# 找出特有reads
comm -23 <(sort ${sample}.normal.hap1_perfect_reads.list) <(sort common_normal_reads.list) > unique_normal_hap1.list
comm -23 <(sort ${sample}.normal.hap2_perfect_reads.list) <(sort common_normal_reads.list) > unique_normal_hap2.list

samtools view -b -N unique_normal_hap1.list ${normal_HP1} | samtools sort -@ 25 -O BAM - > ${sample}.normal.unique.perfect.hap1.bam && samtools index -@ 25 ${sample}.normal.unique.perfect.hap1.bam
samtools view -b -N unique_normal_hap2.list ${normal_HP2} | samtools sort -@ 25 -O BAM - > ${sample}.normal.unique.perfect.hap2.bam && samtools index -@ 25 ${sample}.normal.unique.perfect.hap2.bam
samtools view ${tumor_HP1} | cut -f1 | sort | uniq > ${sample}.tumor.hap1_perfect_reads.list
samtools view ${tumor_HP2} | cut -f1 | sort | uniq > ${sample}.tumor.hap2_perfect_reads.list

# 找出共有的reads
comm -12 <(sort ${sample}.tumor.hap1_perfect_reads.list) <(sort ${sample}.tumor.hap2_perfect_reads.list) > common_tumor_reads.list
# 找出特有reads
comm -23 <(sort ${sample}.tumor.hap1_perfect_reads.list) <(sort common_tumor_reads.list) > unique_tumor_hap1.list
comm -23 <(sort ${sample}.tumor.hap2_perfect_reads.list) <(sort common_tumor_reads.list) > unique_tumor_hap2.list
samtools view -b -N unique_tumor_hap1.list ${tumor_HP1} | samtools sort -@ 25 -O BAM - > ${sample}.tumor.unique.perfect.hap1.bam && samtools index -@ 25 ${sample}.tumor.unique.perfect.hap1.bam
samtools view -b -N unique_tumor_hap2.list ${tumor_HP2} | samtools sort -@ 25 -O BAM - > ${sample}.tumor.unique.perfect.hap2.bam && samtools index -@ 25 ${sample}.tumor.unique.perfect.hap2.bam

samtools depth -a ${normal_HP1} > ${sample}.normal.hap1.assembly.depth.txt
samtools depth -a ${sample}.normal.unique.perfect.hap2.bam > ${sample}.normal.hap2.assembly.depth.txt
samtools depth -a ${tumor_HP1} > ${sample}.tumor.hap1.assembly.depth.txt
samtools depth -a ${sample}.tumor.unique.perfect.hap2.bam > ${sample}.tumor.hap2.assembly.depth.txt


python 06_01_merge_bed_depth.py --bed_file ${HP1_coords} --depth_file ${sample}.normal.hap1.assembly.depth.txt --output_file ${sample}.normal.hap1.bed
python 06_01_merge_bed_depth.py --bed_file ${HP2_coords} --depth_file ${sample}.normal.hap2.assembly.depth.txt --output_file ${sample}.normal.hap2.uniq.bed
python 06_01_merge_bed_depth.py --bed_file ${HP1_coords} --depth_file ${sample}.tumor.hap1.assembly.depth.txt --output_file ${sample}.tumor.hap1.bed
python 06_01_merge_bed_depth.py --bed_file ${HP2_coords} --depth_file ${sample}.tumor.hap2.assembly.depth.txt --output_file ${sample}.tumor.hap2.uniq.bed
python 06_02_merge_bed_files.py --file_hap1 ${sample}.normal.hap1.bed --file_hap2 ${sample}.normal.hap2.uniq.bed --output_file normal.bed
python 06_02_merge_bed_files.py --file_hap1 ${sample}.tumor.hap1.bed --file_hap2 ${sample}.tumor.hap2.uniq.bed --output_file tumor.bed
python 06_03_merge_and_plot.py --start ${start} --end ${end} --regions_to_analyze ${regions_to_analyze}  --file_hp1 normal.bed --file_hp2 tumor.bed --output_file ${sample}.csv --cov_1_png ${sample}.t_n.png --cov_miss_png ${sample}.t_n_miss.png --cov_log_png  ${sample}.t_n_cov_log.png --hp1_col normal_cov --hp2_col tumor_cov --average_size ${average_size}

python 06_04_find_sparse_regions.py --vcf_file ${sample}.nooverlap.overlap.diff.txt --start ${start} --end ${end} --threshold 2 --min_region_size 3000 --output_file hom.region.txt

samtools depth -a ${sample}.normal.unique.perfect.hap1.bam > ${sample}.normal.hap1.uniq.depth.txt
samtools depth -a ${sample}.normal.unique.perfect.hap2.bam > ${sample}.normal.hap2.uniq.depth.txt
samtools depth -a ${sample}.tumor.unique.perfect.hap1.bam > ${sample}.tumor.hap1.uniq.depth.txt
samtools depth -a ${sample}.tumor.unique.perfect.hap2.bam > ${sample}.tumor.hap2.uniq.depth.txt


python 06_01_merge_bed_depth.py --bed_file ${HP1_coords} --depth_file ${sample}.normal.hap1.uniq.depth.txt --output_file ${sample}.normal.hap1.uniq.bed
python 06_01_merge_bed_depth.py --bed_file ${HP1_coords} --depth_file ${sample}.tumor.hap1.uniq.depth.txt --output_file ${sample}.tumor.hap1.uniq.bed
python 06_01_merge_bed_depth.py --bed_file ${HP2_coords} --depth_file ${sample}.normal.hap2.uniq.depth.txt --output_file ${sample}.normal.hap2.uniq.bed
python 06_01_merge_bed_depth.py --bed_file ${HP2_coords} --depth_file ${sample}.tumor.hap2.uniq.depth.txt --output_file ${sample}.tumor.hap2.uniq.bed

python 06_05_merge_and_plot_remove_region.py --start ${start} --end ${end} --regions_to_analyze ${regions_to_analyze}  --file_hp1 ${sample}.tumor.hap1.uniq.bed --file_hp2 ${sample}.tumor.hap2.uniq.bed --output_file ${sample}.tumor.csv --cov_1_png ${sample}.tumor.png --cov_miss_png ${sample}.tumor_miss.png --cov_log_png  ${sample}.tumor_cov_log.png --hp1_col hap1_tumor_cov --hp2_col hap2_tumor_cov --average_size ${average_size}
python 06_05_merge_and_plot_remove_region.py --start ${start} --end ${end} --regions_to_analyze ${regions_to_analyze}  --file_hp1 ${sample}.normal.hap1.uniq.bed --file_hp2 ${sample}.normal.hap2.uniq.bed --output_file ${sample}.normal.csv --cov_1_png ${sample}.normal.png --cov_miss_png ${sample}.normal_miss.png --cov_log_png  ${sample}.normal_cov_log.png --hp1_col hap1_normal_cov --hp2_col hap2_normal_cov --average_size ${average_size}



python 06_06_ploidy_sum_pic_10kbp.py --sample_id ${sample} --purity_csv ${purity_csv} --window_size ${window_size}  --start ${start} --end ${end} --regions_to_analyze ${regions_to_analyze}



######two_deletion
python 06_07_two_deletion_sum_low.py --sample  ${sample}  --start ${start} --end ${end} --regions_to_analyze ${regions_to_analyze} --min_length ${min_length} --min_valid_points ${min_valid_points} --window_size ${window_size}

python 06_08_two_deletion_sum.pic.py --sample  ${sample} --window_size ${window_size} --smooth_size ${smooth_size} --purity_csv /home/users/bczhang/HLALOH/tumor_purity.csv --start ${start} --end ${end} --regions_to_analyze ${regions_to_analyze}



python 06_09_filter_copynumber.py \
    --id ${sample} \
    --region_file ${sample}.low.sum.$((window_size / 500))kbp.txt \
    --copynumber_file /home/users/bczhang/HLALOH/${sample}/${sample}.mhc.copynumber_loess_$((window_size / 500))kbp_sum.csv \
    --output_file ${sample}.mhc.copynumber_loess_$((window_size / 500))kbp_sum_no_del.csv





