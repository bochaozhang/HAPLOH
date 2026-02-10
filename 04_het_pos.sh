#!/bin/bash
if [ "$#" -ne 5 ]; then
    echo ""
    echo "Usage:"
    echo "  $0 <workDir> <HP1> <HP2> <sample> <ref>"
    echo ""
    echo "Arguments:"
    echo "  workDir   Output working directory. All intermediate and final results"
    echo "            will be generated under this directory."
    echo ""
    echo "  HP1       FASTA file of haplotype 1 assembly (HP1)."
    echo ""
    echo "  HP2       FASTA file of haplotype 2 assembly (HP2)."
    echo ""
    echo "  sample    Sample name used as the prefix for output files."
    echo ""
    echo "  ref       Reference FASTA of the hg38 MHC region"
    echo "            (chr6:28903952-33268517)."
    echo ""
    exit 1
fi


# 解析输入参数
workDir=$1
HP1=$2
HP2=$3
sample=$4
ref=$5
cd "$workDir"
mkdir -p HP1_dnadiff
cd HP1_dnadiff
####比对
dnadiff ${ref} ${HP1}
cd ..
mkdir -p HP2_dnadiff
cd HP2_dnadiff
dnadiff ${ref} ${HP2}
cd ..
####重叠region根据idendity进行保留，保证1-1
python 04_01_remove_overlap.py ${workDir} ${HP1} ${HP2} ${ref}

mkdir -p HP1_dnadiff_nooverlap
cd HP1_dnadiff_nooverlap
show-coords -THrcl ../HP1_dnadiff/out_noOverlap.1delta > out_noOverlap.1coords
dnadiff -d ../HP1_dnadiff/out_noOverlap.1delta
####remap到hg38
for i in `grep '>' ${HP1} | cut -d' ' -f1 | cut -d'>' -f2`;do show-aligns  ../HP1_dnadiff/out_noOverlap.1delta  chr6_28903952-33268517 $i | perl -le 'while(<>){chomp;next if /^(\s+)|^$/;next if $.<=3;@F=split/\s+/,$_;if(/-/){print}else{print ;};}' | perl -le '$n=$m=0;while(<>){chomp;chomp($G=<>);@A=@B=@AA=@BB=();if (/^--/ or /^=/){if(/Alignments between (\S+) and (\S+)/){$ref=$1;$denovo=$2;};if($G=~/BEGIN alignment \[ (.*)\]/){$M=$1;};@MM=split/\s+/,$M;print $_;print $G;}else{@A=split/\s+/,$_;@AA=split//,$A[1];@B=split/\s+/,$G;@BB=split//,$B[1];};$a=$A[0];$b=$B[0];for my $i (0..$#AA){if ($AA[$i] ne "."){if($MM[0] eq "+1"){$aa=$a+$i-$n}elsif($MM[0] eq "-1"){$aa=$a-$i+$n};}else{$aa ="N";$n+=1;};if ($BB[$i] ne "."){if($MM[5] eq "+1"){$bb=$b+$i-$m}elsif($MM[5] eq "-1"){$bb=$b-$i+$m};}else{$bb ="N";$m+=1;};if($aa eq "N"){$chr6 = "N"}else{$chr6=$aa+28903952-1};print "chr6\t$ref\t$denovo\t$chr6\t$aa\t$bb\t$AA[$i]\t$BB[$i]\t$MM[0]\t$MM[5]\t\t---$_\t$G";};$n=$m=0; ;}' >> ${sample}.HP1.new.filter.delta.coords.align.detail.infor;done
perl -le 'print "hg38\tmhc\tdenovo\thg38POS\tmhcPOS\tdenovoPOS";while(<>){chomp;@F=split/\s+/,$_;next if /^--|^=/;if($F[9] eq "+1"){print join("\t",@F[0..5]);}else{$m=$F[5]-1;print join("\t",@F[0..4]),"\t$m";};}'  ${sample}.HP1.new.filter.delta.coords.align.detail.infor|perl -le 'while(<>){chomp;@F=split/\s+/,$_;push @{$hash{$F[3]}},$_;}END{foreach my $key (keys %hash){print "$key\t",join("\t",@{$hash{$key}});};}'|sort -k1|perl -lane  'next if /hg38|\tN/;print join("\t",@F[1..6]) if $#F<=6'|perl -le 'while(<>){chomp;@F=split/\s+/,$_;push @{$hash{"$F[2]\t$F[5]"}},$_;}END{foreach my $key (keys %hash){print "$key\t",join("\t",@{$hash{$key}});};}'|perl -lane  'next if /hg38|\tN/;print join("\t",@F[2..7]) if $#F<=7'|sort -k4 >${sample}.HP1.new.filter.delta.coords.align.detail.infor.bed.unique.pos.new.bed



python 04_02_process_regions.py ${sample}.HP1 ${workDir}/HP1_dnadiff_nooverlap


###转化为hg38的坐标，之前是截取的chr6_28903952-33268517，会从1开始
awk -F'\t' -v OFS='\t' '{
    split($12,a,"-"); split(a[1],b,"_"); $1=$1+b[2]-1; $2=$2+b[2]-1;
    print $0
}' ${sample}.HP1.filtered.region.txt > ${sample}.HP1.filtered.region.hg38.txt
awk -F'\t' -v OFS='\t' '{
    split($11,a,"-"); split(a[1],b,"_"); $1=$1+b[2]-1;
    $4=$12"-"$4;
    print $0
}' out.snps > together.snp_pos.txt
#####去除插入缺失的一些位点
awk '($2 != "." && $3 != ".")' together.snp_pos.txt | sort -k1,1 -n > sorted_filtered_together.snp.txt
awk 'BEGIN {OFS="\t"} {print "chr6", $1, ".", $2, $3, $4}' sorted_filtered_together.snp.txt > together.vcf

###构建vcf文件
python 04_03_process_vcf_and_region.py ${sample}.HP1 ${workDir}/HP1_dnadiff_nooverlap

cd ../

mkdir -p HP2_dnadiff_nooverlap
cd HP2_dnadiff_nooverlap
show-coords -THrcl ../HP2_dnadiff/out_noOverlap.1delta > out_noOverlap.1coords
dnadiff -d ../HP2_dnadiff/out_noOverlap.1delta

for i in `grep '>' ${HP2} | cut -d' ' -f1 | cut -d'>' -f2`;do show-aligns  ../HP2_dnadiff/out_noOverlap.1delta  chr6_28903952-33268517 $i | perl -le 'while(<>){chomp;next if /^(\s+)|^$/;next if $.<=3;@F=split/\s+/,$_;if(/-/){print}else{print ;};}' | perl -le '$n=$m=0;while(<>){chomp;chomp($G=<>);@A=@B=@AA=@BB=();if (/^--/ or /^=/){if(/Alignments between (\S+) and (\S+)/){$ref=$1;$denovo=$2;};if($G=~/BEGIN alignment \[ (.*)\]/){$M=$1;};@MM=split/\s+/,$M;print $_;print $G;}else{@A=split/\s+/,$_;@AA=split//,$A[1];@B=split/\s+/,$G;@BB=split//,$B[1];};$a=$A[0];$b=$B[0];for my $i (0..$#AA){if ($AA[$i] ne "."){if($MM[0] eq "+1"){$aa=$a+$i-$n}elsif($MM[0] eq "-1"){$aa=$a-$i+$n};}else{$aa ="N";$n+=1;};if ($BB[$i] ne "."){if($MM[5] eq "+1"){$bb=$b+$i-$m}elsif($MM[5] eq "-1"){$bb=$b-$i+$m};}else{$bb ="N";$m+=1;};if($aa eq "N"){$chr6 = "N"}else{$chr6=$aa+28903952-1};print "chr6\t$ref\t$denovo\t$chr6\t$aa\t$bb\t$AA[$i]\t$BB[$i]\t$MM[0]\t$MM[5]\t\t---$_\t$G";};$n=$m=0; ;}' >> ${sample}.HP2.new.filter.delta.coords.align.detail.infor;done
perl -le 'print "hg38\tmhc\tdenovo\thg38POS\tmhcPOS\tdenovoPOS";while(<>){chomp;@F=split/\s+/,$_;next if /^--|^=/;if($F[9] eq "+1"){print join("\t",@F[0..5]);}else{$m=$F[5]-1;print join("\t",@F[0..4]),"\t$m";};}'  ${sample}.HP2.new.filter.delta.coords.align.detail.infor|perl -le 'while(<>){chomp;@F=split/\s+/,$_;push @{$hash{$F[3]}},$_;}END{foreach my $key (keys %hash){print "$key\t",join("\t",@{$hash{$key}});};}'|sort -k1|perl -lane  'next if /hg38|\tN/;print join("\t",@F[1..6]) if $#F<=6'|perl -le 'while(<>){chomp;@F=split/\s+/,$_;push @{$hash{"$F[2]\t$F[5]"}},$_;}END{foreach my $key (keys %hash){print "$key\t",join("\t",@{$hash{$key}});};}'|perl -lane  'next if /hg38|\tN/;print join("\t",@F[2..7]) if $#F<=7'|sort -k4 >${sample}.HP2.new.filter.delta.coords.align.detail.infor.bed.unique.pos.new.bed
python 04_02_process_regions.py ${sample}.HP2 ${workDir}/HP2_dnadiff_nooverlap

awk -F'\t' -v OFS='\t' '{
    split($12,a,"-"); split(a[1],b,"_"); $1=$1+b[2]-1; $2=$2+b[2]-1;
    print $0
}' ${sample}.HP2.filtered.region.txt > ${sample}.HP2.filtered.region.hg38.txt
awk -F'\t' -v OFS='\t' '{
    split($11,a,"-"); split(a[1],b,"_"); $1=$1+b[2]-1;
    $4=$12"-"$4;
    print $0
}' out.snps > together.snp_pos.txt

awk '($2 != "." && $3 != ".")' together.snp_pos.txt | sort -k1,1 -n > sorted_filtered_together.snp.txt

awk 'BEGIN {OFS="\t"} {print "chr6", $1, ".", $2, $3, $4}' sorted_filtered_together.snp.txt > together.vcf

###构建vcf文件
python 04_03_process_vcf_and_region.py ${sample}.HP2 ${workDir}/HP2_dnadiff_nooverlap

cd ..


####获取hap1、hap2同时能remap到hg38的region
python 04_04_find_overlap_and_pos.py ${workDir} ${sample}

####根据共同的region进行过
cd HP1_dnadiff_nooverlap
python 04_05_filter_vcf_by_overlap.py ${sample} ${workDir}/HP1_dnadiff_nooverlap
cd ../HP2_dnadiff_nooverlap
python 04_05_filter_vcf_by_overlap.py ${sample} ${workDir}/HP2_dnadiff_nooverlap
cd ..


####构建vcf文件，根据1-1remap的情况，会过滤到一个突变，另一个插入缺失的位点
python 04_06_overlap.diff.snp.py  ${workDir}/HP1_dnadiff_nooverlap/together_all.overlap.filter.vcf ${workDir}/HP2_dnadiff_nooverlap/together_all.overlap.filter.vcf  ${workDir}/HP1_dnadiff_nooverlap/${sample}.HP1.new.filter.delta.coords.align.detail.infor.bed.unique.pos.new.bed ${workDir}/HP2_dnadiff_nooverlap/${sample}.HP2.new.filter.delta.coords.align.detail.infor.bed.unique.pos.new.bed ${workDir}/${sample}.nooverlap.overlap.diff.txt
python 04_07_overlap.all.diff.snp.py ${workDir}/HP1_dnadiff_nooverlap/together_all.overlap.filter.vcf ${workDir}/HP2_dnadiff_nooverlap/together_all.overlap.filter.vcf  ${workDir}/HP1_dnadiff_nooverlap/${sample}.HP1.new.filter.delta.coords.align.detail.infor.bed.unique.pos.new.bed ${workDir}/HP2_dnadiff_nooverlap/${sample}.HP2.new.filter.delta.coords.align.detail.infor.bed.unique.pos.new.bed ${workDir}/${sample}.nooverlap.all.overlap.diff.txt




