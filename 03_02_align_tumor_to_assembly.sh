########################## file information #############################
# Author: Bochao Zhang
# Date: 2024-9-26
# Description:
# Align tumor reads to assembly
#########################################################################

set -e

## Read input arguments
while getopts p:w:i1:i2: flag; do
    case "${flag}" in
        w) WORKING_DIR=${OPTARG};;
        p) RAZERS_PATH=${OPTARG};;
		g) GATK_JAR=${OPTARG};;
		s) SAMTOOLS_PATH=${OPTARG};;
		r1) REF1=${OPTARG};;
		r2) REF2=${OPTARG};;
        i1) FASTQ1=${OPTARG};;
        i2) FASTQ2=${OPTARG};;
    esac
done

if ! [ -d "$WORKING_DIR" ]; then
    echo "$WORKING_DIR isn't a direcory"
    exit 1
fi

if ! [ -f "$RAZERS_PATH" ]; then
    echo "$RAZERS_PATH isn't a valid path"
    exit 1
fi

if ! [ -f "$SAMTOOLS_PATH" ]; then
    echo "$SAMTOOLS_PATH isn't a valid path"
    exit 1
fi

if ! [ -f "$GATK_JAR" ]; then
    echo "$GATK_JAR isn't a valid jar file"
    exit 1
fi

if [![ -f "$REF1" ]] || [![ -f "$REF2" ]]; then
	echo "Reference file(s)  doesn't exist"
    exit 1
fi
if [![ -f "$FASTQ1" ]] || [![ -f "$FASTQ2" ]]; then
    echo "Input FASTQ file(s) not exist"
    exit 1
fi

gatk="java -jar $GATK_JAR"

rzr_opt="--percent-identity 99 --no-gaps --distance-range 0 --max-hits 5 --library-length 440 --library-error 400 --sort-order 1 --thread-count 8"

## Align to hap1
$RAZERS_PATH $rzr_opt -o $WORKING_DIR/tumor.to.assembly1.bam $REF1 $FASTQ1 $FASTQ1
$SAMTOOLS_PATH index $WORKING_DIR/tumor.to.assembly1.bam

$gatk MarkDuplicates \
    -I $WORKING_DIR/tumor.to.assembly1.bam \
    -O $WORKING_DIR/tumor.to.assembly1.dedup.bam \
    -M $WORKING_DIR/tumor.to.assembly1.dedup.txt \
    --REMOVE_DUPLICATES true
$SAMTOOLS_PATH index $WORKING_DIR/tumor.to.assembly1.dedup.bam

## Align to hap2
$RAZERS_PATH $rzr_opt -o $WORKING_DIR/tumor.to.assembly2.bam $REF2 $FASTQ1 $FASTQ1
$SAMTOOLS_PATH index $WORKING_DIR/tumor.to.assembly2.bam

$gatk MarkDuplicates \
    -I $WORKING_DIR/tumor.to.assembly2.bam \
    -O $WORKING_DIR/tumor.to.assembly2.dedup.bam \
    -M $WORKING_DIR/tumor.to.assembly2.dedup.txt \
    --REMOVE_DUPLICATES true
$SAMTOOLS_PATH index $WORKING_DIR/tumor.to.assembly2.dedup.bam


# cleanup
if [ -f $WORKING_DIR/tumor.to.assembly1.dedup.bam ]; then
	rm $WORKING_DIR/tumor.to.assembly1.bam
	rm $WORKING_DIR/tumor.to.assembly1.bam.bai
	rm $WORKING_DIR/tumor.to.assembly1.dedup.txt
fi

if [ -f $WORKING_DIR/tumor.to.assembly2.dedup.bam ]; then
	rm $WORKING_DIR/tumor.to.assembly2.bam
	rm $WORKING_DIR/tumor.to.assembly2.bam.bai
	rm $WORKING_DIR/tumor.to.assembly2.dedup.txt
fi

