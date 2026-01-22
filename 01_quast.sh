########################## file information #############################
# Author: Bochao Zhang
# Date: 2023-2-13
# Description:
# Run Quast to evaluate assembly quality
#########################################################################

set -e

## Read input arguments
while getopts p:r:i:o: flag; do
	case "${flag}" in
		p) QUAST_PATH=${OPTARG};;
		r) RESOURCE_DIR=${OPTARG};;
		i) INPUT_DIR=${OPTARG};;
		o) OUTPUT_DIR=${OPTARG};;
	esac
done

if ! [ -d "$QUAST_PATH" ]; then
	echo "$QUAST_PATH isn't a valid path"
	exit 1
fi

if ! [ -d "$RESOURCE_DIR" ]; then
	echo "$RESOURCE_DIR isn't a direcory"
	exit 1
fi

if ! [ -d "$INPUT_DIR" ]; then
	echo "$INPUT_DIR isn't a direcory"
	exit 1
fi

if ! [ -d "$OUTPUT_DIR" ]; then
	echo "$OUTPUT_DIR isn't a direcory"
	exit 1
fi

## Check if assembly files exist
if [![ -f $INPUT_DIR/HP1.1.fasta ]] || [![ -f $INPUT_DIR/HP1.2.fasta ]] || [![ -f $INPUT_DIR/HP2.1.fasta ]] || [![ -f $INPUT_DIR/HP2.2.fasta ]]; then
	echo "Please name your input assembly files to HP1.1.fasta, HP1.2.fasta, HP2.2.fasta, HP2.2.fasta respectively"
	exit 1
fi

## Evaluate assemblies
$QUAST_PATH/quast.py \
    -o $OUTPUT_DIR/quast_hp1_results \
    -r $RESOURCE_DIR/hg38.chr6.fasta \
    -g $RESOURCE_DIR/hla_genes.bed \
    $INPUT_DIR/HP1.1.fasta $INPUT_DIR/HP1.2.fasta
$QUAST_PATH/quast.py \
    -o $OUTPUT_DIR/quast_hp2_results \
    -r $RESOURCE_DIR/hg38.chr6.fasta \
    -g $RESOURCE_DIR/hla_genes.bed \
    $INPUT_DIR/HP2.1.fasta $INPUT_DIR/HP2.2.fasta
