#!/bin/bash
#
#
# Replace ACCOUNT with your account name before submitting.
#
#SBATCH --account=cgl      	 # Replace ACCOUNT with your group account name
#SBATCH --job-name=assign_reads    # The job name
#SBATCH --time=0-6:00            # The time the job will take to run in D-HH:MM
#SBATCH --mem=20G
#SBATCH -o "batch-logs/7-assign-reads-logs/%u-%x-%j.out"

INPUT_DIR=$1
FILE_LIST=$2
EXON_BED=$3
GENE_BED=$4
SCRIPTS_DIR=$5
OUTPUT_DIR=$6
BIN=$7
DATAMASH=$8

BEDTOOLS=$BIN/bedtools

cat $FILE_LIST | while read FILE; do
    PCR_WELL=`basename "$FILE" .bed`

    $BEDTOOLS map \
        -a $INPUT_DIR/$FILE \
        -b $EXON_BED \
        -nonamecheck -s -f 0.95 -c 7 -o distinct -delim '|' \
    | $BEDTOOLS map \
        -a - -b $GENE_BED \
        -nonamecheck -s -f 0.95 -c 4 -o distinct -delim '|' \
    | sort -k4,4 -k2,2n -k3,3n -S 3G \
    | $DATAMASH \
        -g 4 first 1 first 2 last 3 first 5 first 6 collapse 7 collapse 8 \
    | $SCRIPTS_DIR/8-assign-reads-to-genes.py $GENE_BED \
    >$OUTPUT_DIR/$PCR_WELL

    echo "Processed $FILE"
done

