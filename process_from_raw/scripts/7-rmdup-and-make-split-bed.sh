#!/bin/bash
#
#
# Replace ACCOUNT with your account name before submitting.
#
#SBATCH --account=cgl      	 # Replace ACCOUNT with your group account name
#SBATCH --job-name=make_bed      # The job name
#SBATCH --time=0-4:00            # The time the job will take to run in D-HH:MM
#SBATCH --mem=10G
#SBATCH -o "batch-logs/6-make-bed-logs/%u-%x-%j.out"

INPUT_DIR=$1
FILE_LIST=$2
SCRIPTS_DIR=$3
OUTPUT_DIR=$4
BIN=$5

SAMTOOLS=$BIN/samtools
BEDTOOLS=$BIN/bedtools

# this makes "sort" case sensitive
export LC_ALL=C
START_TIME=$SECONDS
cat $FILE_LIST | while read FILE; do
    PCR_WELL=`basename "$FILE" .bam`

    $SAMTOOLS view -h $INPUT_DIR/$FILE \
    | awk -f $SCRIPTS_DIR/7-rmdup.awk \
    | $SAMTOOLS view -bh \
    | $BEDTOOLS bamtobed -i - -split \
    | sort -k1,1 -k2,2n -k3,3n -S 3G \
    >$OUTPUT_DIR/$PCR_WELL.bed

    ELAPSED_TIME=$(($SECONDS - $START_TIME))
    echo "Processed $FILE in $ELAPSED_TIME seconds"

done

