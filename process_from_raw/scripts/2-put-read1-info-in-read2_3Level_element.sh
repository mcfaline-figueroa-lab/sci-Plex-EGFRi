#!/bin/bash
#
#
# Replace ACCOUNT with your account name before submitting.
#
#SBATCH --account=cgl      	 # Replace ACCOUNT with your group account name
#SBATCH --job-name=Read1Read2     # The job name
#SBATCH --time=0-6:00            # The time the job will take to run in D-HH:MM
#SBATCH --mem=25G         # The memory the job will use per cpu core
#SBATCH -o "batch-logs/1-R1-info-R2-logs/%u-%x-%j.out"


INPUT_DIR=$1
R1_FILE_LIST=$2
SCRIPTS_DIR=$3
RT_OLIGO_LIST=$4
LIG_OLIGO_LIST=$5
INDEXING_KEY=$6
OUTPUT_DIR=$7

BATCH_ID=`basename "$R1_FILE_LIST"`

cat $R1_FILE_LIST | while read R1_FILE_RAW; do
    R1_FILE="${R1_FILE_RAW}/${R1_FILE_RAW}_R1.fastq.gz"
    R2_FILE=`echo "$R1_FILE" | sed 's/R1/R2/'`
    PCR_COMBO=`echo "$R1_FILE_RAW"`

    paste \
        <(gunzip <$INPUT_DIR/$R1_FILE) \
        <(gunzip <$INPUT_DIR/$R2_FILE) \
    | awk -f $SCRIPTS_DIR/2-put-read1-info-in-read2_3Level.awk -v PCR_COMBO="$PCR_COMBO" \
        $RT_OLIGO_LIST $LIG_OLIGO_LIST $INDEXING_KEY - \
    | gzip >$OUTPUT_DIR/$PCR_COMBO.fastq.gz

    echo "Processed $PCR_COMBO"
done 

