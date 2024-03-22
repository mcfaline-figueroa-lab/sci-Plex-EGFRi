#!/bin/bash
#
#
# Replace ACCOUNT with your account name before submitting.
#
#SBATCH --account=cgl      	 # Replace ACCOUNT with your group account name
#SBATCH --job-name=sam_sort     # The job name
#SBATCH -c 12                     # The number of cores to request
#SBATCH --time=0-4:00            # The time the job will take to run in D-HH:MM
#SBATCH --mem-per-cpu=6G
#SBATCH -o "batch-logs/4-sam-sort-logs/%u-%x-%j.out"


INPUT_DIR=$1
FILE_LIST=$2
OUTPUT_DIR=$3
BIN=$4

cat $FILE_LIST | while read FILE; do
    SAMPLE=`basename "$FILE" .Aligned.out.bam`

    $BIN/samtools view -bh -q 30 -F 4 $INPUT_DIR/$FILE \
    | $BIN/samtools sort -@ 12 - \
    >$OUTPUT_DIR/$SAMPLE.bam

    echo "Processed $FILE"
done

