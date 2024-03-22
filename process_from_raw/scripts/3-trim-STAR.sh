#!/bin/bash
#
#
# Replace ACCOUNT with your account name before submitting.
#
#SBATCH --account=cgl      	 # Replace ACCOUNT with your group account name
#SBATCH --job-name=TrimAAA     # The job name
#SBATCH --time=0-3:00            # The time the job will take to run in D-HH:MM
#SBATCH --mem=50G         # The memory the job will use
#SBATCH -o "batch-logs/2-trim-logs/%u-%x-%j.out"

INPUT_DIR=$1
FILE_LIST=$2
OUTPUT_DIR=$3
BIN=$4

cat $FILE_LIST | while read FILE; do
    $BIN/trim_galore $INPUT_DIR/$FILE -a AAAAAAAA --three_prime_clip_R1 1 --gzip -o $OUTPUT_DIR --path_to_cutadapt $BIN/cutadapt
    echo "Processed $FILE"
done
