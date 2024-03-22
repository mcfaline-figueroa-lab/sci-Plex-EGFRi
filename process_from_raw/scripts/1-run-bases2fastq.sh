#!/bin/bash
#
#
# Replace ACCOUNT with your account name before submitting.
#
#SBATCH --account=cgl      	 # Replace ACCOUNT with your group account name
#SBATCH --job-name=bases2fastq   # The job name
#SBATCH --time=0-8:00            # The time the job will take to run in D-HH:MM
#SBATCH -c 24                    # The number of cores requested
#SBATCH -o "batch-logs/%u-%x-%j.out"

BIN=$1
SAMPLE_SHEET=$2
INPUT_DIR=$3
OUTPUT_DIR=$4

$BIN/bases2fastq --num-threads 24 \
                 -r $SAMPLE_SHEET \
                    $INPUT_DIR    \
                    $OUTPUT_DIR

