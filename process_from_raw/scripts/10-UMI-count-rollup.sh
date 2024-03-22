#!/bin/bash
#
#
# Replace ACCOUNT with your account name before submitting.
#
#SBATCH --account=cgl      	 # Replace ACCOUNT with your group account name
#SBATCH --job-name=UMI_rollup    # The job name
#SBATCH --time=0-08:00            # The time the job will take to run in D-HH:MM
#SBATCH --mem=4G
#SBATCH -o "batch-logs/10-UMI-rollup-logs/%u-%x-%j.out"


INPUT_DIR=$1
FILE_LIST=$2
OUTPUT_DIR=$3
DATAMASH=$4

cat $FILE_LIST | while read FILE; do
    awk '$3 == "exonic" || $3 == "intronic" {
        split($1, arr, "|");
        printf "%s|%s_%s_%s\t%s\n",
            arr[2], arr[3], arr[4], arr[5], $2;
    }' $INPUT_DIR/$FILE \
    | sort -k1,1 -k2,2 -S 2G \
    | $DATAMASH -g 1,2 count 2 \
    >$OUTPUT_DIR/$FILE

    echo "Processed $FILE"
done

