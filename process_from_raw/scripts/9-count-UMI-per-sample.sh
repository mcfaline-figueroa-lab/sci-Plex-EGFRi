#!/bin/bash
#
#
# Replace ACCOUNT with your account name before submitting.
#
#SBATCH --account=cgl      	 # Replace ACCOUNT with your group account name
#SBATCH --job-name=UMI_per_sample    # The job name
#SBATCH --time=0-10:00            # The time the job will take to run in D-HH:MM
#SBATCH --mem=20G
#SBATCH -o "batch-logs/8-UMI-per-sample-logs/%u-%x-%j.out"

FILTERED_READS_DIR=$1
RMDUP_SPLIT_BED_DIR=$2
BATCH=$3
OUTPUT_DIR=$4
BIN=$5
DATAMASH=$6

SAMTOOLS=$BIN/samtools

cat $BATCH | while read PCR_WELL; do
    awk '{
        split($4, arr, "|");
        if (!seen[arr[1]]) {
            seen[arr[1]] = 1;
            count[arr[2]]++;
        }
    } END {
        for (sample in count)
            print sample "\t" count[sample];
    }' $RMDUP_SPLIT_BED_DIR/$PCR_WELL.bed \
    | sort -k1,1 \
    >$OUTPUT_DIR/$PCR_WELL.UMI.count

    $SAMTOOLS view $FILTERED_READS_DIR/$PCR_WELL.bam \
    | cut -d '|' -f 2 \
    | $DATAMASH -g 1 count 1 \
    | sort -k1,1 -S 2G \
    | $DATAMASH -g 1 sum 2 \
    >$OUTPUT_DIR/$PCR_WELL.read.count

    echo "Processed $PCR_WELL"
done

