#!/bin/bash
#
#
# Replace ACCOUNT with your account name before submitting.
#
#SBATCH --account=cgl      	 # Replace ACCOUNT with your group account name
#SBATCH --job-name=count_rRNA    # The job name
#SBATCH --time=0-4:00            # The time the job will take to run in D-HH:MM
#SBATCH --mem=10G
#SBATCH -o "batch-logs/5-count-rRNA-logs/%u-%x-%j.out"

INPUT_DIR=$1
FILE_LIST=$2
RRNA_BED=$3
OUTPUT_DIR=$4
BIN=$5

BEDTOOLS=$BIN/bedtools

cat $FILE_LIST | while read FILE; do
    PCR_WELL=`basename "$FILE" .Aligned.out.bam`

    $BEDTOOLS intersect -a $INPUT_DIR/$FILE -b $RRNA_BED -c -nonamecheck -bed \
    | awk '{
        split($4, arr, "|");
        if (!seen[arr[1]]) {
            if ($NF > 0)
                rrna_count[arr[2]]++;
            total_count[arr[2]]++;
            seen[arr[1]] = 1;
        }
    } END {
        for (sample in total_count)
            printf "%s\t%d\t%d\n",
                sample, rrna_count[sample], total_count[sample];
    }' \
    >$OUTPUT_DIR/$PCR_WELL

    echo "Processed $FILE"
done

