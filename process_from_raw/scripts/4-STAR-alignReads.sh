#!/bin/bash
#
#
# Replace ACCOUNT with your account name before submitting.
#
#SBATCH --account=cgl      	 # Replace ACCOUNT with your group account name
#SBATCH --job-name=STAR     # The job name
#SBATCH -c 16                     # The number of cores to request
#SBATCH --time=0-24:00            # The time the job will take to run in D-HH:MM
#SBATCH --mem-per-cpu=6G
#SBATCH -o "batch-logs/3-STAR-logs/%u-%x-%j.out"

module load STAR/2.7.9a

INPUT=$1
INDEX=$2
OUTPUT=$3


ls $INPUT | grep "[.]fq[.]gz$" | while read FILE; do
    SAMPLE=`basename "$FILE" _trimmed.fq.gz`

    STAR \
        --runThreadN 16 \
        --genomeDir $INDEX \
        --genomeLoad NoSharedMemory \
        --readFilesIn $INPUT/$FILE \
        --readFilesCommand zcat \
        --outFileNamePrefix $OUTPUT/$SAMPLE. \
        --outSAMtype BAM Unsorted \
        --outSAMstrandField intronMotif
done


