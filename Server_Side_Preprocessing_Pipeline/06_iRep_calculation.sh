#!/bin/bash
# Step 5.6: In situ replication rate (iRep) calculation
# This script estimates the growth rates of core strains during the W08 burst.

# 1. Mapping clean reads back to the MAGs index
# Based on history record #877
for r1 in ./clean_fastq/*_1.fastq.gz; do
    prefix=$(basename $r1 | sed 's/_1.fastq.gz//')
    r2="./clean_fastq/${prefix}_2.fastq.gz"
    echo "Processing sample: ${prefix} ..."
    
    # Using bowtie2 for precise mapping
    bowtie2 -x MAGs_index -1 $r1 -2 $r2 -p 16 --reorder --no-unal 2> mapping_log.txt | \
    samtools view -bS - > ./iRep_mapping/${prefix}.bam
done

# 2. Calculating iRep for the "Top 15" core strains
# Based on history record #893 and #914
# We focused on the Top 15 MAGs to validate their cellular division activity.
iRep -f ./top15_test/*.fa \
     -s ./iRep_mapping/*.bam \
     -o Top15_W08_Activity_Results \
     -t 16 \
     -c 3 -r2 0.7

echo "iRep calculation completed. Results saved in Top15_W08_Activity_Results."
