#!/bin/bash
# Step 5.1: Quality Control using fastp (v0.21.0)

# Define input and output
IN_R1="fastq/ERR13594515_1.fastq.gz"
IN_R2="fastq/ERR13594515_2.fastq.gz"
OUT_R1="clean_fastq/ERR13594515_1.clean.fq.gz"
OUT_R2="clean_fastq/ERR13594515_2.clean.fq.gz"

# Run fastp
fastp -i ${IN_R1} -I ${IN_R2} \
      -o ${OUT_R1} -O ${OUT_R2} \
      --thread 16 \
      -h clean_fastq/report.html \
      -j clean_fastq/report.json

echo "QC completed for sample ERR13594515"
