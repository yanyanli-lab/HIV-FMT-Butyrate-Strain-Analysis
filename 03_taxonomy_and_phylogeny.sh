#!/bin/bash
# Step 5.3: Taxonomic identification and Phylogeny construction

# 1. Taxonomic classification using GTDB-Tk (v2.1.1 or latest)
# This command runs the full workflow: identify, align, and classify
gtdbtk classify_wf \
    --genome_dir ./bins \
    --out_dir ./final_gtdbtk_results \
    --extension fasta \
    --cpus 40

# 2. Extracting the Multiple Sequence Alignment (MSA)
# GTDB-Tk generates a filtered MSA of 120 marker genes
gunzip -c ./final_gtdbtk_results/align/gtdbtk.bac120.user_msa.fasta.gz > gtdbtk.bac120.user_msa.fasta

# 3. De novo Phylogenetic Tree construction using FastTree
# Based on your history record #820 and #824
# This tree file was used for iTOL visualization
FastTree -nt gtdbtk.bac120.user_msa.fasta > My_279_MAGs_Clean.tree

echo "Taxonomy and Tree construction completed."
