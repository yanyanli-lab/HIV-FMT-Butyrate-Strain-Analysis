#!/bin/bash
# Step 5.4: Functional annotation using eggNOG-mapper (v2.1.9)

# Running emapper on predicted protein sequences (from Prokka/Prodigal)
# Based on history record #355
emapper.py -i all_279_proteins_FIXED.faa \
           --output final_results \
           --output_dir ./eggnog_annotation_v5/ \
           -m diamond \
           --cpu 16 \
           --override

# The resulting '.emapper.annotations' file contains the KEGG Orthology (KO) assignments
# used for subsequent butyrate pathway analysis.
