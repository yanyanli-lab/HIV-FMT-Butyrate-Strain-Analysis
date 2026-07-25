#!/bin/bash
# Step 5.5: Calculating "Genomic Load" (Functional Potential) for each MAG
# This script counts the copy numbers of core butyrate-producing enzymes

ANNOTATION_FILE="eggnog_annotation_v5/final_results.emapper.annotations"
OUTPUT_CSV="all_bins_genomic_load.csv"

echo "MAG_ID,K01034,K00929,K01035" > $OUTPUT_CSV

# Based on history record #681 and #749
# Iterating through each MAG to count specific KEGG IDs
for faa in prodigal_results/*.faa; do
    mag=$(basename $faa .faa)
    
    # Counting the key terminal enzymes defined in Methods 5.4
    c1=$(grep -c "K01034" $ANNOTATION_FILE | grep "$mag" | wc -l)
    c2=$(grep -c "K00929" $ANNOTATION_FILE | grep "$mag" | wc -l)
    c3=$(grep -c "K01035" $ANNOTATION_FILE | grep "$mag" | wc -l)
    
    echo "$mag,$c1,$c2,$c3" >> $OUTPUT_CSV
done

echo "Genomic Load calculation completed. Results saved in $OUTPUT_CSV"
