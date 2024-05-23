#!/bin/bash
# Final step for testing enrichment of gene sets in the annotated GWAS summ stats

# Directory containing .raw files (the annotated GWAS summ stats)
Raw_Files_Dir="/data1/meaneylab/eamon/MAGMA_annotated/20_20_window"

# Set annotation file path
Set_Annot_File="/data1/meaneylab/eamon/AM_microRNA_targets/MAGMA/MAGMA_input.txt"

# Output prefix
Output_Prefix_Base="AM_mir_targets"

# Loop through each .raw file in the directory
# This saves me manually running each analysis one by one and just runs them all one after another
for Gene_Results_File in "$Raw_Files_Dir"/*.raw; do
    # Extract the base name for the current .raw file to use in the output prefix
    File_Base_Name=$(basename "$Gene_Results_File" .raw)
    
    # Specify the output prefix including the base name of the current file
    Output_Prefix="${Output_Prefix_Base}_${File_Base_Name}"

    # Run enrichment
    magma --gene-results "$Gene_Results_File" --set-annot "$Set_Annot_File" col=1,2 --out "$Output_Prefix"
done
