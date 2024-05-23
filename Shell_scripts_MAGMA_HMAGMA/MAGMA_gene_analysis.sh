#!/bin/bash

# This is the second step and it annotates the GWAS summary statisitcs to the gene level, giving each gene a P-value
# This is pretty slow so I usually run it in parrallel

mkdir temp_annot # Makes a temporary directory to host the intermediate files

# File from MAGMA website
Data_File="/data1/meaneylab/eamon/MAGMA/MAGMA_aux_files/1000_genomes_euro/g1000_eur"

# File generated in previous step
Annot_File="/data1/meaneylab/eamon/MAGMA_annotated/20_20_window/MAGMA_20k_20k.genes.annot"

# GWAS summ stats, MAGMA expects a specific format of summ stat so sometimes you have to do some formatting
SNP_Pval_File="/data1/meaneylab/eamon/GWAS_sum_stats_MAGMA/MDD_Als_et_al_exl_23nMe.txt"

Output_Prefix="MDD_Als_et_al_exl_23nMe"

# Run the annotation in parallel, 10 threads in this case
parallel magma --batch {} 10  --bfile $Data_File --gene-annot $Annot_File --gene-model snp-wise=mean --pval $SNP_Pval_File ncol=total_N --out temp_annot/$Output_Prefix ::: {1..10}

# Merge all intermediate files generated under the temp_annot files and output a single file 

magma --merge temp_annot/$Output_Prefix --out temp_annot/$Output_Prefix

# Extract merged files for subsequent analysis

cp temp_annot/$Output_Prefix.genes.* .

# Remove the temporary directory

rm -r temp_annot