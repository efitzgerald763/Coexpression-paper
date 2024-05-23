#!/bin/bash

# First step, you only have to do this once oer threshold and then you can use the output file from this to annotate your GWAS

# SNP file from the MAGMA website
SNP_Loc_File="/data1/meaneylab/eamon/MAGMA/MAGMA_aux_files/1000_genomes_euro/g1000_eur.bim"

# Gene file from the MAGMA website
Gene_Loc_File="/data1/meaneylab/eamon/MAGMA/MAGMA_aux_files/SNPs/NCBI38.gene.loc"

# Name the output file
Output_Prefix="MAGMA_20k_20k"

# This is the actual MAGMA function, which creates associates SNPs within a certain distance of a gene (here 20Kb upstream and 20Kb downstream)
# to that gene, which is then used to map GWAS SNPs to the gene. Setting the distance threshold is a bit arbitrary with no real right or
# wrong value. 20Kb usually works well for me and is commonly used in the literature, but you can experiment with it.
magma --annotate window=20,20 --snp-loc $SNP_Loc_File  --gene-loc $Gene_Loc_File  --out $Output_Prefix