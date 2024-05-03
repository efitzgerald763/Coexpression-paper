#!/bin/bash

mkdir temp_annot # make a temporary directory to host the intermediate files

Data_File="/data1/meaneylab/eamon/MAGMA/MAGMA_aux_files/1000_genomes_euro/g1000_eur"

Annot_File="/data1/meaneylab/eamon/MAGMA/MAGMA_aux_files/H-MAGMA_aux_files/HMAGMA_Protocol/Annotation_Files/Adultbrain.transcript.annot"

SNP_Pval_File="/data1/meaneylab/eamon/GWAS_sum_stats_MAGMA/Alzheimers_Bellenguez_2022.txt"

Output_Prefix="Alzheimers_Bellenguez_2022"


parallel magma --batch {} 10  --bfile $Data_File --gene-annot $Annot_File --gene-model snp-wise=mean --pval $SNP_Pval_File ncol=total_N --out temp_annot/$Output_Prefix ::: {1..10}


magma --merge temp_annot/$Output_Prefix --out temp_annot/$Output_Prefix


cp temp_annot/$Output_Prefix.genes.* .


rm -r temp_annot