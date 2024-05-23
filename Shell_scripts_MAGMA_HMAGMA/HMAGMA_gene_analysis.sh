#!/bin/bash

mkdir temp_annot # make a temporary directory to host the intermediate files

Data_File="/data1/meaneylab/eamon/MAGMA/MAGMA_aux_files/1000_genomes_euro/g1000_eur"

Annot_File="/data1/meaneylab/eamon/MAGMA/MAGMA_aux_files/H-MAGMA_aux_files/HMAGMA_Protocol/Annotation_Files/Adultbrain.transcript.annot"

SNP_Pval_File="/data1/meaneylab/eamon/GWAS_sum_stats_MAGMA/MDD_Als_et_al_exl_23nMe.txt"

Output_Prefix="MDD_Als_et_al_exl_23nMe"

# run magma in parallel, 4 threads in this case
parallel magma --batch {} 10  --bfile $Data_File --gene-annot $Annot_File --gene-model snp-wise=mean --pval $SNP_Pval_File ncol=total_N --out temp_annot/$Output_Prefix ::: {1..10}

# merge all intermediate files generated under the temp_annot files

# and send out for one single file set

magma --merge temp_annot/$Output_Prefix --out temp_annot/$Output_Prefix

# extract merged files for subsequent analysis

cp temp_annot/$Output_Prefix.genes.* .

# remove the temporary directory

rm -r temp_annot