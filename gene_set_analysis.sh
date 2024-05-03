#!/bin/bash

Gene_Results_File="GWAS.genes.raw"

Set_Annot_File="Genesets_for_enrichment.txt"

Output_Prefix="GWAS_enrichment"

magma --gene-results $Gene_Results_File  --set-annot $Set_Annot_File col=1,2  --out $Output_Prefix