#!/bin/bash

Gene_Results_File="/data1/meaneylab/eamon/H_MAGMA_adultbrain/Subjective_well_being_2016.genes.raw"

Set_Annot_File="/data1/meaneylab/eamon/MAGMA/All_networks_MAGMA.txt"

Output_Prefix="Subjective_well_being_2016"

magma --gene-results $Gene_Results_File  --set-annot $Set_Annot_File col=1,2  --out $Output_Prefix