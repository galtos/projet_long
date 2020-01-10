#!/bin/bash

path_fasta="/home/sdv/m2bi/gollitrault/M2BI/projet_long/data/fasta"
path_pssm="/home/sdv/m2bi/gollitrault/M2BI/projet_long/data/pssm"

cd $path_fasta


for file in $(ls *_1.fasta);do
    struct_1=$file
    struct_2=$(echo $struct_1|sed "s:_1:_2:g")
    echo $struct_1
    echo $struct_2
    psiblast -query $struct_1 -db ../db_nr/uniref50.fasta -num_iterations 3 -out_ascii_pssm $struct_1.pssm
    psiblast -query $struct_2 -db ../db_nr/uniref50.fasta -num_iterations 3 -out_ascii_pssm $struct_2.pssm
done
