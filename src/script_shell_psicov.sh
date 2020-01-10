#!/bin/bash

path_fasta="/home/sdv/m2bi/gollitrault/M2BI/projet_long/data/fasta"

cd $path_fasta

for file in $(ls *_1.fasta);do
    struct_1=$file
    struct_2=$(echo $struct_1|sed "s:_1:_2:g")
    echo $struct_1
    echo $struct_2
    hhblits -i $struct_1 -maxfilt 100000 -realign_max 100000 -d ../uniclust30_2018_08/uniclust30_2018_08 -all -B 100000 -Z 100000 -n 3 -e 0.001 -oa3m $struct_1.a3m
    
    hhblits -i $struct_2 -maxfilt 100000 -realign_max 100000 -d ../uniclust30_2018_08/uniclust30_2018_08 -all -B 100000 -Z 100000 -n 3 -e 0.001 -oa3m $struct_2.a3m
    
    egrep -v "^>" $struct_1.a3m | sed 's/[a-z]//g' | sort -u > $struct_1.aln
    egrep -v "^>" $struct_2.a3m | sed 's/[a-z]//g' | sort -u > $struct_2.aln
    
    ./../software/psicov/src/psicov -p -r 0.001 $struct_1.aln > $struct_1.psicov
    ./../software/psicov/src/psicov -p -r 0.001 $struct_2.aln > $struct_2.psicov
done
