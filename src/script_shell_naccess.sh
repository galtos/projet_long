#!/bin/bash

path_templates="/home/sdv/m2bi/gollitrault/M2BI/projet_long/data/data_struct3d_bound/templates"
path_templates_bind="/home/sdv/m2bi/gollitrault/M2BI/projet_long/data/data_struct3d_bound/templates_bind"

cd $path_templates


for file in $(ls *_1.pdb);do
    struct_1=$file
    struct_2=$(echo $struct_1|sed "s:_1:_2:g")
    echo $struct_1
    echo $struct_2
    ./../../../software/Naccess/Naccess-correct-1/naccess $struct_1
    ./../../../software/Naccess/Naccess-correct-1/naccess $struct_2
done

cd $path_templates_bind

for file in $(ls *.pdb);do
    struct=$file
    ./../../../software/Naccess/Naccess-correct-1/naccess $struct
done


