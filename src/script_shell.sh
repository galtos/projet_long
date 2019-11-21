#!/bin/bash
#mkdir ../data/templates_bind

path_templates="../data/data_struct3d_bound/templates"
path_templates_bind="../data/data_struct3d_bound/templates_bind"
for file in $(ls $path_templates/*_1.pdb);do
    struct_1=$file
    struct_2=$(echo $struct_1|sed "s:_1:_2:g")
    struct=$(echo $struct_1|sed "s:_1::g"|sed "s:$path_templates:$path_templates_bind:g")
    echo $struct_1
    echo $struct_2
    echo $struct
    cat $struct_1 $struct_2 > $struct 
done
