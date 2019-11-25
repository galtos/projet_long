#!/usr/bin/env python
# coding: utf-8

# ##projet long prediction of prot prot interactions
# python test to download pdb file, extract prot sequence, set the resolution of the structure

from Bio.PDB import * # set ref
import numpy as np
import freesasa


# get pdb file names

def get_structure(name, path, str_structure_id):
    parser = PDBParser()
    return(parser.get_structure(name, path + str_structure_id))

# get the sequence
def get_sequence(structure):
    ppb = CaPPBuilder()
    for pp in ppb.build_peptides(structure):
        polypeptide = pp.get_sequence()
    return(list(polypeptide))

def get_surface_residue(structure, path_structure, rsa_threshold, acc_array = "Sander"):
    list_surface_residue = []
    dssp = DSSP(structure, path_structure, acc_array = acc_array)
    len_chain = []
    sum_len_chain = 0
    for chain in structure:
        len_chain.append(len(chain))
        
    list_dssp = list(dssp)
    print(len(len_chain), "__")
    for i in range(len(len_chain)):
        vec_surface_residue = []
        for j in range(sum_len_chain, sum_len_chain + len_chain[i]):
            print(j, sum_len_chain, len_chain[i])
            if list_dssp[j][3] >= rsa_threshold:
                vec_surface_residue.append([list_dssp[j][0], list_dssp[j][3]])
                
        sum_len_chain =  sum_len_chain + len_chain[i]
        list_surface_residue.append(vec_surface_residue)
        
    return list_surface_residue


def residue_distance(residue_1, residue_2):
    try:
        distance = residue_1['CA'] - residue_2['CA']
    except KeyError:
        print("no CA atoms")
    return distance

def get_interface_residue(structure_1, structure_2, distance_threshold = 8):
    residues_1 = [r for r in structure_1.get_residues()]
    residues_2 = [r for r in structure_2.get_residues()]
    list_interface_residue = []
    vector_interface_1 = np.zeros(len(residues_1))
    vector_interface_2 = np.zeros(len(residues_2))
    print(len(residues_1), len(residues_2))
    for i in range(len(residues_1)):
        for j in range(len(residues_2)):
            distance = residue_distance(residues_1[i], residues_2[j])
            if distance <= distance_threshold:
                list_interface_residue.append([i, j, distance])
                vector_interface_1[i] = 1
                vector_interface_2[j] = 1
    return(list_interface_residue, vector_interface_1, vector_interface_2)

#### MAIN ####
if __name__ == "__main__":
    path_bound = "../data/data_struct3d_bound"
    path_list_bound_pdb_file = path_bound + "/full_structures.0.90.txt"
    
    with open(path_list_bound_pdb_file) as file:
        list_bound_pdb_file = file.readlines()

    print(list_unbound_pdb_file[0], list_bound_pdb_file[0])

    structure_1 = []
    structure_2 = []
    structure_12 = []
    structure_1.append(get_structure('test_bound_1', path_bound + "/templates/" + list_bound_pdb_file[0][0:-1], '_1.pdb')
    structure_2.append(get_structure('test_bound_2', path_bound + "/templates/" + list_bound_pdb_file[0][0:-1], '_2.pdb'))
    structure_12.append(get_structure('test_bound_12', path_bound + "/templates_bind/" + list_bound_pdb_file[0][0:-1],'.pdb'))
    
    test, vec1, vec2 = get_interface_residue(structure_1[0], structure_2[0])
    
    list_surface_residue = get_surface_residue(structure_1[0][0], path_bound + "/templates/" + list_bound_pdb_file[0][0:-1] + '_1.pdb', 0.1)

    list_surface_residue_12 = get_surface_residue(structure_12[0][0], path_bound + "/templates/" + list_bound_pdb_file[0][0:-1] + '_1.pdb', 0.1)

    list_surface_residue_12[1]

    res = [value for value in list_surface_residue if value not in list_surface_residue_12[1]]

    for i in range(len(list_surface_residue)):
        if list_surface_residue[0][0][1] != list_surface_residue_12[0][i][1]:
            print(i)

    list_surface_residue_12[0][0][1]

    list_surface_residue[0][0][1]
    

