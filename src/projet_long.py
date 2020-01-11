#!/usr/bin/env python
# coding: utf-8

# ##projet long prediction of prot prot interactions
# python test to download pdb file, extract prot sequence, set the resolution of the structure

from Bio.PDB import * # set ref
#import DSSP
import numpy as np
import aaindex
import os
import re
from subprocess import call
import sys
import unittest
import operator
from DeepNN_model import DeepNN_model_build
from keras.utils import to_categorical
#pssm
from Bio import AlignIO
#mutual information
#from MI import Statistics

									
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
"""
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
"""

def residue_distance(residue_1, residue_2):
    try:
        distance = residue_1['CA'] - residue_2['CA']
    except KeyError:
        print("no CA atoms")
    return distance
"""
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
"""

def get_neighbor_residues(structure, residue_middle, k_threshold=9):
    residues = [r for r in structure.get_residues()]
    residue_distance = []
    for i in range(len(residues)):
        distance = residue_middle["CA"] - residues[i]["CA"]
        residue_distance.append((i,distance))
    return(sorted(residue_distance, key=operator.itemgetter(1))[0:k_threshold])

def get_rsa_relative(path_file_rsa):
    list_relative_rsa = []
    for line in open(path_file_rsa):
        list = line.split()
        id = list[0]
        if id == 'RES':
            relative_rsa = list[5]
            list_relative_rsa.append(float(relative_rsa))
    return(list_relative_rsa)
def get_asa_value(path_file_asa):
    list_asa_value = []
    for line in open(path_file_asa):
        list = line.split()
        id = list[0]
        if id == 'RES':
            asa_value = list[4]
            list_asa_value.append(float(asa_value))
    return(list_asa_value)
    
def get_rsa_relative_bind(path_file_rsa):
    list_relative_rsa = [[]]
    struct_flag = 2
    for line in open(path_file_rsa):
        list = line.split()
        id = list[0]
        
        if id == 'RES':
            res_number = float(list[3])
            if struct_flag == 2:
                struct_flag = 0
            elif struct_flag == 0 and res_number != (res_number_previous+1):
                struct_flag = 1
                list_relative_rsa.append([])
            relative_rsa = list[5]
            list_relative_rsa[struct_flag].append(float(relative_rsa))
            res_number_previous = float(list[3])
    return((list_relative_rsa[0],list_relative_rsa[1]))
    
def get_surface_residue(list_relative_rsa, threshold_rsa = 0.25):
    list_surface_residue = []
    for i in range(len(list_relative_rsa)):
        if list_relative_rsa[i] >= threshold_rsa:
            list_surface_residue.append(i)
    return(list_surface_residue)
    
def get_interface_residue(list_relative_rsa, list_relative_rsa_bind):
    list_interface_residue = []
    for i in range(len(list_relative_rsa)):
        if list_relative_rsa[i] != list_relative_rsa_bind[i]:
            list_interface_residue.append(i)
    return(list_interface_residue)

def count_intersect_surface_interface_residue(list_surface_residue, list_interface_residue):
    k = 0
    for i in range(len(list_surface_residue)):
        for j in range(len(list_interface_residue)):
            if list_surface_residue[i] == list_interface_residue[j]:
                k = k+1
    print(k)
def get_pssm_value(path_pssm):
    list_ipp_value = []
    list_rwgrmtp_value = []
    for line in open(path_pssm):
        line = line.rstrip()
        line = re.sub(' +',' ',line)
        line = line.split(' ')
        if len(line) >= 44:
            list_ipp_value.append(float(line[-2]))
            list_rwgrmtp_value.append( float(line[-1]) if line[-1] != "inf" else "inf" )
    return((list_ipp_value,list_rwgrmtp_value))

def get_SS(structure, path_structure):
    list_SS = []
    dssp = DSSP(structure, path_structure)
    list_dssp = list(dssp)
    print(list_dssp)
"""
def set_propka():
    pdbs = ['1FTJ-Chain-A',
                '1HPX',
                '4DFR']

        test_dir = os.path.dirname(__file__)
        base_dir = os.path.dirname(test_dir)

        executable = os.path.join(base_dir, "scripts", "propka31.py")

        env = { "PYTHONPATH" : base_dir }

        for pdb in pdbs:
            input_filename = os.path.join(test_dir, "pdb", pdb + ".pdb")
            output_filename = os.path.join(test_dir, pdb + ".out")

            output_file = open(output_filename, "w")
            call([sys.executable, executable, input_filename],
                    stdout=output_file, env=env)
            output_file.close()
"""
def get_aaindex_features(AA):
    aaindex.init(path='../data/aaindex')
    Hydrophobicity = aaindex.get('ARGP820101')
    Hydrophilicity = aaindex.get('HOPT810101')
    Polarity = aaindex.get('GRAR740102')
    Polarizability = aaindex.get('CHAM820101')
    Propensity = aaindex.get('WERD780101')
    Average_accessible_surface_area = aaindex.get('JANJ780101')
    Radius_of_gyration_of_side_chain = aaindex.get('LEVM760105')
    Side_chain_volume = aaindex.get('KRIW790103')
    Charge = aaindex.get('KLEP840101')
    Number_of_hydrogen_bond_donors = aaindex.get('FAUJ880109')
    Molecular_weight = aaindex.get('FASG760101')
    Electron_ion_interaction_potential = aaindex.get('VELV850101')
#code to get features
    aaindex_feature = [Hydrophobicity.get(AA),\
                       Hydrophilicity.get(AA),\
                       Polarity.get(AA),\
                       Polarizability.get(AA),\
                       Propensity.get(AA),\
                       Average_accessible_surface_area.get(AA),\
                       Radius_of_gyration_of_side_chain.get(AA),\
                       Side_chain_volume.get(AA),\
                       Charge.get(AA),\
                       Number_of_hydrogen_bond_donors.get(AA),\
                       Molecular_weight.get(AA),\
                       Electron_ion_interaction_potential.get(AA)]
    return(aaindex_feature)
    
def get_pseudo_hydrophobicity(AA, Hydrophobicity, Charge):
    if Charge.get(AA) < 0:
        return(Hydrophobicity.get(AA)*Charge.get(AA))
    return(Hydrophobicity.get(AA))

def get_max_psicov(path_psicov, size_seq):
    list_max_coev = np.zeros((size_seq, size_seq))
    for line in open(path_psicov):
        line = line.rstrip()
        line = line.split(' ')
        list_max_coev[int(line[0])-1][int(line[1])-1] = float(line[4])
    return(np.amax(list_max_coev, axis=1))

def get_vector(structure, path_file_pssm, path_aaindex, path_file_rsa, path_file_asa):
    aaindex.init(path_aaindex)

    Hydrophobicity = aaindex.get('ARGP820101')
    Hydrophilicity = aaindex.get('HOPT810101')
    Polarity = aaindex.get('GRAR740102')
    Polarizability = aaindex.get('CHAM820101')
    Propensity = aaindex.get('WERD780101')
    Average_accessible_surface_area = aaindex.get('JANJ780101')
    Radius_of_gyration_of_side_chain = aaindex.get('LEVM760105')
    Side_chain_volume = aaindex.get('KRIW790103')
    Charge = aaindex.get('KLEP840101')
    Number_of_hydrogen_bond_donors = aaindex.get('FAUJ880109')
    Molecular_weight = aaindex.get('FASG760101')
    Electron_ion_interaction_potential = aaindex.get('VELV850101')
    
    list_vector = []
    """
    residues = structure_1[0].get_residues()
    for residue in residues:
        residue.get_resname()
        
    ppb=PPBuilder()
    """
    #asa
    list_asa_value = get_asa_value(path_file_asa)
    #rsa
    list_rsa_value = get_rsa_relative(path_file_rsa)
    #pssm

    list_ipp_value, list_rwgrmtp_value = get_pssm_value(path_file_pssm)
    #QIPI
    QIPI = {
'H':1.147, 'R':1.346, 'K':0.784, 'A':0.841, 'V':0.994, 'I':1.084, 'L':1.144, \
'M':1.451, 'P':1.109, 'F':1.334, 'W':1.284, 'Y':1.368, 'G':0.823, 'C':1.172, \
'S':0.873, 'T':0.966, 'N':0.958, 'Q':0.909, 'D':0.830, 'E':0.805}
    ppb = PPBuilder()
    for pp in ppb.build_peptides(structure):
        sequence = list(pp.get_sequence())
    for i in range(len(sequence)):
        list_vector.append([\
                           #asa
                           list_asa_value[i],\
                           #rsa
                           list_rsa_value[i],\
                           #pssm
                           list_ipp_value[i],\
                           list_rwgrmtp_value[i],\
                           #QIPI
                           QIPI[sequence[i]],\
                           #pseudo hydrophobicity
                           get_pseudo_hydrophobicity(sequence[i],\
                                                     Hydrophobicity,\
                                                     Charge),\
                           #aaindex
                           Hydrophobicity.get(sequence[i]),\
                           Hydrophilicity.get(sequence[i]),\
                           Polarity.get(sequence[i]),\
                           Polarizability.get(sequence[i]),\
                           Propensity.get(sequence[i]),\
                           Average_accessible_surface_area.get(sequence[i]),\
                           Radius_of_gyration_of_side_chain.get(sequence[i]),\
                           Side_chain_volume.get(sequence[i]),\
                           Charge.get(sequence[i]),\
                           Number_of_hydrogen_bond_donors.get(sequence[i]),\
                           Molecular_weight.get(sequence[i]),\
                           Electron_ion_interaction_potential.get(sequence[i])\
                           ])
                           #             
    return(list_vector)
def get_vector_neighbors(structure, list_vector):
    list_vector_neighbors = []
    residues = structure.get_residues()
    for residue in residues:
        list_neighbors = get_neighbor_residues(structure, residue)
        vector_neighbors = []
        for neigbor in range(len(list_neighbors)):
            vector_neighbors = vector_neighbors + list_vector[list_neighbors[neigbor][0]]
        list_vector_neighbors.append(vector_neighbors)
    return(list_vector_neighbors)

def get_array_vector(structure_1, structure_2):
    """
    bind vector structure1 and structure 2
    """
    for pp in ppb.build_peptides(structure_1[0]):
        sequence_1 = list(pp.get_sequence())
    for pp in ppb.build_peptides(structure_2[0]):
        sequence_2 = list(pp.get_sequence())
        
    list_vector_1 = get_vector(structure_1[0], path_file_pssm, "", path_file_rsa, path_file_rsa)
    list_vector_2 = get_vector(structure_2[0], path_file_pssm, "", path_file_rsa, path_file_rsa)
    list_vector_neighbors_1 = get_vector_neighbors(structure_1[0], list_vector_1)
    list_vector_neighbors_2 = get_vector_neighbors(structure_2[0], list_vector_2)
    array_vector = []
    for i in range(len(sequence_1)):
        for j in range(len(sequence_2)):
            array_vector.append(list_vector_neighbors_1[i]+list_vector_neighbors_2[j])
    return(array_vector)

def get_X_Y(path_bound,\
            path_list_bound_pdb_file,\
            path_file_rsa,\
            path_file_rsa_bind,\
            path_file_asa,\
            path_file_asa_bind,\
            path_file_pssm,\
            path_aaindex):
    """
    get X and Y vector for training
    """
    with open(path_list_bound_pdb_file) as file:
        list_bound_pdb_file = file.readlines()
    ## get structures
    parser = PDBParser()
    Y = []
    X = []
    XY = []
    for i in range(10):
        structure_1 = parser.get_structure('test_bound_1', path_bound + "/templates/" + list_bound_pdb_file[i][0:-1]+ '_1.pdb')
        structure_2 = parser.get_structure('test_bound_2', path_bound + "/templates/" + list_bound_pdb_file[i][0:-1]+ '_2.pdb')
        structure_12 = parser.get_structure('test_bound_12', path_bound + "/templates_bind/" + list_bound_pdb_file[i][0:-1]+'.pdb')
        
        residues_1 = [r for r in structure_1.get_residues()]
        residues_2 = [r for r in structure_1.get_residues()]
        if len(residues_1) == len(get_sequence(structure_1)) and len(residues_2) == len(get_sequence(structure_2)):
            
            list_value_rsa_1 = get_rsa_relative(path_file_rsa + list_bound_pdb_file[i][0:-1]+ '_1.rsa')
            list_value_rsa_2 = get_rsa_relative(path_file_rsa + list_bound_pdb_file[i][0:-1]+ '_2.rsa')
            
            list_relative_rsa_bind_1, list_relative_rsa_bind_2 = get_rsa_relative_bind(path_file_rsa_bind + list_bound_pdb_file[i][0:-1]+ '.rsa')
            
            
            list_surface_residue_1 = get_surface_residue(list_value_rsa_1, threshold_rsa = 25)
            list_surface_residue_2 = get_surface_residue(list_value_rsa_2, threshold_rsa = 25)
            
            list_interface_residue_1 = get_interface_residue(list_value_rsa_1, list_relative_rsa_bind_1)
            list_interface_residue_2 = get_interface_residue(list_value_rsa_2, list_relative_rsa_bind_2)

            list_vector_1 = get_vector(structure_1,\
                                       path_file_pssm + list_bound_pdb_file[i][0:-1] + '_1.fasta.pssm',\
                                       path_aaindex,\
                                       path_file_rsa + list_bound_pdb_file[i][0:-1] + '_1.rsa',\
                                       path_file_rsa + list_bound_pdb_file[i][0:-1] + '_1.rsa')
            list_vector_2 = get_vector(structure_2,\
                                       path_file_pssm + list_bound_pdb_file[i][0:-1] + '_2.fasta.pssm',\
                                       path_aaindex,\
                                       path_file_rsa + list_bound_pdb_file[i][0:-1] + '_2.rsa',\
                                       path_file_rsa + list_bound_pdb_file[i][0:-1] + '_2.rsa')
            list_vector_neighbors_1 = get_vector_neighbors(structure_1, list_vector_1)
            list_vector_neighbors_2 = get_vector_neighbors(structure_2, list_vector_2)

            list_vector_neighbors_1_surface = []
            list_vector_neighbors_2_surface = []
            
            for i in list_surface_residue_1:
                list_vector_neighbors_1_surface.append(list_vector_neighbors_1[i])
                if i in list_interface_residue_1:
                    Y.append(1)
                else:
                    Y.append(0)

            for i in list_surface_residue_2:
                list_vector_neighbors_2_surface.append(list_vector_neighbors_2[i])
                if i in list_interface_residue_2:
                    Y.append(1)
                else:
                    Y.append(0)
            X = X + list_vector_neighbors_1_surface + list_vector_neighbors_2_surface
            XY.append([len(X),len(Y)])
        else:
            print("BAD number of residues do not correspond dont know why")
    Y = to_categorical(Y, num_classes=2)
    X = np.asarray(X)
    return(X,Y)

def train_DeepNN_model(X, Y):
    model = DeepNN_model_build.build()
    model.compile(loss = "binary_crossentropy",optimizer="adam", metrics=['accuracy'])
    model.fit(X,Y,epochs=1000, batch_size = 50)
    print(model.evaluate(X, Y))
    
def make_fasta(list_bound_pdb_file, path_bound, path_fasta = "../data/fasta/"):
    parser = PDBParser()
    for i in range(len(list_bound_pdb_file)):
        structure_1 = parser.get_structure(list_bound_pdb_file[i][0:-1]+ '_1', path_bound + "/templates/" + list_bound_pdb_file[i][0:-1]+ '_1.pdb')
        structure_2 = parser.get_structure(list_bound_pdb_file[i][0:-1]+ '_2', path_bound + "/templates/" + list_bound_pdb_file[i][0:-1]+ '_2.pdb')
        
        file_fasta_1 = open(path_fasta + list_bound_pdb_file[i][0:-1]+'_1.fasta', 'w')
        file_fasta_1.write('>'+list_bound_pdb_file[i][0:-1]+ '_1'+"\n")
        file_fasta_1.write(''.join(list(get_sequence(structure_1))))
        file_fasta_1.close()
        file_fasta_2 = open(path_fasta + list_bound_pdb_file[i][0:-1]+'_2.fasta', 'w')
        file_fasta_2.write('>'+list_bound_pdb_file[i][0:-1]+ '_2'+"\n")
        file_fasta_2.write(''.join(list(get_sequence(structure_2))))
        file_fasta_2.close()
def make_pssm(list_bound_pdb_file, path_bound, path_pssm = "../data/pssmm/", path_aln = "../data/example.aln"):
    align = AlignIO.read("../data/example.a3m", "fasta-m10")
    print(align)

#### MAIN ####
if __name__ == "__main__":

    #PATHS
    path_bound = "../data/data_struct3d_bound"
    path_list_bound_pdb_file = path_bound + "/full_structures.0.90.txt"

    #get file names
    with open(path_list_bound_pdb_file) as file:
        list_bound_pdb_file = file.readlines()
    #RSA
    path_file_rsa = path_bound + "/templates/"
    path_file_rsa_bind = path_bound + "/templates_bind/"
    #ASA
    path_file_asa = path_bound + "/templates/"
    path_file_asa_bind = path_bound + "/templates_bind/"
    #PSSM
    path_file_pssm = "../data/pssm/"
    #aaindex
    path_aaindex='../data/aaindex'
    
    X, Y = get_X_Y(path_bound,\
                   path_list_bound_pdb_file,\
                   path_file_rsa,\
                   path_file_rsa_bind,\
                   path_file_asa,\
                   path_file_asa_bind,\
                   path_file_pssm,\
                   path_aaindex)
    train_DeepNN_model(X, Y)








