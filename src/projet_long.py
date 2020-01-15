#!/usr/bin/env python
# coding: utf-8

# projet long prediction of prot prot interactions

from Bio.PDB import * # set ref
from Bio.PDB.DSSP import DSSP
import argparse
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
from keras.utils import  np_utils
from sklearn.metrics import confusion_matrix
from sklearn.metrics import roc_curve
from sklearn.metrics import auc
#pssm
from Bio import AlignIO
#mutual information
#from MI import Statistics
#numpy to npy
from numpy import asarray
from numpy import save
from numpy import load
#pyplot
import matplotlib.pyplot as plt
#Voxelization
from moleculekit.molecule import Molecule
from moleculekit.tools.voxeldescriptors import getVoxelDescriptors,\
                                               viewVoxelFeatures
from moleculekit.tools.atomtyper import prepareProteinForAtomtyping
from moleculekit.smallmol.smallmol import SmallMol
from moleculekit.home import home
import os
#
import torch
import torch.nn as nn
import torch.nn.functional as F
#
from sklearn.utils import class_weight

									
# get pdb file names

def get_structure(name, path, str_structure_id):
    """
    take a path to a pdb file and return bioparser architecture
    """
    parser = PDBParser()
    return(parser.get_structure(name, path + str_structure_id))

# get the sequence
def get_sequence(structure):
    """
    take a bioparser structure and return its sequence
    """
    ppb = CaPPBuilder()
    for pp in ppb.build_peptides(structure):
        polypeptide = pp.get_sequence()
    return(list(polypeptide))


def residue_distance(residue_1, residue_2):
    """
    Compute the distance between the CA of 2 residues
    """
    try:
        distance = residue_1['CA'] - residue_2['CA']
    except KeyError:
        print("no CA atoms")
    return distance

def get_neighbor_residues(structure, residue_middle, k_threshold=9):
    residues = [r for r in structure.get_residues()]
    residue_distance = []
    for i in range(len(residues)):
        distance = residue_middle["CA"] - residues[i]["CA"]
        residue_distance.append((i,distance))
    return(sorted(residue_distance, key=operator.itemgetter(1))[0:k_threshold])

def get_rsa_relative(path_file_rsa):
    """
    get a list of the RSA values for one structure given the Naccess .rsa file
    """
    list_relative_rsa = []
    for line in open(path_file_rsa):
        list = line.split()
        id = list[0]
        if id == 'RES':
            relative_rsa = list[5]
            list_relative_rsa.append(float(relative_rsa))
    return(list_relative_rsa)
def get_asa_value(path_file_asa):
    """
    get a list of the ASA values for one srtuctures given the Naccess .rsa file
    """
    list_asa_value = []
    for line in open(path_file_asa):
        list = line.split()
        id = list[0]
        if id == 'RES':
            asa_value = list[4]
            list_asa_value.append(float(asa_value))
    return(list_asa_value)
    
def get_rsa_relative_bind(path_file_rsa):
    """
    get a list of the RSA values for two srtuctures in their bind state
     given the Naccess .rsa file. return 2 lists of their RSA values
    """
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
    """
    given a list of RSA values get the surface residues according a threshold
    """
    list_surface_residue = []
    for i in range(len(list_relative_rsa)):
        if list_relative_rsa[i] >= threshold_rsa:
            list_surface_residue.append(i)
    return(list_surface_residue)
    
def get_interface_residue(list_relative_rsa, list_relative_rsa_bind):
    """
    given a list of the RSA values of the structure alone and bind return
    a list of interface residues
    """
    list_interface_residue = []
    for i in range(len(list_relative_rsa)):
        if list_relative_rsa[i] != list_relative_rsa_bind[i]:
            list_interface_residue.append(i)
    return(list_interface_residue)

def count_intersect_surface_interface_residue(list_surface_residue, list_interface_residue):
    """
    count the number of surface residues in the interface
    """
    k = 0
    for i in range(len(list_surface_residue)):
        for j in range(len(list_interface_residue)):
            if list_surface_residue[i] == list_interface_residue[j]:
                k = k+1
    print(k)
def get_pssm_value(path_pssm):
    """
    parse pssm file to get a list of the ipp values and rwgrmtp values
    """
    list_ipp_value = []
    list_rwgrmtp_value = []
    for line in open(path_pssm):
        line = line.rstrip()
        line = re.sub(' +',' ',line)
        line = line.split(' ')
        if len(line) >= 44:
            list_ipp_value.append(float(line[-2]))
            list_rwgrmtp_value.append(\
            float(line[-1]) if line[-1] != "inf" else "inf" )
    return((list_ipp_value,list_rwgrmtp_value))

def get_SS(structure, path_structure):
    """
    given a structure
    return a list of 3 features coresspondig to the secondary structure 
    for each residue
    Helix : (0,0,1)
    Sheet : (0,1,0)
    Coil : (1,0,0)
    """
    list_SS = []
    dssp = DSSP(structure[0], path_structure)
    list_dssp = list(dssp)
    list_dssp_features = []

    for i in range(len(list_dssp)):
        value_dssp = list(list_dssp[i])
        if value_dssp[2] == "G" or\
           value_dssp[2] == "H" or\
           value_dssp[2] == "I":
            list_dssp_features.append([1,0,0])
        elif value_dssp[2] or value_dssp[2] == "B" or value_dssp[2] =="E":
            list_dssp_features.append([0,1,0])
        elif value_dssp[2] or value_dssp[2] == "T" or\
       value_dssp[2] =="S" or value_dssp[2] =="-":
            list_dssp_features.append([0,1,0])
    return(list_dssp_features)

def get_aaindex_features(AA):
    """
    given an amino acid return its physicochemical features
    """
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
    """
    for an amino acid return its pseudo hydrophobicity
    """
    if Charge.get(AA) < 0:
        return(Hydrophobicity.get(AA)*Charge.get(AA))
    return(Hydrophobicity.get(AA))

def get_max_psicov(path_psicov, size_seq):
    """
    get the max value of psicov data
    """
    list_max_coev = np.zeros((size_seq, size_seq))
    for line in open(path_psicov):
        line = line.rstrip()
        line = line.split(' ')
        list_max_coev[int(line[0])-1][int(line[1])-1] = float(line[4])
    return(np.amax(list_max_coev, axis=1))

def get_vector(structure,\
               path_file_pssm,\
               path_aaindex,\
               path_file_rsa,\
               path_file_asa,\
               list_dssp_features):
    """
    get the features vector for a structure, returna list of all the features
    for each residue
    """
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
                           ] + list_dssp_features[i])
                           #             
    return(list_vector)
def get_vector_neighbors(structure, list_vector):
    """
    compute the vector with neigbohrs features for each residue
    """
    list_vector_neighbors = []
    residues = structure.get_residues()
    for residue in residues:
        list_neighbors = get_neighbor_residues(structure, residue)
        vector_neighbors = []
        for neigbor in range(len(list_neighbors)):
            vector_neighbors = vector_neighbors +\
                               list_vector[list_neighbors[neigbor][0]]
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
        
    list_vector_1 = get_vector(structure_1[0],\
                               path_file_pssm,\
                                "",\
                               path_file_rsa,\
                               path_file_rsa)
    list_vector_2 = get_vector(structure_2[0],\
                               path_file_pssm,\
                                "",\
                               path_file_rsa,\
                               path_file_rsa)
    list_vector_neighbors_1 = get_vector_neighbors(structure_1[0],\
                                                   list_vector_1)
    list_vector_neighbors_2 = get_vector_neighbors(structure_2[0],\
                                                   list_vector_2)
    array_vector = []
    for i in range(len(sequence_1)):
        for j in range(len(sequence_2)):
            array_vector.append(list_vector_neighbors_1[i]+\
                                list_vector_neighbors_2[j])
    return(array_vector)
def get_voxel_data(path_bound,\
                   list_bound_pdb_file, k, voxel_size):
    """
    Compute a voxel around each resiude of the structure and compute 
    the voxelization
    """
    tut_data = home(dataDir='/home/sdv/m2bi/gollitrault/M2BI/projet_long/src')
    boxsize = [voxel_size,voxel_size,voxel_size]
    list_prot_vox = []
    parser = PDBParser()
    for i in range(k):

        structure_1 = parser.get_structure('test_bound_1',\
                                            path_bound + \
                                            "/templates/" +\
                                            list_bound_pdb_file[i][0:-1]+\
                                            '_1.pdb')
        structure_2 = parser.get_structure('test_bound_2',\
                                            path_bound +\
                                            "/templates/" +\
                                            list_bound_pdb_file[i][0:-1]+\
                                            '_2.pdb')
        
        residues_1 = [r for r in structure_1.get_residues()]
        residues_2 = [r for r in structure_1.get_residues()]
        if len(residues_1) == len(get_sequence(structure_1)) and\
           len(residues_2) == len(get_sequence(structure_2)):
            #vox
            try:
                prot_1 = Molecule(os.path.join(tut_data,\
                                               path_bound +\
                                               "/templates/" +\
                                               list_bound_pdb_file[i][0:-1]+\
                                               '_1.pdb'))
            except:
                return(list_prot_vox)
            try:
                prot_2 = Molecule(os.path.join(tut_data,\
                                               path_bound +\
                                               "/templates/" +\
                                               list_bound_pdb_file[i][0:-1]+\
                                               '_2.pdb'))
            except:
                return(list_prot_vox)
            
            try:
                prot_1 = prepareProteinForAtomtyping(prot_1)
            except:
                return(list_prot_vox)
            
            try:
                prot_2 = prepareProteinForAtomtyping(prot_2)
            except:
                return(list_prot_vox)
            #prot.view(guessBonds=False)

            for i in range(len(residues_1)):#len(residues_1)
                try:
                    prot_vox, prot_centers, prot_N = getVoxelDescriptors(\
                              prot_1,\
                              boxsize = boxsize,\
                              center = list(residues_1[i]["CA"].get_vector()),\
                              validitychecks=False)
                except:
                    return(list_prot_vox)
                list_prot_vox.append(prot_vox)
                
            for i in range(len(residues_1)):#len(residues_2)
                try:
                    prot_vox, prot_centers, prot_N = getVoxelDescriptors(\
                              prot_2,\
                              boxsize = boxsize,\
                              center = list(residues_2[i]["CA"].get_vector()),\
                              validitychecks=False)
                except:
                    return(list_prot_vox)
                list_prot_vox.append(prot_vox)
            make_voxel_npy_data(list_prot_vox)
        else:
            print("BAD number of residues do not correspond dont know why")
    
    return(list_prot_vox)
       
def get_X_Y(path_bound,\
            path_list_bound_pdb_file,\
            path_file_rsa,\
            path_file_rsa_bind,\
            path_file_asa,\
            path_file_asa_bind,\
            path_file_pssm,\
            path_aaindex,\
            k,\
            X_voxel):
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
    index_voxel_add = 0
    list_index_voxel_add = []
    for i in range(k):
        s=0
        structure_1 = parser.get_structure('test_bound_1',\
                                            path_bound +\
                                           "/templates/" +\
                                            list_bound_pdb_file[i][0:-1]+\
                                           '_1.pdb')
        structure_2 = parser.get_structure('test_bound_2',\
                                            path_bound +\
                                           "/templates/" +\
                                            list_bound_pdb_file[i][0:-1]+\
                                           '_2.pdb')
        structure_12 = parser.get_structure('test_bound_12',\
                                             path_bound +\
                                            "/templates_bind/" +\
                                             list_bound_pdb_file[i][0:-1]+\
                                            '.pdb')
        
        residues_1 = [r for r in structure_1.get_residues()]
        residues_2 = [r for r in structure_1.get_residues()]
        if len(residues_1) == len(get_sequence(structure_1)) and\
           len(residues_2) == len(get_sequence(structure_2)):
            
            list_value_rsa_1 = get_rsa_relative(path_file_rsa +\
                                                list_bound_pdb_file[i][0:-1]+\
                                                '_1.rsa')
            list_value_rsa_2 = get_rsa_relative(path_file_rsa +\
                                                list_bound_pdb_file[i][0:-1]+\
                                                '_2.rsa')
            
            list_relative_rsa_bind_1, list_relative_rsa_bind_2 =\
                                    get_rsa_relative_bind(path_file_rsa_bind +\
                                    list_bound_pdb_file[i][0:-1]+\
                                    '.rsa')
            
            list_dssp_features_1 = get_SS(structure_1,\
                                          path_bound +\
                                          "/templates/" +\
                                          list_bound_pdb_file[i][0:-1]+\
                                          '_1.pdb')
            list_dssp_features_2 = get_SS(structure_2,\
                                          path_bound +\
                                          "/templates/" +\
                                          list_bound_pdb_file[i][0:-1]+\
                                          '_2.pdb')

            list_surface_residue_1 = get_surface_residue(list_value_rsa_1,\
                                                         threshold_rsa = 30)
            list_surface_residue_2 = get_surface_residue(list_value_rsa_2,\
                                                         threshold_rsa = 30)

            list_interface_residue_1 = get_interface_residue(list_value_rsa_1,\
                                                      list_relative_rsa_bind_1)
            list_interface_residue_2 = get_interface_residue(list_value_rsa_2,\
                                                      list_relative_rsa_bind_2)

            list_vector_1 = get_vector(structure_1,\
                                       path_file_pssm +\
                                       list_bound_pdb_file[i][0:-1] +\
                                       '_1.fasta.pssm',\
                                       path_aaindex,\
                                       path_file_rsa +\
                                       list_bound_pdb_file[i][0:-1] +\
                                       '_1.rsa',\
                                       path_file_rsa +\
                                       list_bound_pdb_file[i][0:-1] +\
                                       '_1.rsa',\
                                       list_dssp_features_1)
            list_vector_2 = get_vector(structure_2,\
                                       path_file_pssm +\
                                       list_bound_pdb_file[i][0:-1] +\
                                       '_2.fasta.pssm',\
                                       path_aaindex,\
                                       path_file_rsa +\
                                       list_bound_pdb_file[i][0:-1] +\
                                       '_2.rsa',\
                                       path_file_rsa +\
                                       list_bound_pdb_file[i][0:-1] +\
                                       '_2.rsa',\
                                       list_dssp_features_2)

            list_vector_neighbors_1 = get_vector_neighbors(structure_1,\
                                                           list_vector_1)
            list_vector_neighbors_2 = get_vector_neighbors(structure_2,\
                                                           list_vector_2)

            list_vector_neighbors_1_surface = []
            list_vector_neighbors_2_surface = []
            
            X_voxel = data = load('../data/voxel_data/voxel_data.npy')

            for i in range(len(residues_1)):
                if i in list_interface_residue_1:
                    list_vector_neighbors_1_surface.append(\
                                                    list_vector_neighbors_1[i])
                    Y.append(1)
                    list_index_voxel_add.append(index_voxel_add)

                else:
                    if i in list_surface_residue_1:
                        list_vector_neighbors_1_surface.append(\
                                                    list_vector_neighbors_1[i])
                        Y.append(0)
                        list_index_voxel_add.append(index_voxel_add)
                index_voxel_add += 1


            for i in range(len(residues_2)):
                if i in list_interface_residue_2:
                    list_vector_neighbors_2_surface.append(\
                                                    list_vector_neighbors_2[i])
                    Y.append(1)
                    list_index_voxel_add.append(index_voxel_add)

                else:
                    if i in list_surface_residue_2:
                        list_vector_neighbors_2_surface.append(\
                                                    list_vector_neighbors_2[i])
                        Y.append(0)
                        list_index_voxel_add.append(index_voxel_add)
                index_voxel_add += 1
            X = X + list_vector_neighbors_1_surface +\
                    list_vector_neighbors_2_surface

        else:
            print("BAD number of residues, does not correspond dont know why")
    X_voxel_final = []
    for i in list_index_voxel_add:
        X_voxel_final.append(X_voxel[i])
    Y = to_categorical(Y, num_classes=2)
    X = np.asarray(X)
    
    X_voxel_final = np.asarray(X_voxel_final)
    print(len(X_voxel_final))
    print(len(X))
    print(len(Y))
    return([X, X_voxel_final],Y)

def evaluate_model(X_test, Y_test, model):
    Y_proba = model.predict(X_test)
    Y_pred = []
    Y_test_true = []
    threshold = 0.5
    for i in range(len(Y_test)):
        Y_test_true.append(Y_test[i][1])
        if Y_proba[i][1] > threshold:
            Y_pred.append(1)
        else:
            Y_pred.append(0)

    TN, FP, FN, TP = confusion_matrix(Y_test_true, Y_pred).ravel()
    # Sensitivity, hit rate, recall, or true positive rate
    TPR = TP/(TP+FN)
    # Specificity or true negative rate
    TNR = TN/(TN+FP) 
    # Precision or positive predictive value
    PPV = TP/(TP+FP)
    # Negative predictive value
    NPV = TN/(TN+FN)
    # Fall out or false positive rate
    FPR = FP/(FP+TN)
    # False negative rate
    FNR = FN/(TP+FN)
    # False discovery rate
    FDR = FP/(TP+FP)

    # Overall accuracy
    ACC = (TP+TN)/(TP+FP+FN+TN)
    return(ACC, TPR, TNR, FPR, FNR, confusion_matrix(Y_test_true, Y_pred))


    
def train_DeepNN_model(X_train, Y_train):
    model = DeepNN_model_build.build()
    model.compile(loss = "binary_crossentropy",optimizer="adam",\
                                               metrics=['accuracy'])
    model.fit(X_train, Y_train,epochs=20, batch_size = 64)
    return(model)
    
def data_roc_curve(X_test, Y_test, model, args):
    Y_proba = model.predict(X_test)
    Y_pred = []
    Y_test_true = []
    for i in range(len(Y_test)):
        Y_test_true.append(Y_test[i][1])
        Y_pred.append(Y_proba[i][1])

    Y_proba = model.predict(X_test)
    fpr, tpr, thresholds = roc_curve(Y_test_true, Y_pred)
    fpr_random, tpr_random, thresholds = roc_curve(Y_test_true,\
                                                   np.random.randint(2,\
                                                    size=len(Y_test_true)))
    model_auc = auc(fpr, tpr)
    args = True
    #if args == True:
    plt.plot(fpr,tpr, label='model')
    plt.plot(fpr_random,tpr_random, label='random')
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.legend()
    plt.savefig('roc_curve.png')
    plt.show()
    return(fpr, tpr, model_auc)

def make_fasta(list_bound_pdb_file, path_bound, path_fasta = "../data/fasta/"):
    parser = PDBParser()
    for i in range(len(list_bound_pdb_file)):
        structure_1 = parser.get_structure(list_bound_pdb_file[i][0:-1]+\
                                           '_1',\
                                           path_bound +\
                                           "/templates/" +\
                                           list_bound_pdb_file[i][0:-1]+\
                                           '_1.pdb')
        structure_2 = parser.get_structure(list_bound_pdb_file[i][0:-1]+\
                                           '_2',\
                                           path_bound +\
                                           "/templates/" +\
                                           list_bound_pdb_file[i][0:-1]+\
                                           '_2.pdb')
        
        file_fasta_1 = open(path_fasta + list_bound_pdb_file[i][0:-1]+\
                            '_1.fasta', 'w')
        file_fasta_1.write('>'+list_bound_pdb_file[i][0:-1]+ '_1'+"\n")
        file_fasta_1.write(''.join(list(get_sequence(structure_1))))
        file_fasta_1.close()
        file_fasta_2 = open(path_fasta + list_bound_pdb_file[i][0:-1]+\
                            '_2.fasta', 'w')
        file_fasta_2.write('>'+list_bound_pdb_file[i][0:-1]+ '_2'+"\n")
        file_fasta_2.write(''.join(list(get_sequence(structure_2))))
        file_fasta_2.close()

def make_voxel_npy_data(list_prot_vox):
    data = asarray(list_prot_vox)
    save('../data/voxel_data/voxel_data.npy', data)

def get_voxel_test(test_pdb, k, voxel_size):
    """
    Compute a voxel around each resiude of the structure and compute 
    the voxelization
    """
    tut_data = home(dataDir='/home/sdv/m2bi/gollitrault/M2BI/projet_long/src')
    boxsize = [voxel_size,voxel_size,voxel_size]
    list_prot_vox = []
    parser = PDBParser()
    structure = parser.get_structure('test', test_pdb)
    residues = [r for r in structure.get_residues()]
    if len(residues) == len(get_sequence(structure)):
    
        prot = Molecule(os.path.join(tut_data, test_pdb))
        prot = prepareProteinForAtomtyping(prot)
        for i in range(len(residues)):#len(residues_1)
            prot_vox, prot_centers, prot_N = getVoxelDescriptors(\
                          prot,\
                          boxsize = boxsize,\
                          center = list(residues[i]["CA"].get_vector()),\
                          validitychecks=False)
            list_prot_vox.append(prot_vox)
    else:
        print("BAD number of residues do not correspond dont know why")
    return(list_prot_vox)

#### MAIN ####
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    
    parser.add_argument("-path_bound", nargs="?",\
        help="path to the dockground folder benchmark",\
        const="../data/data_struct3d_bound")
        
    parser.add_argument("-path_list_bound_pdb_file", nargs="?",\
        help="path to the file/folder list of bound structures of the"+\
             " dockground benchmark",\
        const="../data/data_struct3d_bound/full_structures.0.90.txt")
        
    parser.add_argument("-path_file_rsa", nargs="?",\
        help="path to the file/folder with the rsa values for each residues"+\
             " for the structure alone  given by NACCESS",\
        const="../data/data_struct3d_bound/templates")
    parser.add_argument("-path_file_rsa_bind", nargs="?",\
        help="path to the file/folder with the rsa values for each residues"+\
             " for the bind structure given by NACCESS",\
        const="../data/data_struct3d_bound/templates_bind/")
        
    parser.add_argument("-path_file_asa", nargs="?",\
        help="path to the file/folder with the asa values for each residues"+\
             " for the structure alone  given by NACCESS",\
        const="../data/data_struct3d_bound/templates")
    parser.add_argument("-path_file_asa_bind", nargs="?",\
        help="path to the file/folder with the rsa values for each residues"+\
             " for the bind structure given by NACCESS",\
        const="../data/data_struct3d_bound/templates_bind/")
        
    parser.add_argument("-path_file_pssm", nargs="?",\
        help="path to the file/folder wit the pssm generated by psi blast",\
        const="../data/pssm/")
        
    parser.add_argument("-path_aaindex", nargs="?",\
        help="path to the folder wit the aaindex files",\
        const="../data/aaindex")
        
    parser.add_argument("-roc_plot", \
        action="store_true",\
        help="plot the roc curve and save it as roc_curve.png")

    parser.add_argument("-train",\
        type=int,\
        help="train the model on k number of proteins and output its "+\
             "carachteristics, k should be the same as k of voxel or lower")
    parser.add_argument("-test", nargs="?",\
        const="store_true",\
        help="test a protein")
    parser.add_argument("-voxel",\
        type=int,\
        help="compute the voxelization data for k proteins")
    parser.add_argument("-voxel_size", nargs="?",\
        help="chosse the box size around the CA atom",\
        default=10)
    parser.add_argument("-model_name", nargs="?",\
        help="chose the model name",\
        default="final_model.h5")
        
    args = parser.parse_args()
    
    
    #PATHS
    if args.path_bound:
        path_bound = args.path_bound
    else:
        path_bound = "../data/data_struct3d_bound"
    if args.path_list_bound_pdb_file:
        path_list_bound_pdb_file = args.path_list_bound_pdb_file
    else:
        path_list_bound_pdb_file = path_bound + "/full_structures.0.90.txt"
    if args.path_file_rsa:
        path_file_rsa = args.path_file_rsa
    else:
        path_file_rsa = path_bound + "/templates/"
    if args.path_file_rsa_bind:
        path_file_rsa_bind = args.path_file_rsa_bind
    else:
        path_file_rsa_bind = path_bound + "/templates_bind/"
    if args.path_file_asa:
        path_file_asa = args.path_file_asa
    else:
        path_file_asa = path_bound + "/templates/"
    if args.path_file_asa_bind:
        path_file_asa_bind = args.path_file_asa_bind
    else:
        path_file_asa_bind = path_bound + "/templates_bind/"
    if args.path_file_pssm:
        path_file_pssm = args.path_file_pssm
    else:
        path_file_pssm = "../data/pssm/"
    if args.path_aaindex:
        path_aaindex = args.path_aaindex
    else:
        path_aaindex = "../data/pssm/"
    if args.voxel_size:
        voxel_size = args.voxel_size
    else:
        voxel_size = 10

    #create voxel data
    if args.voxel:
        with open(path_list_bound_pdb_file) as file:
            list_bound_pdb_file = file.readlines()
        list_prot_vox = get_voxel_data(path_bound,\
                                       list_bound_pdb_file,\
                                       args.voxel,\
                                       voxel_size)
        make_voxel_npy_data(list_prot_vox)

    #

    if args.train:
        with open(path_list_bound_pdb_file) as file:
            list_bound_pdb_file = file.readlines()
        X_voxel = load('../data/voxel_data/Test.npy')
        X, Y = get_X_Y(path_bound,\
                       path_list_bound_pdb_file,\
                       path_file_rsa,\
                       path_file_rsa_bind,\
                       path_file_asa,\
                       path_file_asa_bind,\
                       path_file_pssm,\
                       path_aaindex,\
                       args.train,
                       X_voxel)
        nchannels = X[1].shape[2]
        X[1] = X[1].transpose().reshape([len(X[1]),\
                                        nchannels,\
                                        voxel_size,\
                                        voxel_size,\
                                        voxel_size])
        X[1] = np.asarray(X[1])
        
        X_train = [X[0][:int(len(X[0])*2/3)], X[1][:int(len(X[0])*2/3)]]
        X_test = [X[0][int(len(X[0])*2/3):], X[1][int(len(X[0])*2/3):]]
        Y_train = Y[:int(len(Y)*2/3)]
        Y_test = Y[int(len(Y)*2/3):]
        
        model = train_DeepNN_model(X_train, Y_train)
        ACC, TPR, TNR, FPR, FNR, conf_mat = evaluate_model(X_test,\
                                                           Y_test,\
                                                           model)
        ACC_train,\
        TPR_train,\
        TNR_train,\
        FPR_train,\
        FNR_train,\
        conf_mat_train = evaluate_model(X_train, Y_train, model)
        fpr, tpr, auc = data_roc_curve(X_test, Y_test, model, True)
        model.save_weights(filepath=args.model_name)
        
        n_1 = 0
        for i in range(len(Y)):
            if Y[i][1] == 1:
                n_1 += 1
        #wrtie file
        myfile = open("stats.txt", 'w')
        myfile.write("Statistics of the model:\n")
        myfile.write("Number of interface residues\n")
        myfile.write("\nN Total interface\n")
        myfile.write(str(len(Y)))
        myfile.write("\nn_Y_1\n")
        myfile.write(str(n_1))
        myfile.write("\nSize X_train:\n")
        myfile.write(str(len(X_train[0])))
        myfile.write("\nSize X_test:\n")
        myfile.write(str(len(X_test[0])))
        myfile.write("\nSize Y_train:\n")
        myfile.write(str(len(Y_train)))
        myfile.write("\nSize Y_test:\n")
        myfile.write(str(len(Y_test)))
        myfile.write("\nAccuracy on Test:\n")
        myfile.write("ACC: "+str(ACC)+" TPR:"+str(TPR)+" TNR: "+str(TNR)+\
                     " FPR: "+str(FPR)+" FNR: "+str(FNR)+"\n")
        myfile.write(str(conf_mat))
        myfile.write("\nAccuracy on Train:\n")
        myfile.write("ACC: "+str(ACC_train)+" TPR:"+str(TPR_train)+\
                     " TNR: "+str(TNR_train)+" FPR: "+str(FPR_train)+\
                     " FNR: "+str(FNR_train)+"\n")
        myfile.write(str(conf_mat_train))
        myfile.write("\n")
        myfile.close()
        #
    test_pdb = "../data/test_file/12as0A0B_1.pdb"
    test_rsa = "../data/test_file/12as0A0B_1.rsa"
    test_pssm = "../data/test_file/12as0A0B_1.fasta.pssm"
    if args.test:
        list_prot_vox = get_voxel_test(test_pdb, args.voxel, voxel_size)
        list_prot_vox = asarray(list_prot_vox)
        save('../data/test_file/voxel_test.npy', list_prot_vox)
        structure = parser.get_structure('test', test_pdb)
        list_dssp_features = get_SS(structure, test_pdb)
        list_vector = get_vector(structure,\
                                   test_pssm,\
                                   path_aaindex,\
                                   test_rsa,\
                                   test_rsa,\
                                   list_dssp_features)

        list_vector_neighbors = get_vector_neighbors(structure, list_vector)
        model = DeepNN_model_build.build()
        model.load_weights("final_model.h5")
        X = [list_vector_neighbors,list_prot_vox]
        nchannels = X[1].shape[2]
        X[1] = X[1].transpose().reshape([len(X[1]), nchannels, voxel_size,\
         voxel_size, voxel_size])
        X[1] = np.asarray(X[1])
        print(len(X[0]))
        print(len(X[1]))
        print("yep")
    




