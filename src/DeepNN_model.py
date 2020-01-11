# -*- coding: utf-8 -*-
"""
Created on Sat Dec 28 17:38:58 2019

@author: Guillaume
"""

import numpy as np
import pandas 
from numpy import loadtxt
from keras.models import Sequential
from keras.layers import Dense, Dropout, Activation
from keras.utils import np_utils
from keras.utils import to_categorical
from keras.layers.normalization import BatchNormalization
#scikit-learn
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import KFold
from sklearn.preprocessing import LabelEncoder

class DeepNN_model_build:
    def build():
        model = Sequential()
        model.add(Dense(162, input_dim = 162, activation = "relu"))
        model.add(BatchNormalization())
        model.add(Dropout(0.3))
        model.add(Dense(20,activation ="relu"))
        model.add(Dense(2, activation="relu"))
        
        return(model)
