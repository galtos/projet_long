# -*- coding: utf-8 -*-
"""
Created on Sat Dec 28 17:38:58 2019

@author: Guillaume
"""

import numpy as np
import pandas 
from numpy import loadtxt
from keras.models import Sequential
from keras.layers import Dense, Dropout, Activation, Flatten
from keras.layers.convolutional import Convolution3D
from keras.layers.convolutional import MaxPooling3D
from keras.utils import np_utils
from keras.utils import to_categorical
from keras.layers.normalization import BatchNormalization
from keras.layers import Input
from keras.layers import concatenate
from keras.models import Model
#scikit-learn
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import KFold
from sklearn.preprocessing import LabelEncoder

class DeepNN_model_build:
    def build():
        inputA = Input(shape=(189,))
        inputB = Input(shape=(8, 10, 10,10,))
        
        x = Dense(162, activation = "relu")(inputA)
        x = BatchNormalization()(x)
        x = Dense(20, activation = "relu")(x)
        x = Model(inputs=inputA, outputs=x)
        
        y = Convolution3D(\
            filters=16,\
            kernel_size=2,\
            padding='valid',\
            data_format='channels_first',\
            )(inputB)
        y = Dropout(0.2)(y)
        y = MaxPooling3D(\
            pool_size=(2,2,2),\
            strides=None,\
            padding='valid',\
            data_format='channels_first'
            )(y)
        y = Flatten()(y)
        y = Dense(20, activation = "relu")(y)
        y = Model(inputs=inputB, outputs=y)
        
        
        combined = concatenate([x.output, y.output])
        
        z = Dense(4, activation="relu")(combined)
        z = Dense(2, activation="softmax")(z)
        
        model = Model(inputs=[x.input, y.input], outputs=z)
        return(model)
        
        
        
        
        
        
        
        
