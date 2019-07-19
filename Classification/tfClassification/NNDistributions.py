#!/usr/bin/env python


#import tensorflow as tf

#print(tf.version.VERSION)



from __future__ import absolute_import, division, print_function, unicode_literals

# TensorFlow and tf.keras
import tensorflow as tf
from tensorflow import keras

# Helper libraries
import numpy as np
import matplotlib.pyplot as plt

import csv

# read in the quantiles of the distributions

import os

dir_path = os.path.dirname(os.path.realpath(__file__))
filename = dir_path + "/olympics2016.csv"

features = tf.placeholder(tf.int32, shape=[3], name='features')
country = tf.placeholder(tf.string, name='country')
total = tf.reduce_sum(features, name='total')



path = "/home/willy/PredictingProteinInteractions/data/animals/"
fName = "proj_n_500_m_100_q_10.csv"
#readInCsv(path+fName)
































































