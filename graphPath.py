#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 21:42:23 2020

@author: diogenes
"""

#%%
import os
import pickle
import pandas as pd

import matplotlib.pyplot as plt


#%%
path_to_files = '/home/diogenes/Projects/brainsignals/results/info_files/'
mice_list = sorted(os.listdir(path_to_files))

all_mice = pd.DataFrame()
for mouse in mice_list:
    temp = pickle.load(open(path_to_files + mouse, 'rb'))
    
    all_mice = all_mice.append(temp['path'], ignore_index=True)
    #print(temp['path'])

all_mice = all_mice.dropna()   
print(all_mice)
print("Shape {}".format(all_mice.shape))    


#%%

plt.hist(all_mice['dwell'])