#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import time
import numpy as np
import wavelets as wl
from scipy import signal
from matplotlib import pyplot as plt


# In[7]:


# Downsampling parameters: final resolution = 1000 Hz
fs = 30000.0
final_fs  = 1000.0
ds_factor = fs // final_fs

npoints = 180000
n_samples = np.int(npoints / ds_factor)

# Reducing epoch
start = 60000
stop = 90000


# In[3]:


IDs_WT = ['SERT1597', 'SERT1659'] #, 'SERT1678', 'SERT1908', 'SERT1984', 'SERT1985', 'SERT2014']
# IDs_KO = ['SERT1668', 'SERT1665', 'SERT2018', 'SERT2024', 'SERT2013'] 

# allFigs_dir = '/home/maspe/filer/SERT/ALL/figs/'

mPFC_WT = dict()
print('Processing WT')
iteration = 0
for ID in IDs_WT:
    clock = time.clock()
    npys_dir = '/home/maspe/filer/SERT/' + ID + '/npys/'
        
    mPFC_WT[ID] = np.load(npys_dir + 'mPFC_epochs.npy', allow_pickle=True)
    print('{} loaded!'.format(ID))
    
    
    print('Begin downsampling')
    mPFC_WT[ID] = signal.resample(x=mPFC_WT[ID], num=n_samples, axis=1)

    print('{} processed in {}'.format(ID, time.clock() - start))
    print('{} matrix shape = {}'.format(ID, mPFC_WT[ID].shape))
    
    #NAC_epochs  = np.mean(np.load(npys_dir + 'NAC_epochs.npy', allow_pickle=True), axis=2)
    #BLA_epochs  = np.mean(np.load(npys_dir + 'BLA_epochs.npy', allow_pickle=True), axis=2)
    #vHip_epochs = np.mean(np.load(npys_dir + 'vHip_epochs.npy', allow_pickle=True), axis=2)

