#!/usr/bin/env python
# coding: utf-8

# # Import data from .npys preprocessed files and time-transform them

# In[ ]:


'''
Import continuos data from .openephys raw file (same as continuos_epochs), but save the Morlet transformed variables.




'''


# #### Notice that must take the npys arrays given by 'continuous_epochs' and with epochs removed!!
# #### Import required modules

# In[1]:


# Import required modules
#import glob
#import sys
import numpy as np
#import pandas as pd
#import physig as ps
from scipy import signal
#from matplotlib import pyplot as plt
import wavelets as wl
import time
#import pickle


# In[2]:


ID = 'SERT1597'


# In[4]:


npys_dir  = '/home/maspe/filer/SERT/' + ID + '/npys/'
data = np.load(npys_dir + 'mPFC.npy')


# In[6]:


data.shape


# ## Downsampling and band-pass parameters

# <ol>
#     <li>Downsample to a resolution of 1 KHz.</li>
#     <li>Butter low-pass at 300 Hz, N = 9</li>
# </ol>

# In[ ]:


# Create filter
def butter_bandpass(highcut, fs, order=5):
    nyq  = 0.5 * fs
    high = highcut / nyq
    b, a = signal.butter(order, high)
    return b, a

# Filter parameters
highcut = 300.0
N       = 9
b, a    = butter_bandpass(highcut, fs, order=N)

### Parameters
# Downsampling parameters: final resolution = 1000 Hz
fs = 30000
final_fs  = 1000.0
OF_points = fs * 60 * 10
#ds_factor = fs // final_fs


# In[ ]:


# Morlet parameters
dt = 1 / final_fs
time_windows = np.arange(0, 600, dt) # 600 por numero segundos en los 10 min de OF
frequencies = np.arange(1, 100, 1)
periods = 1 / (frequencies * dt)
scales = periods / wl.Morlet.fourierwl
n_frequencies = frequencies.shape[0]
time_points = time_windows.shape[0]


# In[ ]:


##### Main loop #####
# Loop for loading and low-pass all channels of this mice

structures = ['mPFC', 'NAC', 'BLA', 'vHip']

iteration = 0
for this_file in files:
    data = ps.loadContinuous(this_file)

    #print('Low-pass filtering (order = {}) at {} Hz...'.format(N, highcut))
    #data = channel['data'] #[start_OF - points_before : stop_OF]
    #data = signal.filtfilt(b=b, a=a, x=data - np.mean(data),
    #                       axis=-1, padtype='odd', padlen=None, method='pad', irlen=None)


    if iteration == 0:
        data_matrix = np.empty((n_channels, len(data)))
              
    data_matrix[iteration, :] = data    
    
    clock = time.time()
    print('Downsampling...')
    data = signal.resample(x=data[start_OF : start_OF + OF_points], num=time_points)
    print('Downsampled in {:.2f} min.'.format((time.time() - clock) / 60))

    clock = time.time()
    print('Morlet transform...')
    if iteration == 0:
        morlet_matrix = np.empty((n_frequencies, len(data), n_channels))
    
    transformed = wl.Morlet(data, scales=scales).getnormpower()
    morlet_matrix[:, :, iteration] = transformed
    print('Transformed in {:.2f} min.\n'.format((time.time() - clock) / 60))
    
    iteration += 1

    
    
print('\nCollecting all channels and Morlets by structure...')    

print('mPFC')
mPFC_morlet = morlet_matrix[:, :, mPFC_indexes]
print('Saving variable (this takes some time)...')
np.save(npys_dir + 'mPFC_morlet', mPFC_morlet)
print('mPFC: Done!')


print('mPFC')
NAC_morlet = morlet_matrix[:, :, NAC_indexes]
print('Saving variable (this takes some time)...')
np.save(npys_dir + 'NAC_morlet', NAC_morlet)
print('NAC: Done!')


print('BLA')
BLA_morlet = morlet_matrix[:, :, BLA_indexes]
print('Saving variable (this takes some time)...')
np.save(npys_dir + 'BLA_morlet', BLA_morlet)
print('BLA: Done!')


print('vHip')
vHip_morlet = morlet_matrix[:, :, vHip_indexes]
print('Saving variable (this takes some time)...')
np.save(npys_dir + 'vHip_morlet', vHip_morlet)
print('vHip: Done!')

del [data, transformed, data_matrix, morlet_matrix]  

print('Mouse processed in {:.2f} min.\n'.format((time.time() - master_clock) / 60))
print('Done!\n\n')

