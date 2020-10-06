#!/usr/bin/env python
# coding: utf-8

# # Import data from <i>.npys</i> preprocessed files and time-transform them

# In[ ]:


'''
Import continuos data from created .npy raw file (same as continuos_epochs), but save the Morlet transformed variables.

Input:
    ('/SERT/SERTXXXX/npys/structure.npy')
    
    
Output:
    np.save('/SERT/SERTXXXX/npys/structure_morlet.npy')


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
import pickle


# In[ ]:


print(wl.__file__)


# ## Downsampling and band-pass parameters

# <ol>
#     <li>Downsample to a resolution of 1 KHz.</li>
#     <li>Butter low-pass at 300 Hz, N = 9</li>
# </ol>

# In[2]:


### Parameters
# Downsampling parameters: final resolution = 1000 Hz
fs = 30000
final_fs  = 1000.0
task_length = 10 # min
OF_points = fs * task_length * 60
#ds_factor = fs // final_fs


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


# In[3]:


# Morlet parameters
dt = 1 / final_fs
time_windows = np.arange(0, 600, dt) # 600 por numero segundos en los 10 min de OF
frequencies = np.arange(1, 100, 1)
periods = 1 / (frequencies * dt)
scales = periods / wl.Morlet.fourierwl
n_frequencies = frequencies.shape[0]
time_points = time_windows.shape[0]


# In[4]:


ID = 'SERT1597'
npys_dir  = '/home/maspe/filer/SERT/' + ID + '/npys/'
#data = np.load(npys_dir + 'mPFC.npy')

with open(npys_dir + ID + '.info', 'rb') as f:
        info = pickle.load(f, encoding='latin1')

        
start_OF = info['startOF']
stop_OF = info['stopOF']

#OF_points = fs * 60 * 10


# In[18]:


##### Main loop #####
# Loop for loading and low-pass all channels of this mice

structures = ['mPFC'] #, 'NAC', 'BLA', 'vHip']

for structure in structures:
#iteration = 0
#for this_file in files:
    clock = time.time()
    print('Loading structure...')
    data = np.load(npys_dir + structure + '.npy')
    #print('Low-pass filtering (order = {}) at {} Hz...'.format(N, highcut))
    #data = channel['data'] #[start_OF - points_before : stop_OF]
    #data = signal.filtfilt(b=b, a=a, x=data - np.mean(data),
    #                       axis=-1, padtype='odd', padlen=None, method='pad', irlen=None)
    #if iteration == 0:
    #    data_matrix = np.empty((n_channels, len(data)))         
    #data_matrix[iteration, :] = data    
    print('Loaded in {:.2f} min.'.format((time.time() - clock) / 60))
    
    
    clock = time.time()
    print('Downsampling...')
    data = signal.resample(x=data[:, start_OF : start_OF + OF_points], num=time_points, axis=1)
    print('Downsampled in {:.2f} min.'.format((time.time() - clock) / 60))

    clock = time.time()
    print('Morlet transform...')
    #if iteration == 0:
    #    morlet_matrix = np.empty((n_frequencies, len(data), n_channels))
    
    n_channels = data.shape[0]
    morlet_matrix = np.empty((n_frequencies, time_points, n_channels))
    for channel in range(n_channels):
        morlet_matrix[:, :, channel] = wl.Morlet(data[channel, :], scales=scales).getnormpower()
    #morlet_matrix[:, :, iteration] = transformed
    #print('Transformed in {:.2f} min.\n'.format((time.time() - clock) / 60))
    
    #iteration += 1
    
print('Done!')


# Takes as input $\text{data}_{(n\_channels, time\_points)}$ with downsampled $\text{time points} = 1 KHz \times 6 s = 600000$. Then, for each $data[channel, :] = (600000, 1)$ transforms it to the frequency domain with wl.Morlet, getting a <i>(n_frequencies, time_points)</i> morlet_matrix for each channel. 

# In[14]:


#data[channel, :].shape


# In[15]:


#d=wl.Morlet(data[channel, :], scales=scales).getnormpower()


# In[16]:


#d.shape


# In[13]:


#morlet_matrix[channel, :, :].shape


# In[ ]:


#print('\nCollecting all channels and Morlets by structure...')    

#print('mPFC')
#mPFC_morlet = morlet_matrix[:, :, mPFC_indexes]
#print('Saving variable (this takes some time)...')
#np.save(npys_dir + 'mPFC_morlet', mPFC_morlet)
#print('mPFC: Done!')


#print('mPFC')
#NAC_morlet = morlet_matrix[:, :, NAC_indexes]
#print('Saving variable (this takes some time)...')
#np.save(npys_dir + 'NAC_morlet', NAC_morlet)
#print('NAC: Done!')


#print('BLA')
#BLA_morlet = morlet_matrix[:, :, BLA_indexes]
#print('Saving variable (this takes some time)...')
#np.save(npys_dir + 'BLA_morlet', BLA_morlet)
#print('BLA: Done!')


#print('vHip')
#vHip_morlet = morlet_matrix[:, :, vHip_indexes]
#print('Saving variable (this takes some time)...')
#np.save(npys_dir + 'vHip_morlet', vHip_morlet)
#print('vHip: Done!')

#del [data, transformed, data_matrix, morlet_matrix]  

#print('Mouse processed in {:.2f} min.\n'.format((time.time() - master_clock) / 60))
#print('Done!\n\n')

