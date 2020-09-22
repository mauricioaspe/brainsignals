#!/usr/bin/env python
# coding: utf-8

# # Import continuous data from .continuous Openephys raw files
# 
# ### Ideas
# <ul>
#     <li>Crear una clase Mouse donde poner todo por ratón?</li>
#     <li>Definitivamente crear una clase con los filtros</li>
#     </ul>

# In[ ]:


""" The file works from terminal taking the name of the mouse as its only argument.
E.g. Run:

python3 continuous_epochs.py 'SERT1597'

to preprocess the mouse 'SERT1957'


Output: 
    1. Info files
    2. Matrix by structure, with 
        - Channels at 30 KHz
        - All channels with mean subtracted
        - Low-passed at 300 Hz
        - Median subtracted
        - Saved at SERTXXXX/npys_dir + structure.npy
        
""" 


# #### Import the required modules

# In[3]:


# Import required modules
import sys
sys.path.insert(1, '/home/diogenes/Projects/SERT/modules/')

import glob
import numpy as np
import pandas as pd
import physig as ps
from scipy import signal
from matplotlib import pyplot as plt



# #### Create the filter

# Butterpass at 300 Hz, oder 9

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


## Main loop

# The loop</b><br>
#     <li>Load each of the 32 .continuous Openephys files</li>
#     <li>Subtract the mean</li>
#     <li>Low-pass it at 300 Hz</li>
#     <li>Subtract the median by structure</li>
#     <li>Save it as '/home/maspe/filer/SERT/SERTXXXX/npys/mPFC'</li>


##### Main loop #####
# Loop for loading and low-pass all channels of this mice
iteration = 0
for this_file in files:
    channel = ps.loadContinuous(this_file)

    print('Low-pass filtering (order = {}) at {} Hz...'.format(N, highcut))
    data = channel['data'] #[start_OF - points_before : stop_OF]
    data = signal.filtfilt(b=b, a=a, x=data - np.mean(data),
                           axis=-1, padtype='odd', padlen=None, method='pad', irlen=None)

    
    if iteration == 0:
        data_matrix = np.empty((len(channels_locations), len(data)))          
    
    data_matrix[iteration, :] = data    
    iteration += 1

print('\nCollecting all channels by structure...')    
mPFC = data_matrix[mPFC_indexes, :] - np.median(data_matrix[mPFC_indexes, :], axis=0)
print('mPDF: Done!')
NAC  = data_matrix[NAC_indexes, :]  - np.median(data_matrix[NAC_indexes, :], axis=0)
print('NAC: Done!')
BLA  = data_matrix[BLA_indexes, :]  - np.median(data_matrix[BLA_indexes, :], axis=0)
print('BLA: Done!')
vHip = data_matrix[vHip_indexes, :] - np.median(data_matrix[vHip_indexes, :], axis=0)
print('vHip: Done!')

del [iteration, channel, data, data_matrix]  


if save:
    print('Saving variables...')
    np.save(npys_dir + 'mPFC', mPFC)
    np.save(npys_dir + 'NAC', NAC)
    np.save(npys_dir + 'BLA', BLA)
    np.save(npys_dir + 'vHip', vHip)


print('Done!')

