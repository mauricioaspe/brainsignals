#!/usr/bin/env python
# coding: utf-8

'''
The file works from terminal taking the name of the mouse as its only argument.
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

'''

# Import required modules
import sys
import os
import numpy as np
import physig as ps
from scipy import signal
from matplotlib import pyplot as plt

#sys.path.insert(1, '/home/diogenes/Projects/SERT/modules/')


# Create the filter
# Butterpass at 300 Hz, oder 9
def butter_bandpass(highcut, fs, order=5):
    nyq  = 0.5 * fs
    high = highcut / nyq
    b, a = signal.butter(order, high)
    return b, a

# Filter parameters
fs = 30000.0
highcut = 300.0
N       = 9
b, a    = butter_bandpass(highcut, fs, order=N)


# Main loop
# The loop
#     Load each of the 32 .continuous Openephys files</li>
#     Subtract the mean</li>
#     Low-pass it at 300 Hz</li>
#     Subtract the median by structure</li>
#     Save it as '/home/maspe/filer/SERT/SERTXXXX/npys/mPFC'</li>

#pulses_path = '/home/diogenes/Projects/brainsignals/DATA/MICE/1597/continuous/AUX/'

idS = ['SERT1597', 'SERT1659', 'SERT1678', 'SERT1908',
       'SERT1984', 'SERT1985', 'SERT2014', 'SERT1668',
       'SERT1665', 'SERT2018', 'SERT2024', 'SERT2013']


#%%
##### Main loop #####
# Loop for loading and low-pass all channels of this mice
iteration = 0
for this_file in files:
    #pulses = ps.loadContinuous(pulses_file)
    channel = ps.loadContinuous(this_file)

    #print("Extracting pulses...")
    #d = np.diff(np.where(pulses['data'] < -4))[0]
    #on = np.where(d > 250)[0]

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

