#!/usr/bin/env python
# coding: utf-8

### Importing required modules
import os
import sys

### Configuring directories
#os.chdir('/home/maspe/filer/projects/brainsignals')
sys.path.append('/home/maspe/filer/projects/brainsignals/modules/')

import pickle
import time
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
import connectivity_measures as cm

abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)
print(os.getwd)

### List of animals to process
IDs = ['ID1597', 'ID1659']
#, 'SERT1678', 'SERT1908', 'SERT1984', 'SERT1985', 'SERT2014', 'SERT1668', 'SERT1665', 'SERT2018', 'SERT2024', 'SERT2013']
# ### Probar en ventanas de 500 ms!!

### Butterworth filter
fs = 1000.0
def butter_bandpass(lowcut, highcut, fs, order=5):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = signal.butter(order, [low, high], btype='band')
    return b, a

### Filter parameters. You can add more frequency bands (lowcuts and highcuts in Hz)
filter_parameters = {'theta': {'N': 4, 'lowcut': 3.5,  'highcut': 7.5}}

### Creating a filter for each frequency band
butterworths = dict()
for band in filter_parameters.keys():
    lowcut = filter_parameters[band]['lowcut']
    highcut = filter_parameters[band]['highcut']
    N = filter_parameters[band]['N']
    butterworths[band] = butter_bandpass(lowcut, highcut, fs, order=N)


iCohs = {key: {} for key in filter_parameters.keys()}
for ID in IDs: 
    epochs_dir = 'DATA/epochs/'
    print('\nLoading mouse {}...'.format(ID))
    
    ### Loading data
    data = np.load(epochs_dir + ID + '_epochs.npy')
    
    ### Loop
    #filtered_bands = {key: {} for key in filter_parameters.keys()}
    #iterator = 0
    
    # Filters the 32 channels x 180,000 time_points x n epochs matrix along the time_points axis 
    print('Filtering...')    
    for band in iCohs.keys():
        filtered = signal.filtfilt(b=butterworths[band][0], a=butterworths[band][1],
                                   x=data, axis=1)

    #print("Data shape: {}".format(data.shape))
    #print("Filtered shape: {}".format(filtered.shape))

    print('Calculating iCoh for {} band...'.format(band))
    clock = time.time()

    # Loop
    n_epochs = filtered.shape[2]
    time_points = filtered.shape[1]
    for epoch in range(n_epochs):
        if epoch == 0:
            icoh_pre  = cm.icoh(filtered[:,:int(time_points / 2),epoch], average = False)
            icoh_post  = cm.icoh(filtered[:,int(time_points / 2):,epoch], average = False)            
        else:
            icoh_pre  = np.dstack((icoh_pre, cm.icoh(filtered[:,:int(time_points / 2),epoch], average = False)))
            icoh_post  = np.dstack((icoh_pre, cm.icoh(filtered[:,int(time_points / 2):,epoch], average = False)))

    print('iCoh calculated in {:.2f} min'.format(time.time() - clock))
    print("ICoh shape (channels, channels, epochs): {}".format(icoh_pre.shape))

    iCohs[band]['pre'] = icoh_pre
    iCohs[band]['post'] = icoh_post

#pickle.dump(iCohs, open(npys_dir + mouse + '_BLA_vHip.icoh', 'wb'), protocol=2)

print('Done!')  


