#!/usr/bin/env python
# coding: utf-8

# In[1]:


### Importing modules
import time
import pickle
import numpy as np
import wavelets as wl
from scipy import signal


# In[3]:


### Parameters
# Downsampling parameters: final resolution = 1000 Hz
fs = 30000.0
final_fs  = 1000.0
ds_factor = fs // final_fs

# 
npoints = 180000
n_samples = np.int(npoints / ds_factor) # Downsample to 1000 Hz
fs = 1000.0

# Morlet parameters
dt = 1 / fs
time_windows = np.arange(0, 3, dt)
frequencies = np.arange(1, 50, 1)
periods = 1 / (frequencies * dt)
scales = periods / wl.Morlet.fourierwl
n_frequencies = frequencies.shape[0]
time_points = time_windows.shape[0]

# Reducing epoch
start = 1000
stop = 4000
baseline = 2000


# In[ ]:


### Loading data


# In[4]:


IDs_WT = ['SERT1597'] #, 'SERT1659', 'SERT1678', 'SERT1908', 'SERT1984', 'SERT1985', 'SERT2014'] 
IDs_KO = ['SERT1665'] #, 'SERT1668', 'SERT2013', 'SERT2018', 'SERT2024'] 

WT = {'mPFC': {}, 'BLA': {}} #, 'NAC': {}, 'vHip': {}}
KO = {'mPFC': {}, 'BLA': {}} #, 'NAC': {}, 'vHip': {}}
clock = time.time()

print('###################\nDownsampling WTs to 1000 Hz...')
for structure in WT.keys():
    for ID in IDs_WT:
        print('Loading {} from {}...'.format(structure, ID))
        npys_dir = '/home/maspe/filer/SERT/' + ID + '/npys/'
        x = np.load(npys_dir + structure + '_epochs.npy', allow_pickle=True)      
    
        print('Downsampling...\n')
        WT[structure][ID] = signal.resample(x=x, num=n_samples, axis=1)

        
print('###################\nDownsampling KOs to 1000 Hz...')
for structure in KO.keys():
    for ID in IDs_KO:
        print(ID)
        print('Loading {} from {}...'.format(structure, ID))
        npys_dir = '/home/maspe/filer/SERT/' + ID + '/npys/'        
        x = np.load(npys_dir + structure + '_epochs.npy', allow_pickle=True)
    
        print('Downsampling...\n')
        KO[structure][ID] = signal.resample(x=x, num=n_samples, axis=1)

        
print('All mice downsampled in {:.2f} s.'.format(time.time() - clock))


# In[10]:


SRS_WT['BLA']['SERT1597'].shape


# In[8]:


# Morlet transform
SRS_WT = {'mPFC':{}, 'BLA':{}} #, 'NAC':{}, 'vHip':{}}
SRS_KO = {'mPFC':{}, 'BLA':{}} #, 'NAC':{}, 'vHip':{}}

# For WT mice
for structure in WT.keys():
    for mouse in WT[structure].keys():
        clock = time.time()
        print('Morlet wavelet for {} from {} started!'.format(structure, mouse))
    
        n_channels = WT[structure][mouse].shape[0]
        n_epochs = WT[structure][mouse].shape[2]

        SRS_WT[structure][mouse] = np.zeros((n_frequencies, time_points, n_channels, n_epochs))
    
        for epoch in range(n_epochs):
            for channel in range(n_channels): # SRP.shape[0]):
                ### DOCS says data: data in array to transform, length must be power of 2 !!!!
                wavel1 = wl.Morlet(WT[structure][mouse][channel, start:stop, epoch], scales=scales)
                SRS_WT[structure][mouse][:, :, channel, epoch] = wavel1.getnormpower()
            
        #print('Transforming to z-score...\n')        
        #baseline_mean = np.mean(SRS_WT[structure][mouse][:, baseline:, :, :], axis=1)
        #baseline_sd = np.std(SRS_WT[structure][mouse][:, baseline:, :, :], axis=1)
        #SRS_WT[structure][mouse] = (SRS_WT[structure][mouse] - baseline_mean[:, None, :, :]) / baseline_sd[:, None, :, :]

        
#print('Saving Morlets for WT...\n')
#pickle.dump(SRS_WT, open('/home/maspe/filer/SERT/ALL/npys/SRS_WT.dict', 'wb'), protocol=2)


# For KO mice
for structure in KO.keys():
    for mouse in KO[structure].keys():
        clock = time.time()
        print('Morlet wavelet for {} from {} started!'.format(structure, mouse))
    
        n_channels = KO[structure][mouse].shape[0]
        n_epochs = KO[structure][mouse].shape[2]

        SRS_KO[structure][mouse] = np.zeros((n_frequencies, time_points, n_channels, n_epochs))
    
        for epoch in range(n_epochs):
            for channel in range(n_channels): # SRP.shape[0]):
                ### DOCS says data: data in array to transform, length must be power of 2 !!!!
                wavel1 = wl.Morlet(KO[structure][mouse][channel, start:stop, epoch], scales=scales)
                SRS_KO[structure][mouse][:, :, channel, epoch] = wavel1.getnormpower()
            
        #print('Transforming to z-score...')        
        #baseline_mean = np.mean(SRS_KO[structure][mouse][:, baseline:, :, :], axis=1)
        #baseline_sd = np.std(SRS_KO[structure][mouse][:, baseline:, :, :], axis=1)
        #SRS_KO[structure][mouse] = (SRS_KO[structure][mouse] - baseline_mean[:, None, :, :]) / baseline_sd[:, None, :, :]

        
#print('Saving Morlets for KO...\n')
#pickle.dump(SRS_KO, open('/home/maspe/filer/SERT/ALL/npys/SRS_KO.dict', 'wb'), protocol=2)


print('All mice Fourier transformed in {:.2f} min.\n'.format((time.time() - clock) / 60))        
print('Done!')

