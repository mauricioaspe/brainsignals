#!/usr/bin/env python
# coding: utf-8

# # Import data from <i>.npys</i> epoched files and time-transform them

# In[ ]:


'''
Import epoched data from created .npy raw file (same as continuos_epochs) (??)

Input:
    * Epoched data stored in
    'SERT/SERTXXXX/npys/baseline_epochs.npy'

    
Output:
    np.save('/SERT/SERTXXXX/npys/morlets_epochs.npy')


'''


# #### Notice that must take the npys arrays given by 'continuous_epochs' and with epochs removed!!
# #### Import required modules

# In[1]:


#!/usr/bin/python
# -*- coding: utf-8 -*-

# Import required modules
import numpy as np
from scipy import signal
from matplotlib import pyplot as plt
import wavelets as wl
import time
import pickle


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


# ### Revisar que duracion baseline vs epocas no interfiera en la transformada!

# In[3]:


# Morlet parameters
epoch_length = 6 # In seconds
dt = 1 / final_fs
time_windows = np.arange(0, epoch_length, dt)
frequencies = np.arange(1, 100, 1)
periods = 1 / (frequencies * dt)
scales = periods / wl.Morlet.fourierwl
n_frequencies = frequencies.shape[0]


# In[4]:


IDs = {'SERT1597': {},
       'SERT1659': {},
       'SERT1665': {},
       'SERT1668': {},
       'SERT1678': {},
       'SERT1908': {},
       'SERT1984': {},
       'SERT1985': {},
       'SERT2013': {},
       'SERT2014': {},
       'SERT2018': {},
       'SERT2024': {}}


# ## Main loop

# Takes as input $\text{data}_{(n\_channels, time\_points)}$ with downsampled $\text{time points} = 1 KHz \times 6 s = 600000$. Then, for each $data[channel, :] = (600000, 1)$ transforms it to the frequency domain with wl.Morlet, getting a <i>(n_frequencies, time_points)</i> morlet_matrix for each channel. 

# In[11]:


structures = {'mPFC': {}, 'NAC': {}, 'BLA': {}, 'vHip': {}}

iteration = 1
data_dict = {}
master_clock = time.time()
for mouse in IDs.keys():
    clock = time.time()
    print('Loading mouse {} (#{})...'.format(mouse, iteration))
    npys_dir  = '/home/maspe/filer/SERT/' + mouse + '/npys/'

    with open(npys_dir + mouse + '.info', 'rb') as f:
        info = pickle.load(f, encoding='latin1')

    ### Loading data    
    all_data = pickle.load(open(npys_dir + mouse + '.epochs', 'rb'), encoding="latin1")              

    
    for structure in structures.keys():
        print('Loading {}...'.format(structure))
        for condition in all_data[structure].keys(): # Iterates in baselines and epochs
            if condition == 'epochs':
                print('Processing epochs...')
                time_points = time_windows.shape[0]
            else:
                print('Processing baselines...')
                time_points = time_windows.shape[0] // 2
            
            ### Loading the data
            data = all_data[structure][condition]
        
            ### Downsampling
            print('Downsampling...')
            data = signal.resample(x=data, num=time_points, axis=1)

            ### Morlet transform
            print('Morlet transform...')
            n_channels = data.shape[0]
            n_epochs = data.shape[2]
            morlet_matrix = np.empty((n_frequencies, time_points, n_channels, n_epochs))  
            for epoch in range(n_epochs):
                for channel in range(n_channels):
                    morlet_matrix[:, :, channel, epoch] = wl.Morlet(data[channel, :, epoch], scales=scales).getnormpower()
            
            
            if condition == 'epochs':
                morlets_epochs    = morlet_matrix
            else:
                morlets_baselines = morlet_matrix
  
        
        ### Obtaining z-scores
        print('Transforming to z-score...\n')
        
        ### Getting average and sd on time dimension
        baseline_mean = np.mean(morlets_baselines, axis=(1,3))
        baseline_sd = np.std(morlets_baselines, axis=(1,3))
        morlets_epochs = (morlets_epochs - baseline_mean[:, None, :, None]) / baseline_sd[:, None, :, None]
    
        structures[structure] = np.mean(morlets_epochs, axis=3)
    
    
    iteration += 1
    
    print('Saving dictionary...')    
    pickle.dump(structures, open(npys_dir + 'morlets.epochs', 'wb'), protocol=2)
    print('Mouse processed in {:.2f} min.\n'.format((time.time() - clock) / 60))

print('All mice processed in {:.2f} min.\n'.format((time.time() - master_clock) / 60))    
print('Done!')


# ################
