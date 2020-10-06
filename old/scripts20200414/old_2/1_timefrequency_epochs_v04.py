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

# In[8]:


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

# In[9]:


### Parameters
# Downsampling parameters: final resolution = 1000 Hz
fs = 30000
final_fs  = 1000.0


# ### Revisar que duracion baseline vs epocas no interfiera en la transformada!

# In[10]:


# Morlet parameters
epoch_length = 6 # In seconds
dt = 1 / final_fs
time_windows = np.arange(0, epoch_length, dt)
frequencies = np.arange(1, 100, 1)
periods = 1 / (frequencies * dt)
scales = periods / wl.Morlet.fourierwl
n_frequencies = frequencies.shape[0]


# In[11]:


IDs = {'SERT1597': {}, 'SERT1659': {} #,
 #'SERT1665': {},
 #'SERT1668': {},
 #'SERT1678': {},
 #'SERT1908': {},
 #'SERT1984': {},
 #'SERT1985': {},
 #'SERT2013': {},
 #'SERT2014': {},
 #'SERT2018': {},
 #'SERT2024': {}
      }


# ## Main loop

# In[ ]:


structures = {'mPFC': {}, 'NAC': {}} #, 'BLA': {}, 'vHip': {}}

for mouse in IDs.keys():
    clock = time.time()
    print('Loading mouse {}...'.format(mouse))
    npys_dir  = '/home/maspe/filer/SERT/' + mouse + '/npys/'

    with open(npys_dir + mouse + '.info', 'rb') as f:
        info = pickle.load(f, encoding='latin1')

        
    all_data = np.load(npys_dir + 'baseline_epochs.npy', allow_pickle=True)[()]

    
    for structure in structures.keys():
        for condition in all_data[structure].keys(): # Iterates in baselines and epochs
            if condition == 'epochs':
                print('Processing epochs...')
                time_points = time_windows.shape[0]
            else:
                print('Processing baselines...')
                time_points = time_windows.shape[0] // 2
            
            #print(time_points)
            data = all_data[structure][condition]
        
            # Downsampling
            print('Downsampling...')
            #print(type(data))
            #print(type(time_points))
            data = signal.resample(x=data, num=time_points, axis=1)

            # Morlet transform
            print('Morlet transform...')
            n_channels = data.shape[0]
            n_epochs = data.shape[2]
            morlet_matrix = np.empty((n_frequencies, time_points, n_channels, n_epochs))  
            for epoch in range(n_epochs):
                for channel in range(n_channels):
                    morlet_matrix[:, :, channel, epoch] = wl.Morlet(data[channel, :, epoch], scales=scales).getnormpower()
            
            structures[structure][condition] = np.mean(morlet_matrix, axis=3)
            
    print('Mouse processed in {:.2f} min.\n'.format((time.time() - clock) / 60))

    
    print('Saving dictionary...')    
    pickle.dump(structures, open(npys_dir + 'morlets.epochs', 'wb'), protocol=2)   
    
print('Done!')


# ################

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


# Takes as input $\text{data}_{(n\_channels, time\_points)}$ with downsampled $\text{time points} = 1 KHz \times 6 s = 600000$. Then, for each $data[channel, :] = (600000, 1)$ transforms it to the frequency domain with wl.Morlet, getting a <i>(n_frequencies, time_points)</i> morlet_matrix for each channel. 
