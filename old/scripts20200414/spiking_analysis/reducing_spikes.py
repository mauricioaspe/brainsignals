#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import pickle
import time
import numpy as np
from matplotlib import pyplot as plt


# #### Loading spikes, reduce it to 1 Hz and save the reduced dictionary

# In[ ]:


npys_dir = '/home/maspe/filer/SERT/ALL/npys/'
IDs_WT = ['SERT1597', 'SERT1659', 'SERT1678', 'SERT1908', 'SERT1984', 'SERT1985', 'SERT2014']
IDs_KO = ['SERT1668', 'SERT1665', 'SERT2018', 'SERT2024', 'SERT2013'] 


# In[ ]:


sample_rate = 30000

### For WTs
print('Processing WTs...')
reduced_spikes_WT = {}
for ID in IDs_WT:
    print('Importing {}...'.format(ID))
    clock = time.time()
    
    # Get the raw data, sampled at 30 kHz
    data = np.load(npys_dir + ID + '_spikes_WT.npy', allow_pickle=True)
    
    # Reduce it to a sample rate of 1 Hz
    reduced_spikes_WT[ID] = np.add.reduceat(data, range(0, data.shape[1], sample_rate), axis=1)
        
    print('{} imported in {} min.'.format(ID, (time.time() - clock) / 60))

    
print('Saving the data for WTs...')
with open(npys_dir + 'reduced_WT.spikes', 'wb') as f:
    pickle.dump(reduced_spikes_WT, f, protocol=2)
print('Done!')


### For KOs
print('Processing KOs...')
reduced_spikes_KO = {}
for ID in IDs_KO:
    print('Importing {}...'.format(ID))
    clock = time.time()
    
    # Get the raw data, sampled at 30 kHz
    data = np.load(npys_dir + ID + '_spikes_KO.npy', allow_pickle=True)
    
    # Reduce it to a sample rate of 1 Hz
    reduced_spikes_KO[ID] = np.add.reduceat(data, range(0, data.shape[1], sample_rate), axis=1)
        
    print('{} imported in {} min.'.format(ID, (time.time() - clock) / 60))

    
print('Saving the data for KOs...')
with open(npys_dir + 'reduced_KO.spikes', 'wb') as f:
    pickle.dump(reduced_spikes_KO, f, protocol=2)
print('Done!')

