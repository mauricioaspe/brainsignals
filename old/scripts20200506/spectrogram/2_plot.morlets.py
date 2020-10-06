#!/usr/bin/env python
# coding: utf-8

# # Create morlets from <i>.npys</i> epoched by structure 

# In[ ]:


'''
Input: 
    'SERT/SERTXXXX/npys/morlets_epochs.npy'
    
    created from: timefrequency_epochs.py()
    
    Contains:
    * OF data epoched from 3s before to after each center entrance. 
    * Baseline data epoched from xy.dict of 3 s activity in periphery (at 1 KHz??) and multiplied by 30 to 30 KHz.

    Structured as:
    * Dictionary with baselines and OF epochs sorted by structure.
    * Containing (frequencies, channels, time) matrices with:
        - Channels downsampled to 1 KHz
        - All channels with mean subtracted
        - Low-passed at 300 Hz
        - Median subtracted
        - Time/frequency transformed to 99 frequencies
        
    * Saved at 'SERT/SERTXXXX/npys_dir/morlets_epochs.npy

Output:
    '/home/maspe/filer/SERT/ALL/npys/SRS_KO.dict'
'''


# In[ ]:


### Importing modules
import time
import pickle
import numpy as np
#import wavelets as wl
#from scipy import signal


# ### Parameters

# In[ ]:


### ROI
#start = 1500
#stop = 4500
#baseline = 3000


# ## Main loop

# In[ ]:


IDs = {'SERT1597': {}, 'SERT1659': {}, 'SERT1678': {}, 'SERT1908': {}, 'SERT1984': {}, 'SERT1985': {}, 'SERT2014': {},
       'SERT1665': {}, 'SERT1668': {}, 'SERT2013': {}, 'SERT2018': {}, 'SERT2024': {}} 

structures = {'mPFC': {}, 'BLA': {}, 'NAC': {}, 'vHip': {}}


clock = time.time()
for mouse in IDs.keys():
    print('Loading {}...'.format(mouse))
    
    ### Loading data
    npys_dir = '/home/maspe/filer/SERT/' + mouse + '/npys/'
    data = pickle.load(open(npys_dir + 'morlets.epochs', 'rb'))      
    
    for structure in structures.keys():    
        print('Transforming to z-score...\n')        
        
        ### Getting average and sd on time dimension
        baseline_mean = np.mean(data[structure]['baselines'], axis=(1,3))
        baseline_sd   = np.std(data[structure]['baselines'], axis=(1,3))
        
        IDs[mouse][structure] = ((data['mPFC']['epochs'] - baseline_mean[:, None, :, None])
                                 / baseline_sd[:, None, :, None])
        
        
print('All mice processed in {:.2f} s.'.format(time.time() - clock))
print('Done!')


# In[ ]:


#data['mPFC']['baselines'].shape


# In[ ]:





# In[ ]:




