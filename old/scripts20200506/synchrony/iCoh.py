#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pickle
import time
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
import connectivity_measures as cm


# In[2]:


### Contructing filters
### Butterworth filter
def butter_bandpass(lowcut, highcut, fs, order=5):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = signal.butter(order, [low, high], btype='band')
    return b, a

fs      = 1000.0

filter_parameters = {#'delta':    {'N': 3, 'lowcut': 0.5,  'highcut': 3},
                    'theta':    {'N': 4, 'lowcut': 3.5,  'highcut': 7.5},
                    #'alpha':     {'N': 4, 'lowcut': 8.0,  'highcut': 13.0},
                    'beta':      {'N': 4, 'lowcut': 13.0, 'highcut': 25.0},
                    'low_gamma': {'N': 6, 'lowcut': 25.0, 'highcut': 60}}

butterworths = dict()
for band in filter_parameters.keys():
    lowcut = filter_parameters[band]['lowcut']
    highcut = filter_parameters[band]['highcut']
    N = filter_parameters[band]['N']
    butterworths[band] = butter_bandpass(lowcut, highcut, fs, order=N)


# In[4]:


IDs = ['SERT1597', 'SERT1659', 'SERT1678', 'SERT1908', 'SERT1984',
       'SERT1985', 'SERT2014', 'SERT1668', 'SERT1665', 'SERT2018',
       'SERT2024', 'SERT2013']


# ### Probar en ventanas de 500 ms!!

# In[ ]:


condition = 'epochs'

iCohs = {'theta': {}, 'beta': {}, 'low_gamma': {}}
for mouse in IDs: 
    npys_dir = '/home/maspe/filer/SERT/' + mouse + '/npys/'
    print('\nLoading mouse {}...'.format(mouse))
    
    ### Loading data
    data = pickle.load(open(npys_dir + mouse + '.epochs', 'rb'), encoding='latin1')
    
    ### Loop
    filtered_bands = {'theta': {}, 'beta': {}, 'low_gamma': {}}
    iterator = 0
    
    # For this mouse, collect all structures
    # in a (structure_n_channels * n_structures, time_points, n_epochs) matrix
    for structure in ['mPFC', 'NAC', 'BLA', 'vHip']: 
        print('Loading ' + structure + '...')
        if iterator == 0:
            all_structures = data[structure][condition]
        else:
            all_structures = np.vstack((all_structures, data[structure][condition]))
        
        iterator += 1

    # Filters the 32 channels x 180,000 time_points x n epochs matrix along the time_points axis 
    print('Filtering...')    
    for band in filtered_bands.keys():
        filtered = signal.filtfilt(b=butterworths[band][0], a=butterworths[band][1],
                                   x=all_structures, axis=1)


# #### Backup code

# In[6]:


condition = 'epochs'

iCohs = {'theta': {}, 'beta': {}, 'low_gamma': {}}
for mouse in IDs: 
    npys_dir = '/home/maspe/filer/SERT/' + mouse + '/npys/'
    print('\nLoading mouse {}...'.format(mouse))
    
    ### Loading data
    data = pickle.load(open(npys_dir + mouse + '.epochs', 'rb'), encoding='latin1')
    
    ### Loop
    filtered_bands = {'theta': {}, 'beta': {}, 'low_gamma': {}}
    iterator = 0
    
    # For this mouse, collect all structures
    # in a (structure_n_channels * n_structures, time_points, n_epochs) matrix
    for structure in ['mPFC', 'NAC', 'BLA', 'vHip']: 
        print('Loading ' + structure + '...')
        if iterator == 0:
            all_structures = data[structure][condition]
        else:
            all_structures = np.vstack((all_structures, data[structure][condition]))
        
        iterator += 1

    # Filters the 32 channels x 180,000 time_points x n epochs matrix along the time_points axis 
    print('Filtering...')    
    for band in filtered_bands.keys():
        filtered = signal.filtfilt(b=butterworths[band][0], a=butterworths[band][1],
                                   x=all_structures, axis=1)
        
        
        print('Calculating iCoh for {} band...'.format(band))
        clock = time.time()
            
        # Defining ROIs in seconds    
        roi_pre = np.array([-1.5, -1])
        roi_post = np.array([-0.5, 1.5])

        # Converting ROIs to data points
        pre = ((roi_pre + 3) * 30000).astype(int)
        post = ((roi_post + 3) * 30000).astype(int)

        # Loop
        n_epochs = filtered.shape[2]
        for epoch in range(n_epochs):
            if epoch == 0:
                icoh_pre  = cm.icoh(filtered[:,pre[0]:pre[1],epoch], average = False)
                icoh_post = cm.icoh(filtered[:,post[0]:post[1],epoch], average = False)
            else:
                icoh_pre  = np.dstack((icoh_pre, cm.icoh(filtered[:,pre[0]:pre[1],epoch], average = False)))
                icoh_post = np.dstack((icoh_post, cm.icoh(filtered[:,post[0]:post[1],epoch], average = False)))
        
        print('iCoh calculated in {} min'.format(time.time() - clock))
        
        iCohs[band]['pre'] = icoh_pre
        iCohs[band]['post'] = icoh_post

    
    pickle.dump(iCohs, open(npys_dir + mouse + '_BLA_vHip.icoh', 'wb'), protocol=2)

    
print('Done!')  
print('Mouse processed in {} min.'.format((time.time() - clock) / 60))

