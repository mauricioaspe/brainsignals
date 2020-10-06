#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pickle
import numpy as np
from scipy import signal
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

filter_parameters = {'theta':    {'N': 4, 'lowcut': 3.5,  'highcut': 7.5},
                    'alpha':     {'N': 4, 'lowcut': 8.0,  'highcut': 13.0},
                    'beta':      {'N': 4, 'lowcut': 13.0, 'highcut': 25.0},
                    'low_gamma': {'N': 6, 'lowcut': 25.0, 'highcut': 60}}

butterworths = dict()
for band in filter_parameters.keys():
    lowcut = filter_parameters[band]['lowcut']
    highcut = filter_parameters[band]['highcut']
    N = filter_parameters[band]['N']
    butterworths[band] = butter_bandpass(lowcut, highcut, fs, order=N)


# In[ ]:


WT = {'SERT1597': {}, 'SERT1659': {}} #, 'SERT1678': {}, 'SERT1908': {}, 'SERT1984': {}, 'SERT1985': {}, 'SERT2014': {}}
KO = {'SERT1668': {}, 'SERT1665': {}} #, 'SERT2018': {}, 'SERT2024': {}, 'SERT2013': {}}


print('\n#############\nProcessing WT')
filtered = {'theta': [], 'alpha': [], 'beta': [], 'low_gamma': []}
structures = ['mPFC', 'NAC', 'BLA', 'vHip']

for ID in WT.keys():
    npys_dir = '/home/maspe/filer/SERT/' + ID + '/npys/'
    
    print('\nLoading {}...'.format(npys_dir))
       
    for band in filtered.keys():
        iteration = 0
        print('outside iteration: {}'.format(iteration))
        for structure in structures:
            print(iteration)
            data = np.load(npys_dir + structure + '_epochs.npy', allow_pickle=True)
            print('loaded data shape: {}'.format(data.shape))
            print('Filtering and averaging {} band in {}'.format(band, structure))
           
            data = signal.filtfilt(b=butterworths[band][0],
                                   a=butterworths[band][1],
                                   x=data, axis=1)
            
            print('filtered data shape: {}'.format(data.shape))
            if iteration == 0:
                print('Entering if!!')
                mean_data = np.mean(data, axis=2) # Collapses epochs
                print('initial mean data shape: {}'.format(mean_data.shape))
                    
            else:
                print('previous mean shape: {}'.format(np.mean(data, axis=2).shape))
                mean_data = np.vstack((mean_data, np.mean(data, axis=2)))
                print('stacked mean data shape: {}'.format(mean_data.shape))
                
            print('mean data shape: {}'.format(mean_data.shape))
        
            iteration += 1
        
        WT[ID][band] = mean_data

    
print('\nSaving filtered dictionary...')
pickle.dump(WT, open('/home/maspe/filer/SERT/ALL/npys/band_filtered_WT_2.dict', 'wb'), protocol=2)


print('\n#############\nProcessing KO')
filtered = {'theta': [], 'alpha': [], 'beta': [], 'low_gamma': []}
structures = ['mPFC', 'NAC', 'BLA', 'vHip']

for ID in KO.keys():
    npys_dir = '/home/maspe/filer/SERT/' + ID + '/npys/'
    
    print('\nLoading {}...'.format(npys_dir))
       
    for band in filtered.keys():
        iteration = 0
        print('outside iteration: {}'.format(iteration))
        for structure in structures:
            print(iteration)
            data = np.load(npys_dir + structure + '_epochs.npy', allow_pickle=True)
            print('loaded data shape: {}'.format(data.shape))
            print('Filtering and averaging {} band in {}'.format(band, structure))
           
            data = signal.filtfilt(b=butterworths[band][0],
                                   a=butterworths[band][1],
                                   x=data, axis=1)
            
            print('filtered data shape: {}'.format(data.shape))
            if iteration == 0:
                print('Entering if!!')
                mean_data = np.mean(data, axis=2) # Collapses epochs
                print('initial mean data shape: {}'.format(mean_data.shape))
                    
            else:
                print('previous mean shape: {}'.format(np.mean(data, axis=2).shape))
                mean_data = np.vstack((mean_data, np.mean(data, axis=2)))
                print('stacked mean data shape: {}'.format(mean_data.shape))
                
            print('mean data shape: {}'.format(mean_data.shape))
        
            iteration += 1
        
        KO[ID][band] = mean_data

    
print('\nSaving filtered dictionary...')
pickle.dump(KO, open('/home/maspe/filer/SERT/ALL/npys/band_filtered_KO_2.dict', 'wb'), protocol=2)

print('\nDone!')

