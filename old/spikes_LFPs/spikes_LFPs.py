#!/usr/bin/env python
# coding: utf-8

# In[2]:


import pickle
import time
import numpy as np
from scipy import signal
from matplotlib import pyplot as plt


# In[3]:


spikes_epochs = pickle.load(open('/home/maspe/filer/SERT/ALL/npys/spikes2.epochs', 'rb'), encoding='latin1')


# In[6]:


IDs = pickle.load(open('/home/maspe/filer/scripts/preprocessing/IDs.dict', 'rb'))['list_all']
WTs = pickle.load(open('/home/maspe/filer/scripts/preprocessing/IDs.dict', 'rb'))['list_WT']


# In[ ]:


ms = 200
window_length = ms * 30
window = [0, window_length] # Windows in ms
windows = []

for i in range(30):
    windows.append(window)
    window = [window[0] + window_length, window[1] + window_length]


# In[ ]:


mPFC_units = np.arange(0,10)
vHip_channels = np.arange(26,32)
phases, amplitudes = dict(), dict()

clock = time.time()
for mouse in IDs:
    print('Loading mouse {}...'.format(mouse))
    
    npys_dir = '/home/maspe/filer/SERT/' + mouse + '/npys/' 
    LFPs_epochs = pickle.load(open(npys_dir + mouse + '.epochs_filtered', 'rb'), encoding='latin1')
    
    spikes_LFPs = {'theta': [], 'beta': [], 'low_gamma': []}
    #spikes_LFPs_KOs = {'theta': [], 'beta': [], 'low_gamma': []}
    #print('Lenght 1 theta: {}'.format(len(spikes_LFPs['theta'])))

    
    for band in spikes_LFPs.keys():
        print('Shape LFPs_epochs: {}'.format(LFPs_epochs[band].shape))
        print('Applying Hilbert transform to {} band...'.format(band))
        print('How many windows? {}'.format(len(windows)))

        z = signal.hilbert(LFPs_epochs[band])

        phases[band]     = np.angle(z)
        amplitudes[band] = np.absolute(z)
   
     
        print("Collecting spikes...")
        all_spikes = spikes_epochs[mouse]
        n_epochs = phases[band].shape[2]
        print('Number of epochs: {}'.format(n_epochs))

        for window in windows:
            for unit in mPFC_units:
                this_unit = []
                for channel in vHip_channels: # vHip channels
                    for epoch in range(n_epochs):
                        spikes = all_spikes[unit][epoch]
                        windows_spikes = spikes[(spikes > window[0]) & (spikes < window[1])]
                        this_unit.extend(phases[band][channel,:,epoch][windows_spikes])
                        #print('Unit length: {}'.format(len(this_unit)))
    
            #print('Tuple to append: ')
            spikes_LFPs[band].append((np.mean(this_unit), np.std(this_unit)))
            #print('Lenght 2 theta: {}'.format(len(spikes_LFPs['theta'])))

    print('Saving...')
    pickle.dump(spikes_LFPs, open(npys_dir + mouse + '.spikes_LFPs', 'wb'), protocol = 2)
    
    print('Mouse processed in {} min!'.format((time.time() - clock) / 60))


print('Done!')


# ### Backup !!!

# In[ ]:


#collapsed_units = []
#
#    for unit in range(n_units):
#        collapsed_units.append([item for sublist in spikes_epochs[mouse][unit] for item in sublist])

