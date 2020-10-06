#!/usr/bin/env python
# coding: utf-8

# In[7]:


import h5py
import pickle
import time
import numpy as np
from datetime import datetime
from datetime import date

#import pandas as pd
#import wavelets as wl
#from scipy import signal
#from matplotlib import pyplot as plt


# #### Get a dictionary with mice's experiment information and another with the spiking activity in h5py format  

# In[2]:


# List of mice, classified by condition
IDs_WT = ['SERT1597', 'SERT1659', 'SERT1678', 'SERT1908', 'SERT1984', 'SERT1985', 'SERT2014']
IDs_KO = ['SERT1668', 'SERT1665', 'SERT2013', 'SERT2018', 'SERT2024']

# Dictionaries to collect info about the experiment and spiking activity
all_units_WT = {}
all_channels_WT = {}
all_units_KO = {}
all_channels_KO = {}

all_info_WT = {}
all_mice_WT = {}
all_info_KO = {}
all_mice_KO = {}

n_channels = 32
    
print('Processing wild-types...')
for ID in IDs_WT:
    npys_dir = '/home/maspe/filer/SERT/' + ID + '/npys/'
    spikes_dir = '/home/maspe/filer/SERT/' + ID + '/spikes/results/'
    
    print('Loading ' + ID)
    with open(npys_dir + ID + '.info', 'rb') as f:
        info = pickle.load(f, encoding='latin1')
    
    # Put the info of this mouse into de info dictionary
    all_info_WT[ID] = info
    channels = info['channels_list']
    
    fit = []
    for channel in channels:
        path = spikes_dir + channel + '.result.hdf5'
    
        # Get the spiking activity from the h5py file
        fit.append(h5py.File(path, 'r'))
        
     
    # Put the spiking activity of this mouse into the spiking dictionary
    all_mice_WT[ID] = fit
    
    
    ########################
    units = []
    channels_id = []
    
    iteration = 0
    for channel in range(n_channels):
        for unit in all_mice_WT[ID][channel]['spiketimes'].keys():
            units.append(all_mice_WT[ID][channel]['spiketimes'][unit][()]) # Final "[()]" is to import values from h5py 
      
            channels_id.append(all_info_WT[ID]['channels_locs'][iteration])
        
        iteration += 1
        
            
    all_units_WT[ID] = units
    all_channels_WT[ID] = channels_id
    
    
### Same for KO mice
#print('Processing knock-out')
print('Processing knock-outs...')
for ID in IDs_KO:
    npys_dir = '/home/maspe/filer/SERT/' + ID + '/npys/'
    spikes_dir = '/home/maspe/filer/SERT/' + ID + '/spikes/results/'
    
    print('Loading ' + ID)
    with open(npys_dir + ID + '.info', 'rb') as f:
        info = pickle.load(f, encoding='latin1')

    all_info_KO[ID] = info
    channels = info['channels_list']
    
    fit = []
    for channel in channels:
        path = spikes_dir + channel + '.result.hdf5'
    
        fit.append(h5py.File(path, 'r'))
                
    all_mice_KO[ID] = fit
    
    
    ########################
    units = []
    channels_id = []
    
    iteration = 0
    for channel in range(n_channels):
        for unit in all_mice_KO[ID][channel]['spiketimes'].keys():
            units.append(all_mice_KO[ID][channel]['spiketimes'][unit][()]) # Final "[()]" is to import values from h5py 
      
            channels_id.append(all_info_KO[ID]['channels_locs'][iteration])
        
        iteration += 1
        
            
    all_units_KO[ID] = units
    all_channels_KO[ID] = channels_id
    #########################

# Saving channels info for spiking analysis
print('Saving channels info for spiking analysis...')

with open('/home/maspe/filer/SERT/ALL/npys/channels_by_spikes_WT.info', 'wb') as f:
    pickle.dump(all_channels_WT, f, protocol=2)
with open('/home/maspe/filer/SERT/ALL/npys/channels_by_spikes_KO.info', 'wb') as f:
    pickle.dump(all_channels_KO, f, protocol=2)



print('Done!')
    
# all_mice_WT['SERT1597'][0]['spiketimes']['temp_0'][()].shape
# all_units_WT['SERT1597']


# #### Create the windows of interest

# We are currently using a time windows of 2 minutes that finishes 15 seconds before the beginning of the OF, for the HC condition; and a time windows that begins 15 seconds after the beginning of the OF, and extends to minute 10 (total length: 9.75 min).

# In[ ]:


### THIS CAN BE MOVED TO THE CREATION OF INFO FILES!!! ###
### Extract a windows of n_secs length at the beginning and end of the task
n_min = 2
exclude_sec = 15
sample_rate = 30000
#window = int(sampleRate * secs)
n_points = 30000 * 60 * 10 # Length of the OF
window = sample_rate * n_min * 60 # 2 min windows
exclude_window = sample_rate * exclude_sec # 15 sec windows to exclude at the HC-OF transition

all_epochs_WT = {}
all_perispikes_WT = {}
for mouse in all_mice_WT.keys():
    stopHC = np.int(all_info_WT[mouse]['startOF']) - exclude_window
    startHC = stopHC - window
    
    startOF = np.int(all_info_WT[mouse]['startOF']) + exclude_window
    stopOF = startOF + n_points - exclude_window
    
    baseline = np.arange(startHC, stopHC, 1)
    task = np.arange(startOF, stopOF, 1)
    
    complete_window = np.concatenate([baseline, task])
    all_epochs_WT[mouse] = complete_window
    

all_epochs_KO = {}
all_perispikes_KO = {}
for mouse in all_mice_KO.keys():
    stopHC = np.int(all_info_KO[mouse]['startOF']) - exclude_window
    startHC = stopHC - window
    
    startOF = np.int(all_info_KO[mouse]['startOF']) + exclude_window
    stopOF = startOF + n_points - exclude_window
    
    baseline = np.arange(startHC, stopHC, 1)
    task = np.arange(startOF, stopOF, 1)
    
    complete_window = np.concatenate([baseline, task])
    all_epochs_KO[mouse] = complete_window


task_npoints = complete_window.shape[0]
print('Windows of interest created!')


# #### Extracting the spikes that fall inside the time windows

# In[ ]:


### For WT
all_perispikes_WT = {}
npys_dir = '/home/maspe/filer/SERT/ALL/npys/'

print('Loading WTs')
clock = time.time()
for mouse in all_mice_WT.keys():
    all_spikes = all_units_WT[mouse]
    peristimulus_spikes = np.zeros((len(all_spikes), task_npoints))   
    
    print('Processing mouse %s...' % mouse)
    
    for unit in range(len(all_spikes)):
        peristimulus_spikes[unit, :] = np.isin(all_epochs_WT[mouse], all_spikes[unit])
        
        
    all_perispikes_WT[mouse] = peristimulus_spikes
    np.save(npys_dir + mouse + '_spikes_WT.npy', peristimulus_spikes)
    
print('WTs done in {} min.!'.format((time.time() - clock) / 60))

### For KO
all_perispikes_KO = {}
      
print('Loading KOs')
clock = time.time()
for mouse in all_mice_KO.keys():
    all_spikes = all_units_KO[mouse] 
    peristimulus_spikes = np.zeros((len(all_spikes), task_npoints))   
    
    print('Processing mouse %s...' % mouse)

    for unit in range(len(all_spikes)):
        peristimulus_spikes[unit, :] = np.isin(all_epochs_KO[mouse], all_spikes[unit])

        
    all_perispikes_KO[mouse] = peristimulus_spikes
    np.save(npys_dir + mouse + '_spikes_KO.npy', peristimulus_spikes)
    
print('KO done in {} min.!'.format((time.time() - clock) / 60)) 


# In[9]:


now = datetime.now()

current_time = now.strftime("%H:%M:%S")
print("Script finished on {} at {}".format(date.today(), current_time))


# In[ ]:


#np.sum(all_perispikes_KO['SERT1668'])

