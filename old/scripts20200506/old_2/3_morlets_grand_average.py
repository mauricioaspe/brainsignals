#!/usr/bin/env python
# coding: utf-8

# # Plot morlets Grand-averages

# In[ ]:


'''
Input:
    '/home/maspe/filer/SERT/ALL/npys/SRS_KO.dict'
    
Output:
    '/home/maspe/filer/SERT/ALL/figs/SRS/mPFC_WT.png'

'''


# In[ ]:


### Importing modules
import time
import pickle
import numpy as np
import wavelets as wl
from scipy import signal
from matplotlib import pyplot as plt


# In[ ]:


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
frequencies = np.arange(1, 100, 1)
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
print('Loading data...')
clock = time.time()
SRS_WT = pickle.load(open('/home/maspe/filer/SERT/ALL/npys/SRS_WT_nobaseline.dict', 'rb'))
print('WT loaded!')
SRS_KO = pickle.load(open('/home/maspe/filer/SERT/ALL/npys/SRS_KO_nobaseline.dict', 'rb'))
print('KO loaded!')
print('Data loaded in {} min'.format((time.time() - clock) / 60))


# #### Grand-average

# In[ ]:


grand_average_WT = {'mPFC': [], 'NAC': [], 'BLA': [], 'vHip': []}
grand_average_KO = {'mPFC': [], 'NAC': [], 'BLA': [], 'vHip': []}

# WT
print('Collecting Morlets for WT...')
for structure in SRS_WT.keys():    
    iteration = 0
  
    n_mice = len(SRS_WT[structure].keys())
    mouse_average = np.empty((49, 3000, n_mice))
    for mouse in SRS_WT[structure].keys():
        print('numpy mean mouse {} structure {} shape: {}'.format(mouse, structure, np.mean(SRS_WT[structure][mouse], axis=(2,3)).shape))
        
        mouse_average[:, :, iteration] = np.mean(SRS_WT[structure][mouse], axis=(2,3))

        print('mouse average shape: {}'.format(mouse_average.shape))
        iteration += 1
        
    
    grand_average_WT[structure] = np.mean(mouse_average, axis=2)


print('Collecting Morlets for KO...')
for structure in SRS_KO.keys():    
    iteration = 0
  
    n_mice = len(SRS_KO[structure].keys())
    mouse_average = np.empty((49, 3000, n_mice))
    for mouse in SRS_KO[structure].keys():
        #print(iteration)
        #print(mouse)
        print('numpy mean mouse {} structure {} shape: {}'.format(mouse, structure, np.mean(SRS_KO[structure][mouse], axis=(2,3)).shape))
        
        mouse_average[:, :, iteration] = np.mean(SRS_KO[structure][mouse], axis=(2,3))
            
        #else:
        #    mouse_average = np.stack((mouse_average, np.mean(SRS_KO[structure][mouse], axis=(2,3))), axis=2)
        
        print('mouse average shape: {}'.format(mouse_average.shape))

        iteration += 1
        
    
    grand_average_KO[structure] = np.mean(mouse_average, axis=2)

    
print('Done!')


# ### Plotting

# #### Wild-type

# In[ ]:


print('Plotting...')
for structure in grand_average_WT.keys():
    print(grand_average_WT[structure])
    pwr1 = grand_average_WT[structure]
    fmin = min(frequencies)
    fmax = max(frequencies)

    plt.figure(1, figsize=(10, 4))
    plt.clf()

    ax1 = plt.subplot2grid((1, 5),(0, 0),colspan=4)
    
    
       
    plt.imshow(pwr1,cmap='RdBu',vmax=np.max(pwr1),vmin=-np.max(pwr1),
               extent=(min(time_windows),max(time_windows),fmin,fmax),
               origin='lower', interpolation='none',aspect='auto')
    plt.colorbar(fraction=0.05,pad=0.02)
    
    plt.axvline(x=2, color='black')
   
    locs, labels = plt.xticks()    
    plt.xticks(locs, ['-2', '-1.5', '-1', '-0.5', '0', '0.5', '1'], fontsize=14)
    plt.yticks(fontsize=14)

    #ax1.set_yscale('log')
    ax1.set_xlabel('Time (s)', fontsize=16)
    ax1.set_ylabel('Frequency (Hz)', fontsize=16)
    plt.title('SRS Grand-average for {} in WT'. format(structure), fontsize=20)
    
    print('Saving figures for WT...')
    plt.savefig('/home/maspe/filer/SERT/ALL/figs/SRS/{}_WT_nobaseline.png'.format(structure), dpi=150, orientation='landscape')


# #### Knock-out

# In[ ]:


for structure in grand_average_KO.keys():
    print(grand_average_WT[structure])
    pwr1 = grand_average_KO[structure]
    fmin = min(frequencies)
    fmax = max(frequencies)

    plt.figure(1, figsize=(10, 4))
    plt.clf()

    ax1 = plt.subplot2grid((1, 5),(0, 0),colspan=4)
    
    
       
    plt.imshow(pwr1,cmap='RdBu',vmax=np.max(pwr1),vmin=-np.max(pwr1),
               extent=(min(time_windows),max(time_windows),fmin,fmax),
               origin='lower', interpolation='none',aspect='auto')
    plt.colorbar(fraction=0.05,pad=0.02)
    
    plt.axvline(x=2, color='black')
   
    locs, labels = plt.xticks()    
    plt.xticks(locs, ['-2', '-1.5', '-1', '-0.5', '0', '0.5', '1'], fontsize=14)
    plt.yticks(fontsize=14)

    #ax1.set_yscale('log')
    ax1.set_xlabel('Time (s)', fontsize=16)
    ax1.set_ylabel('Frequency (Hz)', fontsize=16)
    plt.title('SRS Grand-average for {} in KO'. format(structure), fontsize=20)
    
    print('Saving figures for KO...')
    plt.savefig('/home/maspe/filer/SERT/ALL/figs/SRS/{}_KO_nobaseline.png'.format(structure), dpi=150, orientation='landscape')

