#!/usr/bin/env python
# coding: utf-8

# In[1]:


### Importing modules
import time
import pickle
import numpy as np
import wavelets as wl
from scipy import signal
from matplotlib import pyplot as plt


# In[2]:


# Morlet parameters
dt = 1.0 / 1000
time_windows = np.arange(0, 3, dt)

frequencies = np.arange(1, 50, 1)
#periods = 1 / (frequencies * dt)

#scales = periods / wl.Morlet.fourierwl
#n_frequencies = frequencies.shape[0]
#time_points = time_windows.shape[0]

# Reducing epoch
#start = 1000
#stop = 4000
#baseline = 2000


# #### Grand-average

# In[3]:


print('LOADING!!!')
SRS_WT = pickle.load(open('/home/maspe/filer/SERT/ALL/npys/SRS_WT.dict', 'rb'))
SRS_KO = pickle.load(open('/home/maspe/filer/SERT/ALL/npys/SRS_KO.dict', 'rb'))
print('Done!')


# In[5]:


#n_mice=len(SRS_WT[SRS_WT.keys()[0]].keys())
#n_mice


# In[ ]:


grand_average_WT = {'mPFC': [], 'NAC': [], 'BLA': [], 'vHip': []}
grand_average_KO = {'mPFC': [], 'NAC': [], 'BLA': [], 'vHip': []}

# WT
print('Loading WT...')
n_mice = 7
mouse_average = np.empty((49, 3000, n_mice))
for structure in SRS_WT.keys():    
    iteration = 0
    
    for mouse in SRS_WT[structure].keys():
        print(iteration)
        print('Loading structure {} for mouse {}...'.format(structure, mouse))
        
        mouse_average[:, :, iteration] = np.mean(SRS_WT[structure][mouse], axis=(2,3))
        print('Data shape: {}'.format(mouse_average.shape))
        
        #if iteration == 0:
        #    mouse_average = data
        #    print('Iteration 1 shape: {}'.format(mouse_average.shape))
                   
        #else:
            #print('Original shape average to stack: {}'.format(np.mean(SRS_KO[structure][mouse], axis=(2,3)).shape))
        #    mouse_average = np.stack((mouse_average, data), axis=2)
        
        print('average shape: {}'.format(mouse_average.shape))
        iteration += 1
        
    
    grand_average_WT[structure] = np.mean(mouse_average, axis=2)


print('Print loading KO...')
n_mice = 5
mouse_average = np.empty((49, 3000, n_mice))
for structure in SRS_KO.keys():    
    iteration = 0
  
    for mouse in SRS_KO[structure].keys():
       
        print(iteration)
        print('Loading structure {} for mouse {}...'.format(structure, mouse))
        
        mouse_average[:, :, iteration] = np.mean(SRS_KO[structure][mouse], axis=(2,3))
        print('Data shape: {}'.format(mouse_average.shape))
        
        #if iteration == 0:
        #    mouse_average = data
        #    print('Iteration 1 shape: {}'.format(mouse_average.shape))
                   
        #else:
            #print('Original shape average to stack: {}'.format(np.mean(SRS_KO[structure][mouse], axis=(2,3)).shape))
        #    mouse_average = np.stack((mouse_average, data), axis=2)
        
        print('average shape: {}'.format(mouse_average.shape))
        iteration += 1
    
    grand_average_KO[structure] = np.mean(mouse_average, axis=2)
  

print('Done!')


# ### Plotting

# #### Wild-type

# In[ ]:


# Setting colormap range
mycolormap = {'mPFC': (-3,3), 'NAC': (-8,8), 'BLA': (-4,4), 'vHip': (-4,4)}

for structure in grand_average_WT.keys():
    pwr1 = grand_average_WT[structure]
    fmin = min(frequencies)
    fmax = max(frequencies)

    plt.figure(1, figsize=(10, 5))
    plt.clf()

    ax1 = plt.subplot2grid((1, 5),(0, 0),colspan=4)
       
    plt.imshow(pwr1,cmap='RdBu',vmax=np.max(pwr1),vmin=-np.max(pwr1),
               extent=(min(time_windows),max(time_windows),fmin,fmax),
               origin='lower', interpolation='none',aspect='auto')
    plt.colorbar(fraction=0.05,pad=0.02)
    plt.clim(mycolormap[structure][0], mycolormap[structure][1])
    
    plt.axvline(x=2, color='black')
    plt.axhline(y=3, color='red', linestyle='--')
    plt.axhline(y=8, color='red', linestyle='--')
    plt.axhline(y=13, color='red', linestyle='--')
    plt.axhline(y=25, color='red', linestyle='--')
   
    locs, labels = plt.xticks()    
    plt.xticks(locs, ['-2', '-1.5', '-1', '-0.5', '0', '0.5', '1'], fontsize=14)
    plt.yticks([3, 8, 13, 25], ['3', '8', '13', '25'], fontsize=14)

    ax1.set_xlabel('Time (s)', fontsize=16)
    ax1.set_ylabel('Frequency (Hz)', fontsize=16)
    plt.title('SRS Grand-average for {} in WT'. format(structure), fontsize=20)
    
    plt.savefig('/home/maspe/filer/SERT/ALL/figs/SRS/{}_WT.png'.format(structure), dpi=150, orientation='landscape')


# #### Knock-out

# In[ ]:


for structure in grand_average_KO.keys():
    pwr1 = grand_average_KO[structure]
    fmin = min(frequencies)
    fmax = max(frequencies)

    plt.figure(1, figsize=(10, 5))
    plt.clf()

    ax1 = plt.subplot2grid((1, 5),(0, 0),colspan=4)
       
    plt.imshow(pwr1,cmap='RdBu',vmax=np.max(pwr1),vmin=-np.max(pwr1),
               extent=(min(time_windows),max(time_windows),fmin,fmax),
               origin='lower', interpolation='none',aspect='auto')
    plt.colorbar(fraction=0.05,pad=0.02)
    plt.clim(mycolormap[structure][0], mycolormap[structure][1])
    
    plt.axvline(x=2, color='black')

    plt.axhline(y=3, color='red', linestyle='--')
    plt.axhline(y=8, color='red', linestyle='--')
    plt.axhline(y=13, color='red', linestyle='--')
    plt.axhline(y=25, color='red', linestyle='--')
   
    locs, labels = plt.xticks()    
    plt.xticks(locs, ['-2', '-1.5', '-1', '-0.5', '0', '0.5', '1'], fontsize=14)
    plt.yticks([3, 8, 13, 25], ['3', '8', '13', '25'], fontsize=14)

    ax1.set_xlabel('Time (s)', fontsize=16)
    ax1.set_ylabel('Frequency (Hz)', fontsize=16)
    plt.title('SRS Grand-average for {} in KO'. format(structure), fontsize=20)
    
    plt.savefig('/home/maspe/filer/SERT/ALL/figs/SRS/{}_KO.png'.format(structure), dpi=150, orientation='landscape')

