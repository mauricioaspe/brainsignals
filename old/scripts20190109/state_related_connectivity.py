#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pickle
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
import connectivity_measures as cm


# In[2]:


print('Loading files...')
WT = pickle.load(open('/home/maspe/filer/SERT/ALL/npys/band_filtered_WT.dict', 'rb'))
KO = pickle.load(open('/home/maspe/filer/SERT/ALL/npys/band_filtered_KO.dict', 'rb'))


# In[3]:


start = 60000
event = 90000
stop  = 105000

all_WT = dict()
bands = ['theta', 'beta', 'low_gamma', 'alpha']

print('Collecting bands...')
for band in bands:    
    iteration = 0
    
    for mouse in WT.keys():
        
        if iteration == 0:
            data = WT[mouse][band]
            
        else:
            data = np.dstack((data, WT[mouse][band]))
            
        iteration += 1
    
    all_WT[band] = data
    
    
    
all_KO = dict()
bands = ['theta', 'beta', 'low_gamma', 'alpha']

for band in bands:    
    iteration = 0
    
    for mouse in KO.keys():
        
        if iteration == 0:
            data = KO[mouse][band]
            
        else:
            data = np.dstack((data, KO[mouse][band]))
            
        iteration += 1
    
    all_KO[band] = data


# In[5]:


icoh_WT = {'theta': {}, 'alpha': {}, 'beta': {}, 'low_gamma': {}}
icoh_KO = {'theta': {}, 'alpha': {}, 'beta': {}, 'low_gamma': {}}

print('Beginning iCoh for WT...')
for band in all_WT.keys():
    n_mice = all_WT[band].shape[2]
    
    pre_icoh = np.empty((32, 32, n_mice))
    post_icoh = np.empty((32, 32, n_mice))
    
    for mouse in range(n_mice):
        pre_icoh[:, :, mouse]   = cm.icoh(all_WT[band][:, start : event, mouse], average = False)
        post_icoh[:, :, mouse]  = cm.icoh(all_WT[band][:, event : stop, mouse], average = False)
        
    icoh_WT[band]['pre']  = np.mean(pre_icoh, axis=2)
    icoh_WT[band]['post'] = np.mean(post_icoh, axis=2)
    
print('Saving...')
pickle.dump(icoh_WT, open('/home/maspe/filer/SERT/ALL/npys/icoh_WT.dict', 'wb'), protocol=2)

print('Beginning iCoh for KO...')
for band in all_KO.keys():
    n_mice = all_KO[band].shape[2]
    
    pre_icoh = np.empty((32, 32, n_mice))
    post_icoh = np.empty((32, 32, n_mice))
    
    for mouse in range(n_mice):
        pre_icoh[:, :, mouse]   = cm.icoh(all_WT[band][:, start : event, mouse], average = False)
        post_icoh[:, :, mouse]  = cm.icoh(all_WT[band][:, event : stop, mouse], average = False)
        
    icoh_KO[band]['pre']  = np.mean(pre_icoh, axis=2)
    icoh_KO[band]['post'] = np.mean(post_icoh, axis=2)


print('Saving...')
pickle.dump(icoh_KO, open('/home/maspe/filer/SERT/ALL/npys/icoh_KO.dict', 'wb'), protocol=2)
    
print('Done!')


# In[15]:


# For WT
print('Plotting WT...')
for band in icoh_WT.keys():    
    matrix_pre = icoh_WT[band]['pre']
    matrix_post = icoh_WT[band]['post']

    plt.figure(figsize=(10,20))

    plt.matshow(matrix_pre)
    plt.xticks([0+5, 10+3, 16+4, 25+3], ['mPFC', 'NAC', 'BLA', 'vHip'], rotation=0, fontsize=14)
    plt.yticks([0+5, 10+3, 16+4, 25+3], ['mPFC', 'NAC', 'BLA', 'vHip'], rotation=0, fontsize=14)

    plt.axvline(x=10-0.6, color='white')
    plt.axvline(x=16-0.6, color='white')
    plt.axvline(x=25-0.6, color='white')

    plt.axhline(y=10-0.6, color='white')
    plt.axhline(y=16-0.6, color='white')
    plt.axhline(y=25-0.6, color='white')
    plt.colorbar()

    # Saving pre
    print('Saving...')
    plt.savefig('/home/maspe/filer/SERT/ALL/figs/pre_{}_WT.png'.format(band), dpi=150)
    plt.close()
    
    
    # Figure post
    print('Plotting KO...')
    plt.figure(figsize=(10,20))

    plt.matshow(matrix_post)
    plt.xticks([0+5, 10+3, 16+4, 25+3], ['mPFC', 'NAC', 'BLA', 'vHip'], rotation=0, fontsize=14)
    plt.yticks([0+5, 10+3, 16+4, 25+3], ['mPFC', 'NAC', 'BLA', 'vHip'], rotation=0, fontsize=14)

    plt.axvline(x=10-0.6, color='white')
    plt.axvline(x=16-0.6, color='white')
    plt.axvline(x=25-0.6, color='white')

    plt.axhline(y=10-0.6, color='white')
    plt.axhline(y=16-0.6, color='white')
    plt.axhline(y=25-0.6, color='white')
    plt.colorbar()
    
    # Saving post
    print('Saving')
    plt.savefig('/home/maspe/filer/SERT/ALL/figs/post_{}_WT.png'.format(band), dpi=150)
    plt.close()


# In[ ]:


# For KO
for band in icoh_KO.keys():    
    matrix_pre = icoh_KO[band]['pre']
    matrix_post = icoh_KO[band]['post']

    plt.figure(figsize=(10,20))

    plt.matshow(matrix_pre)
    plt.xticks([0+5, 10+3, 16+4, 25+3], ['mPFC', 'NAC', 'BLA', 'vHip'], rotation=0, fontsize=14)
    plt.yticks([0+5, 10+3, 16+4, 25+3], ['mPFC', 'NAC', 'BLA', 'vHip'], rotation=0, fontsize=14)

    plt.axvline(x=10-0.6, color='white')
    plt.axvline(x=16-0.6, color='white')
    plt.axvline(x=25-0.6, color='white')

    plt.axhline(y=10-0.6, color='white')
    plt.axhline(y=16-0.6, color='white')
    plt.axhline(y=25-0.6, color='white')
    plt.colorbar()

    # Saving pre
    plt.savefig('/home/maspe/filer/SERT/ALL/figs/pre_{}_KO.png'.format(band), dpi=150)
    plt.close()
    
    
    # Figure post
    plt.figure(figsize=(10,20))

    plt.matshow(matrix_post)
    plt.xticks([0+5, 10+3, 16+4, 25+3], ['mPFC', 'NAC', 'BLA', 'vHip'], rotation=0, fontsize=14)
    plt.yticks([0+5, 10+3, 16+4, 25+3], ['mPFC', 'NAC', 'BLA', 'vHip'], rotation=0, fontsize=14)

    plt.axvline(x=10-0.6, color='white')
    plt.axvline(x=16-0.6, color='white')
    plt.axvline(x=25-0.6, color='white')

    plt.axhline(y=10-0.6, color='white')
    plt.axhline(y=16-0.6, color='white')
    plt.axhline(y=25-0.6, color='white')
    plt.colorbar()
    
    # Saving post
    plt.savefig('/home/maspe/filer/SERT/ALL/figs/post_{}_KO.png'.format(band), dpi=150)
    plt.close()


# In[ ]:


print('Done!')

