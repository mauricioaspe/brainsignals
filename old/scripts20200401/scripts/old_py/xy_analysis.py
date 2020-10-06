#!/usr/bin/env python
# coding: utf-8

# In[2]:


### Importing modules
import time
import pickle
import numpy as np
import pandas as pd
from matplotlib import cm
from matplotlib import pyplot as plt
from scipy import signal


# ###############

# In[ ]:


filt_params = {'theta'    : (3,  8),
               'alpha'    : (8,  13),
               'beta'     : (13, 25),
               'low_gamma': (25, 60)}


# In[ ]:


IDs_WT = {'SERT1678': {}, 'SERT1597': {}, 'SERT1659': {}} #, 'SERT1908', 'SERT1984', 'SERT1985', 'SERT2014'] #'SERT1678',
IDs_KO = {'SERT1668': {}, 'SERT1665': {}} #, 'SERT2018', 'SERT2024', 'SERT2013'] 

WT = {'mPFC':{}, 'BLA':{}, 'NAC':{}, 'vHip':{}}
KO = {'mPFC':{}, 'BLA':{}} #, 'NAC':{}, 'vHip':{}}
clock = time.time()


band = 'theta'
print('###################\nLoading continuous...')
for structure in WT.keys():
    iteration = 0
    for mouse in IDs_WT.keys():
        print(iteration)
        print('Loading {} from {}...'.format(structure, mouse))
        npys_dir = '/home/maspe/filer/SERT/' + mouse + '/npys/'
        data = np.load(npys_dir + structure + '_morlet.npy', allow_pickle=True)
        print('original data shape: {}'.format(data.shape))
        
        min_freq = filt_params[band][0]
        max_freq = filt_params[band][1]
        
        print('Collecting {} band...'.format(band))
        if iteration == 0:
            matrix = np.mean(data[min_freq : max_freq, :, :], axis=(2, 0))
            print(matrix.shape)
        else:
            temp = np.mean(data[min_freq : max_freq, :, :], axis=(2, 0))
            matrix = np.vstack((matrix, temp))
            print(matrix.shape)
        
        print('final data shape: {}'.format(matrix.shape))
        
        iteration += 1
        
        WT[structure][band] = matrix    


print('Saving dictionary...')
pickle.dump(WT, open('/home/maspe/filer/SERT/ALL/npys/bands_xy_WT.dict', 'wb'), protocol=2)

print('Done!')     


# In[3]:


WT = pickle.load(open('/home/maspe/filer/SERT/ALL/npys/bands_xy_WT.dict', 'rb'))


# In[5]:


WT['BLA']['theta'].shape


# In[ ]:


#df = pd.read_csv('/home/maspe/filer/xy/SERT1597OF_WT.csv', header=None)
#start = int(df[7][0])
#stop = int(df[7][1])


# In[ ]:


#x = df[1].to_numpy()
#y = df[2].to_numpy()

#plt.plot(x,y)


# In[ ]:


#norm_test = (test - np.min(test)) / (np.max(test) - np.min(test))


# In[ ]:


#plt.scatter(x[1:], y[1:], s=16, c=norm_test, cmap='jet', alpha=0.3)

