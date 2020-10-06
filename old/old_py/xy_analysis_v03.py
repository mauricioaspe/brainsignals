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


# In[2]:


final_fs = 30
dt = 1.0 / final_fs

time_windows = np.arange(0, 600, dt) # 600 por numero segundos en los 10 min de OF
time_points = time_windows.shape[0]


# ###############

# In[3]:


filt_params = {'theta'    : (3,  8),
               'alpha'    : (8,  13),
               'beta'     : (13, 25),
               'low_gamma': (25, 60)}


# In[5]:


IDs_WT = {'SERT1678': {}, 'SERT1597': {}, 'SERT1659': {}} #, 'SERT1908', 'SERT1984', 'SERT1985', 'SERT2014'] #'SERT1678',
IDs_KO = {'SERT1668': {}, 'SERT1665': {}} #, 'SERT2018', 'SERT2024', 'SERT2013'] 

WT = {'mPFC': {'theta': {}, 'alpha': {}, 'beta': {}, 'low_gamma': {}},
      'NAC': {'theta': {}, 'alpha': {}, 'beta': {}, 'low_gamma': {}},
      'BLA': {'theta': {}, 'alpha': {}, 'beta': {}, 'low_gamma': {}},
      'vHip': {'theta': {}, 'alpha': {}, 'beta': {}, 'low_gamma': {}}}

#KO = {'mPFC': {'theta': {}, 'alpha': {}, 'beta': {}, 'low_gamma': {}},
#      'NAC': {'theta': {}, 'alpha': {}, 'beta': {}, 'low_gamma': {}},
#      'BLA': {'theta': {}, 'alpha': {}, 'beta': {}, 'low_gamma': {}},
#      'vHip': {'theta': {}, 'alpha': {}, 'beta': {}, 'low_gamma': {}}}

clock = time.time()
xy = dict()
for band in filt_params.keys():
    iteration = 0
    print(iteration)

    
    print('###################\nLoading continuous...')
    for structure in WT.keys():    
        for mouse in IDs_WT.keys():
            print('Loading {} from {}...'.format(structure, mouse))
            npys_dir = '/home/maspe/filer/SERT/' + mouse
            
            if iteration == 0:
                # Start and end of OF
                df = pd.read_excel(npys_dir + '/continuous/xy.xlsx', sheet_name=0, header=None)
                xy[mouse] = (df[1].to_numpy(), df[2].to_numpy())
                OF_npoints = int(np.floor(df[3][1])) - int(np.ceil(df[3][0]))
                print('OF points: {}'.format(OF_npoints))
                print('xy length: {}'.format(df[1].to_numpy().shape))

            data = np.load(npys_dir + '/npys/' + structure + '_morlet.npy', allow_pickle=True)
            print('original data shape: {}'.format(data.shape))
        
            min_freq = filt_params[band][0]
            max_freq = filt_params[band][1]
        
            # Modificacion
            print('Collecting {} band...'.format(band))
            data = np.mean(data[min_freq : max_freq, :, :], axis=(2, 0))
            
            print('Resampling...')
            WT[structure][band][mouse] = signal.resample(x=data, num=time_points)

            
        iteration += 1


print('Saving dictionaries...')
pickle.dump(xy, open('/home/maspe/filer/SERT/ALL/npys/xy_WT.dict', 'wb'), protocol=2)
pickle.dump(WT, open('/home/maspe/filer/SERT/ALL/npys/bands_xy_WT.dict', 'wb'), protocol=2)

print('Done!')     


# In[3]:


#WT = pickle.load(open('/home/maspe/filer/SERT/ALL/npys/bands_xy_WT.dict', 'rb'))
#xy = pickle.load(open('/home/maspe/filer/SERT/ALL/npys/xy_WT.dict', 'rb'))


# In[20]:


#band = 'alpha'
#x = xy['SERT1678'][0]
#y = xy['SERT1678'][1]
#spectrum = WT['vHip'][band]['SERT1597']
#norm_spectrum = (spectrum - np.min(spectrum)) / (np.max(spectrum) - np.min(spectrum))
#plt.scatter(x[1:], y[1:], s=16, c=norm_spectrum, cmap='jet', alpha=0.1)


# In[ ]:


#plt.scatter(x[1:], y[1:], s=16, c=norm_test, cmap='jet', alpha=0.3)


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

