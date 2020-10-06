#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import pickle
import time
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
import connectivity_measures as cm


# In[ ]:


IDs = ['SERT1597', 'SERT1659', 'SERT1678', 'SERT1908',
       'SERT1984', 'SERT1985', 'SERT2014', 'SERT1668',
       'SERT1665', 'SERT2018', 'SERT2024', 'SERT2013']


# In[ ]:


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


# In[4]:


condition = 'epochs'

wPLI = {'theta': {}, 'low_gamma': {}}
for mouse in IDs: 
    npys_dir = '/home/maspe/filer/SERT/' + mouse + '/npys/'
    print('\nLoading mouse {}...'.format(mouse))
    
    ### Loading data
    data = pickle.load(open(npys_dir + mouse + '.epochs', 'rb') , encoding='latin1')
       
    ### Loop
    filtered_bands = {'theta': {}, 'low_gamma': {}}
    iterator = 0
    
    print('Stacking structures...')
    for structure in data.keys(): 
        if iterator == 0:
            all_structures = data[structure][condition]
        else:
            all_structures = np.vstack((all_structures, data[structure][condition]))
        
        iterator += 1
        
    print('Filtering...')    
    for band in filtered_bands.keys():
        filtered = signal.filtfilt(b=butterworths[band][0], a=butterworths[band][1],
                                   x=all_structures, axis=1)
        
        print('Getting phases for {} band...'.format(band))
        transformed = signal.hilbert(filtered)
        phases = np.angle(transformed)
        
        print('Calculating wPLI for {} band...'.format(band))
        clock = time.time()
        
        roi_pre = np.array([-2, 0])
        roi_post = np.array([2, 2.5])

        pre = ((roi_pre + 3) * 30000).astype(int)
        post = ((roi_post + 3) * 30000).astype(int)

        
        n_epochs = phases.shape[2]
        for epoch in range(n_epochs):
            if epoch == 0:
                wpli_pre  = cm.PLI2(phases[:, pre[0]:pre[1], epoch], average = False)
                wpli_post = cm.PLI2(phases[:, post[0]:post[1], epoch], average = False)
            else:
                wpli_pre  = np.dstack((wpli_pre, cm.PLI2(phases[:,pre[0]:pre[1],epoch], average = False)))
                wpli_post = np.dstack((wpli_post, cm.PLI2(phases[:,post[0]:post[1],epoch], average = False)))
        
        print('wPLI calculated in {} s.'.format(time.time() - clock))
        
        wPLI[band]['pre'] = wpli_pre
        wPLI[band]['post'] = wpli_post
      
    
    pickle.dump(wPLI, open(npys_dir + mouse + '.pli', 'wb'), protocol=2)


print('Done!')        

