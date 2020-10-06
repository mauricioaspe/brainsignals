#!/usr/bin/env python
# coding: utf-8

# In[1]:


import time
import pickle
import numpy as np
import pandas as pd
from scipy import signal


# In[ ]:


IDs_WT = ['SERT1597', 'SERT1659', 'SERT1678', 'SERT1908', 'SERT1984', 'SERT1985', 'SERT2014']
IDs_KO = ['SERT1668', 'SERT2018', 'SERT2024', 'SERT2013'] # 'SERT1665', 

WT = {'mPFC':{}, 'BLA':{}, 'NAC':{}, 'vHip':{}}
KO = {'mPFC':{}, 'BLA':{}, 'NAC':{}, 'vHip':{}}

WT_reshaped = {'mPFC':{}, 'BLA':{}, 'NAC':{}, 'vHip':{}}
KO_reshaped = {'mPFC':{}, 'BLA':{}, 'NAC':{}, 'vHip':{}}

n_samples = 600000
clock = time.time()

print('###################\nLoading WTs...')
for structure in WT.keys():
    for ID in IDs_WT:
        print('Loading {} from {}...'.format(structure, ID))
        npys_dir = '/home/maspe/filer/SERT/' + ID + '/npys/'
        data = np.load(npys_dir + structure + '.npy', allow_pickle=True)  
    
        # Start and end of OF
        files_dir = '/home/maspe/filer/SERT/' + ID + '/continuous/'

        df = pd.read_excel(files_dir + 'xy.xlsx', sheet_name=0, header=None)
        start_OF = int(np.ceil(df[3][0] * 30)) # Times 30 because of sampling at 30 KHz
        stop_OF = int(np.floor(df[3][1] * 30))
    
        print('Downsampling...\n')
        WT[structure][ID] = signal.resample(x=data[:, start_OF : stop_OF], num=n_samples, axis=1)
        
        print('Reshaping...')
        WT_reshaped[structure][ID] = WT[structure][ID].reshape(WT[structure][ID].shape[0], 30000, 20)
        
        
print('Saving dictionaries...')
pickle.dump(WT, open('/home/maspe/filer/SERT/ALL/npys/WT_downsampled.dict', 'wb'), protocol=2)
pickle.dump(WT_reshaped, open('/home/maspe/filer/SERT/ALL/npys/WT_reshaped.dict', 'wb'), protocol=2)
        
        
print('###################\nLoading KOs...')
for structure in KO.keys():
    for ID in IDs_KO:
        print('Loading {} from {}...'.format(structure, ID))
        npys_dir = '/home/maspe/filer/SERT/' + ID + '/npys/'        
        data = np.load(npys_dir + structure + '.npy', allow_pickle=True)
        
        df = pd.read_excel(files_dir + 'xy.xlsx', sheet_name=0, header=None)
        start_OF = int(np.ceil(df[3][0] * 30)) # Times 30 because of sampling at 30 KHz
        stop_OF = int(np.floor(df[3][1] * 30))
    
        print('Downsampling...\n')
        KO[structure][ID] = signal.resample(x=data[:, start_OF : stop_OF], num=n_samples, axis=1)
                
        print('Reshaping...')
        KO_reshaped[structure][ID] = KO[structure][ID].reshape(KO[structure][ID].shape[0], 30000, 20)
        
    
print('\nSaving dictionaries...')
pickle.dump(KO, open('/home/maspe/filer/SERT/ALL/npys/KO_downsampled.dict', 'wb'), protocol=2)
pickle.dump(KO_reshaped, open('/home/maspe/filer/SERT/ALL/npys/KO_reshaped.dict', 'wb'), protocol=2)
        
        
print('All mice downsampled in {:.2f} min.'.format((time.time() - clock) / 60))

