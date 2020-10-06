#!/usr/bin/env python
# coding: utf-8

# In[1]:


### Importing modules
import time
import pickle
import numpy as np
from matplotlib import cm
from matplotlib import pyplot as plt


# In[ ]:


master_clock = time.time()
IDs_WT = ['SERT1597', 'SERT1659', 'SERT1678', 'SERT1908', 'SERT1984', 'SERT1985', 'SERT2014']
IDs_KO = ['SERT1668', 'SERT1665', 'SERT2018', 'SERT2024', 'SERT2013'] 

WT = {'mPFC':{}, 'BLA':{}, 'NAC':{}, 'vHip':{}}
KO = {'mPFC':{}, 'BLA':{}, 'NAC':{}, 'vHip':{}}
clock = time.time()

print('###################\nLoading WTs...')
for structure in WT.keys():
    for ID in IDs_WT:
        print('Loading {} from {}...'.format(structure, ID))
        npys_dir = '/home/maspe/filer/SERT/' + ID + '/npys/'
        WT[structure][ID] = np.load(npys_dir + structure + '_morlet.npy', allow_pickle=True)   
        
        
print('###################\nLoading KOs...')
for structure in KO.keys():
    for ID in IDs_KO:
        print('Loading {} from {}...'.format(structure, ID))
        npys_dir = '/home/maspe/filer/SERT/' + ID + '/npys/'        
        KO[structure][ID] = np.load(npys_dir + structure + '_morlet.npy', allow_pickle=True)
    
        
print('All mice loaded in {:.2f} s.'.format(time.time() - clock))


# In[ ]:


print('Collecting Morlets...')
theta = [5, 12]
#lgamma = [50, 80]
     
theta_morlets_WT = {'mPFC':{}, 'BLA':{}, 'NAC':{}, 'vHip':{}}
theta_morlets_KO = {'mPFC':{}, 'BLA':{}, 'NAC':{}, 'vHip':{}}

# For WT
for structure in theta_morlets_WT.keys():
    for mouse in WT[structure].keys():
        theta_morlets_WT[structure][mouse] = np.mean(WT[structure][mouse][theta[0]:theta[1], :, :], axis=(0, 2))

# For KO
for structure in theta_morlets_KO.keys():
    for mouse in KO[structure].keys():
        theta_morlets_KO[structure][mouse] = np.mean(KO[structure][mouse][theta[0]:theta[1], :, :], axis=(0, 2))


# In[ ]:


print('Saving variables...')
pickle.dump(theta_morlets_WT, open('/home/maspe/filer/SERT/ALL/npys/theta_morlets_WT.dict', 'wb'))
pickle.dump(theta_morlets_KO, open('/home/maspe/filer/SERT/ALL/npys/theta_morlets_KO.dict', 'wb'))

#np.save('/home/maspe/filer/SERT/ALL/npys/theta_morlets_WT.npy', theta_morlets_WT)
#np.save('/home/maspe/filer/SERT/ALL/npys/theta_morlets_KO.npy', theta_morlets_KO)

print('Completed in {:.2f} min.'.format((time.time() - master_clock) / 60))


# In[2]:


#morlets_WT = pickle.load(open('/home/maspe/filer/SERT/ALL/npys/theta_morlets_WT.dict'))


# In[ ]:


#c = cm.jet((c-np.min(c))/(np.max(c)-np.min(c)))
#    ax = plt.gca()
#    for i in np.arange(len(x)-1):
#        ax.plot([x[i],x[i+1]], [y[i],y[i+1]], c=c[i])

