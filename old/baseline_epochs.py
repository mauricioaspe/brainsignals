#!/usr/bin/env python
# coding: utf-8

# # Create epochs for baselines

# In[ ]:


'''
Input:
    Collected epochs in
    'SERT/ALL/npys/xy.dict'
    
    * Sampled at 30 Hz from the start of the 10 min OF.
    
Output:
    Dictionary with n baselines for each mouse 
    
    Saved at
    'SERT/ALL/npys/baselines.dict'
    
    
'''


# In[ ]:


### Importing modules
#import time
import pickle
import numpy as np
#import pandas as pd
#from matplotlib import cm
#from matplotlib import pyplot as plt
#from scipy import signal


# In[ ]:


xy = pickle.load(open('/home/maspe/filer/SERT/ALL/npys/xy.dict', 'rb'))


# In[ ]:


IDs = pickle.load(open('/home/maspe/filer/scripts/preprocessing/IDs.dict'))['dict']
WTs = pickle.load(open('/home/maspe/filer/scripts/preprocessing/IDs.dict'))['list_WT']
KOs = pickle.load(open('/home/maspe/filer/scripts/preprocessing/IDs.dict'))['list_KO']


# The period to extract correspond to 2 seconds with the rat <i>moving</i> in the very border of the OF. This corresponds to the spatial points in $[0, 0.125]$ and $[0.875, 1]$ for each normalised axis.
# 
# Criteria:
# <ul>
#     <li>2 consecutive seconds</li>
#     <li>Mouse always moving</li>
# </ul>

# In[ ]:


n_baselines = 25
baselines = {}
for mouse in IDs.keys():
    x = xy[mouse][0] / np.max(xy[mouse][0])
    y = xy[mouse][1] / np.max(xy[mouse][1])
    
    left  = np.where(x < 0.125)
    right = np.where(x > 0.875)
    up    = np.where(y < 0.125)
    down  = np.where(y > 0.875)
    
    baselines[mouse] = np.random.choice(np.concatenate([left[0], right[0], up[0], down[0]]), n_baselines)


# In[ ]:


pickle.dump(baselines, open('/home/maspe/filer/SERT/ALL/npys/baselines.dict', 'wb'), protocol=2)


# ### Extract baseline epochs

# Extract baseline epochs, defined as 20 time points (at 30 Hz) where the mouse is in the periphery of the OF. It asumes that the animal remains on periphery the next 2 seconds (120 time points). The periphery is defined as a corridor of 0.125 width close the wall of the normalised OF (40 x 40 cm).
# 
# <b>Note:</b> Remove first $x$ sec. 
# 
# Each of these points can be located in the continuous preprocessed .npy file (at 30 KHz) multiplying it by 1000.

# In[ ]:


#windows_length = 2 # sec
#window = windows_length * 60

#baseline_windows = [np.arange(epoch, epoch + window, 1) for epoch in baseline_epochs]

