#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pickle


# In[2]:


IDs = {'list_all': ['SERT1597', 'SERT1659', 'SERT1678', 'SERT1908', 'SERT1984',
                    'SERT1985', 'SERT2014', 'SERT1668', 'SERT1665', 'SERT2018',
                    'SERT2024', 'SERT2013'],
       'list_WT': ['SERT1597', 'SERT1659', 'SERT1678', 'SERT1908', 'SERT1984', 'SERT1985', 'SERT2014'],
       'list_KO': ['SERT1668', 'SERT1665', 'SERT2018', 'SERT2024', 'SERT2013'], 
       'dict': {'SERT1597': {}, 'SERT1659': {}, 'SERT1678': {}, 'SERT1908': {}, 'SERT1984': {},
                'SERT1985': {}, 'SERT2014': {}, 'SERT1668': {}, 'SERT1665': {}, 'SERT2018': {},
                'SERT2024': {}, 'SERT2013': {}}
      }


# In[3]:


pickle.dump(IDs, open('IDs.dict', 'wb'), protocol=2)

