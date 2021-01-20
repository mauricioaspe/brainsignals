#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pickle


# In[2]:


IDs = {'list_all': ['ID1597', 'ID1659', 'ID1678', 'ID1908', 'ID1984',
                    'ID1985', 'ID2014', 'ID1668', 'ID1665', 'ID2018',
                    'ID2024', 'ID2013'],
       'list_WT': ['ID1597', 'ID1659', 'ID1678', 'ID1908', 'ID1984', 'ID1985', 'ID2014'],
       'list_KO': ['ID1668', 'ID1665', 'ID2018', 'ID2024', 'ID2013'], 
       'dict': {'ID1597': {}, 'ID1659': {}, 'ID1678': {}, 'ID1908': {}, 'ID1984': {},
                'ID1985': {}, 'ID2014': {}, 'ID1668': {}, 'ID1665': {}, 'ID2018': {},
                'ID2024': {}, 'ID2013': {}}
      }


# In[3]:


pickle.dump(IDs, open('IDs.dict', 'wb'), protocol=2)

