#!/usr/bin/env python3
# coding: utf-8

import openephys2brainsig2 as o2b



##################
IDs = ['SERT1597', 'SERT1659'] #, 'SERT1678', 'SERT1908', 'SERT1984', 'SERT1985', 'SERT2014', 'SERT1665', 'SERT1668', 'SERT2013', 'SERT2018', 'SERT2024'] 
#IDs_WT = ['SERT1597', 'SERT1659', 'SERT1678', 'SERT1908', 'SERT1984', 'SERT1985', 'SERT2014'] 
#IDs_KO = ['SERT1665', 'SERT1668', 'SERT2013', 'SERT2018', 'SERT2024'] 

#IDs = ['SERT1597']
#structures = ['mPFC', 'NAC', 'BLA', 'vHip']


##########################
#channels_indexes = o2b.getChannIndexes('SERT1597')
#print(channels_indexes)

o2b.openephys_to_npy(IDs=IDs, structures=['all'], filter_order=9, substract='median', downsampling=False, load_triggers=True, load_accelerometers=False, save_data=True)



