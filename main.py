#!/usr/bin/env python3
# coding: utf-8

import time2frequency as t2f



##################
IDs = ['SERT1597', 'SERT1659', 'SERT1678', 'SERT1908', 'SERT1984', 'SERT1985', 'SERT2014', 'SERT1665', 'SERT1668', 'SERT2013', 'SERT2018', 'SERT2024'] 

IDs_WT = ['SERT1597', 'SERT1659', 'SERT1678', 'SERT1908', 'SERT1984', 'SERT1985', 'SERT2014'] 
IDs_KO = ['SERT1665', 'SERT1668', 'SERT2013', 'SERT2018', 'SERT2024'] 



t2f.transform(mice_list=IDs, structures_list=['mPFC'], substract_baseline=False, epoch_length= 6, final_fs=1000.0, save_data=True)

