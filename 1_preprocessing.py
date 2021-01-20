#!/usr/bin/env python3
# coding: utf-8

import openephys2brainsig as o2b


IDs = ['ID1597', 'ID1659', 'ID1678', 'ID1908', 'ID1984', 'ID1985', 'ID2014', 'ID1668', 'ID1665', 'ID2018', 'ID2024', 'ID2013']


only = []
o2b.raw_to_npy(IDs=IDs, only=only, highcut=200, filter_order=9, detrend=True, substract='median', downsampling=True, sanity=True, load_triggers=True, load_accelerometers=False, save_data=True)



