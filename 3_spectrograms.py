#!/usr/bin/env python
# coding: utf-8

'''

'''

### Importing required modules
import spectralActivity as spectra


IDs = ['ID1597', 'ID1659'] #, 'ID1678', 'ID1908', 'ID1984', 'ID1985', 'ID2014', 'ID1668', 'ID1665', 'ID2018', 'ID2024', 'ID2013']


spectra.transform(IDs, n_freq=80, substract_baseline=False, epoch_length=4, final_fs=1000.0, save_data=False)
