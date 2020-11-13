#!/usr/bin/env python3
# coding: utf-8

import spectralAnalysis as spectrograms



##################
# IDs = ['SERT1597', 'SERT1659', 'SERT1678', 'SERT1908', 'SERT1984', 'SERT1985', 'SERT2014', 'SERT1665', 'SERT1668', 'SERT2013', 'SERT2018', 'SERT2024'] 


IDs = ['SERT1597', 'SERT1659', 'SERT1665', 'SERT1668'] 


structures = ['mPFC', 'NAC', 'BLA', 'vHip']

IDs_WT = ['SERT1597', 'SERT1659', 'SERT1678', 'SERT1908', 'SERT1984', 'SERT1985', 'SERT2014'] 
IDs_KO = ['SERT1665', 'SERT1668', 'SERT2013', 'SERT2018', 'SERT2024'] 


##########################
#this_type = 'no-baselines'
#save_figs = True

# Setting colormap range
mycolormap = {'mPFC': (-1.5,1.5), 'NAC': (-1.5,1.5), 'BLA': (-1.5,1.5), 'vHip': (-1.5,1.5)}
color_lim = False
this_type = 'baselines'

grandAverages = spectrograms.groupGenotypes(IDs=IDs, IDs_WT=IDs_WT, IDs_KO=IDs_KO, structures=structures, this_type=this_type)

print(grandAverages.keys())

spectrograms.plotGenotypes(grandAverages, mycolormap, genotype='WT', this_type=this_type, color_lim=True, save_figs=True)
spectrograms.plotDifferences(grandAverages, mycolormap, structures=structures, this_type='baselines', color_lim=False, save_figs=True)
