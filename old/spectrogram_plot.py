#!/usr/bin/env python
# coding: utf-8

'''
Input: 
    'SERT/SERTXXXX/npys/morlets_epochs.npy'
    
    created from: timefrequency_epochs.py()
    
    Contains:
    * OF data epoched from 3s before to after each center entrance. 
    * Baseline data epoched from xy.dict of 3 s activity in periphery (at 1 KHz??) and multiplied by 30 to 30 KHz.

    Structured as:
    * Dictionary with baselines and OF epochs sorted by structure.
    * Containing (frequencies, channels, time) matrices with:
        - Channels downsampled to 1 KHz
        - All channels with mean subtracted
        - Low-passed at 300 Hz
        - Median subtracted
        - Time/frequency transformed to 99 frequencies
        
    * Saved at 'SERT/SERTXXXX/npys_dir/morlets_epochs.npy

Output:
    '/home/maspe/filer/SERT/ALL/npys/SRS_KO.dict'
'''


### Importing modules
#import matplotlib
#matplotlib.use('Agg')
import os
import pickle
import numpy as np

from matplotlib import pyplot as plt


### Morlet parameters
dt = 1.0 / 1000
time_windows = np.arange(0, 4, dt)
frequencies = np.arange(1, 50, 1)


### Loading dictionaries with spectrograms
os.chdir('/home/maspe/filer/projects/brainsignals')


n_frequencies = 49
time_points = 4000


### For WTs
def loadWT():
    print('Loading WT...')
    SRS_WT = pickle.load(open('DATA/ALL/morlets/SRS_WT.dict', 'rb'))


    grand_average_WT = {key: [] for key in SRS_WT.keys()}
    #grand_average_WT = {'mPFC': [], 'NAC': [], 'BLA': [], 'vHip': []}

    n_mice = 7
    mouse_average = np.empty((n_frequencies, time_points, n_mice))
    for structure in SRS_WT.keys():    
        iteration = 0

        for mouse in SRS_WT[structure].keys():
            #print('Loading structure {} for mouse {}...'.format(structure, mouse))
            mouse_average[:, :, iteration] = np.mean(SRS_WT[structure][mouse], axis=(2,3))
            iteration += 1

        grand_average_WT[structure] = np.mean(mouse_average, axis=2)

    print('Done!')
    return grand_average_WT
    

def loadKO():
    ### For KOs
    print('Loading KO...')
    SRS_KO = pickle.load(open('DATA/ALL/morlets/SRS_KO.dict', 'rb'))
    grand_average_KO = {'mPFC': [], 'NAC': [], 'BLA': [], 'vHip': []}

    n_mice = 5
    mouse_average = np.empty((49, 3000, n_mice))
    for structure in SRS_KO.keys():    
        iteration = 0

        for mouse in SRS_KO[structure].keys():
            print('Loading structure {} for mouse {}...'.format(structure, mouse))
            mouse_average[:, :, iteration] = np.mean(SRS_KO[structure][mouse], axis=(2,3))
            iteration += 1

        grand_average_KO[structure] = np.mean(mouse_average, axis=2)

    print('Done!')
    return grand_average_KO



############
### Plotting

def plotting(grand_average, genotype, colormap={'mPFC':(-2,2),'NAC':(-2,2),'BLA':(-2,2),'vHip':(-2,2)}, save_fig=False):

    for structure in grand_average.keys():
        print("Processing {}".format(structure))
        
        pwr1 = grand_average[structure]
        fmin = min(frequencies)
        fmax = max(frequencies)

        plt.figure(1, figsize=(10, 5))
        plt.clf()

        ax1 = plt.subplot2grid((1, 5),(0, 0),colspan=4)

        plt.imshow(pwr1,cmap='RdBu',vmax=np.max(pwr1),vmin=-np.max(pwr1),
                   extent=(min(time_windows),max(time_windows),fmin,fmax),
                   origin='lower', interpolation='none',aspect='auto')
        plt.colorbar(fraction=0.05,pad=0.02)
        plt.clim(colormap[structure][0], colormap[structure][1])

        plt.axvline(x=2, color='black')
        plt.axhline(y=3, color='red', linestyle='--')
        plt.axhline(y=8, color='red', linestyle='--')
        plt.axhline(y=13, color='red', linestyle='--')
        plt.axhline(y=25, color='red', linestyle='--')

        locs, labels = plt.xticks()    
        plt.xticks(locs, ['-2', '-1.5', '-1', '-0.5', '0', '0.5', '1'], fontsize=14)
        plt.yticks([3, 8, 13, 25], ['3', '8', '13', '25'], fontsize=14)

        ax1.set_xlabel('Time (s)', fontsize=16)
        ax1.set_ylabel('Frequency (Hz)', fontsize=16)
        plt.title('SRS Grand-average for {} in {}'. format(structure, genotype), fontsize=20)

        if save_fig:
            plt.savefig('figs/spectrograms/{}_{}.png'.format(structure, genotype), dpi=150, orientation='landscape')

        plt.close()


