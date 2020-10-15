#!/usr/bin/env python
# coding: utf-8

'''

FINAL: TOMA MORLET.EPOCH DENTRO DE NPYS DE CADA RATON
       E IMPRIME ESPECTROGRAMA


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
import os
os.chdir('/home/maspe/filer/projects/brainsignals')
print("Dir: {}".format(os.getcwd()))

import time
import pickle
import numpy as np
from matplotlib import pyplot as plt


IDs_all = pickle.load(open('IDs.dict', 'rb'))

save_figs = True

### Main loop
IDs = IDs_all['dict']


IDs_WT = ['SERT1597', 'SERT1659', 'SERT1678', 'SERT1908', 'SERT1984', 'SERT1985', 'SERT2014'] 
IDs_KO = ['SERT1665', 'SERT1668', 'SERT2013', 'SERT2018', 'SERT2024'] 

grandAverage_WTs = dict()
grandAverage_KOs = dict()
clock = time.time()
iterator_WT = 0
iterator_KO = 0


print("List all: {}".format(IDs_all['list_WT']))


### Group by genotype and structure
for mouse in mouse_list:
    print('Loading {}...'.format(mouse))
    
    ### Loading data
    npys_dir = '/home/maspe/filer/projects/brainsignals/DATA/MICE/' + mouse + '/npys/'
    data = pickle.load(open(npys_dir + 'morlets.epochs', 'rb'))
    print("Data keys: {}".format(data.keys()))
    print("mPFC shape: {}".format(data['mPFC'].shape))
    
    if mouse in IDs_all['list_WT']:
        print("Mouse: {}".format(mouse))
        for structure in data.keys():
            if iterator_WT == 0:
                grandAverage_WTs[structure] = data[structure]
            else:
                grandAverage_WTs[structure] = np.dstack((grandAverage_WTs[structure], data[structure]), axis=2)
     
        iterator_WT += 1

print("Keys: {}".format(grandAverage_WTs.keys()))
print("mPFC shape: {}".format(grandAverage_WTs['mPFC'].shape))

# else:
    #     for structure in data.keys():
    #         if iterator_KO == 0:
    #             grandAverage_KOs[structure] = np.mean(data[structure], axis=2)
    #         else:
    #             grandAverage_KOs[structure] = np.dstack((grandAverage_KOs[structure], np.mean(data[structure], axis=2)))
            
    #     iterator_KO += 1
     
        
print('All mice processed in {:.2f} s.'.format(time.time() - clock))
print('Done!')


grandAverage = {'WTs': grandAverage_WTs} #, 'KOs': grandAverage_KOs}
#for structure in structures.keys():
#    grandAverage[structure] = np.mean(structures[structure], axis=2)
    
print('Grand averages just created!')


#npoints = 180000
#n_samples = np.int(npoints / ds_factor) # Downsample to 1000 Hz
fs = 1000.0

# Morlet parameters
dt = 1 / fs
time_windows = np.arange(0, 6, dt) ###
frequencies = np.arange(1, 80, 1) ###
#periods = 1 / (frequencies * dt)
#scales = periods / wl.Morlet.fourierwl
#n_frequencies = frequencies.shape[0]
#time_points = time_windows.shape[0]

# Reducing epoch
#start = 1000
#stop = 4000
#baseline = 2000


print('Plotting...')
# Setting colormap range
mycolormap = {'mPFC': (-1.5,1.5), 'NAC': (-1.5,1.5), 'BLA': (-1.5,1.5), 'vHip': (-1.5,1.5)}

for structure in grandAverage_WTs.keys():
    print(grandAverage_WTs[structure].shape)
    pwr1 = np.mean(grandAverage_WTs[structure], axis=2)
    fmin = min(frequencies)
    fmax = max(frequencies)

    plt.figure(1, figsize=(10, 4))
    plt.clf()

    ax1 = plt.subplot2grid((1, 5),(0, 0),colspan=4)
    
          
    plt.imshow(pwr1,cmap='RdBu',vmax=np.max(pwr1),vmin=-np.max(pwr1),
               extent=(min(time_windows),max(time_windows),fmin,fmax),
               origin='lower', interpolation='none',aspect='auto')
    plt.colorbar(fraction=0.05,pad=0.02)
    #plt.clim(mycolormap[structure][0], mycolormap[structure][1])
    
    plt.axvline(x=3, color='black')
    plt.axhline(y=3, color='black', linestyle='--', linewidth=1)
    plt.axhline(y=8, color='black', linestyle='--', linewidth=1)
    plt.axhline(y=13, color='black', linestyle='--', linewidth=1)
    plt.axhline(y=25, color='black', linestyle='--', linewidth=1)
    plt.axhline(y=50, color='black', linestyle='--', linewidth=1)

   
    locs, labels = plt.xticks()    
    plt.xticks(locs, ['-3', '-2', '-1', '0', '1', '2', '3'], fontsize=14)
    plt.yticks([1.5, 6, 11, 19, 36, 65, 80], [r'$\delta$', r'$\theta$', r'$\alpha$',
                                    r'$\beta$', r'$\gamma_{low}$', r'$\gamma_{high}$', '80'], fontsize=14)


    #ax1.set_yscale('log')
    ax1.set_xlabel('Time (s)', fontsize=16)
    ax1.set_ylabel('Frequency (Hz)', fontsize=16)
    plt.title('SRS Grand-average for {} in WT'.format(structure), fontsize=20)

    if save_figs:
        print('Saving figures for {}...'.format(structure))
        plt.savefig('figs/spectrograms/{}_WT.png'.format(structure),
                dpi=150, orientation='landscape')


# #######################
# print('Plotting KOs...')
# for structure in grandAverage_KOs.keys():
#     pwr1 = np.mean(grandAverage_KOs[structure], axis=2)
#     fmin = min(frequencies)
#     fmax = max(frequencies)

#     plt.figure(1, figsize=(10, 4))
#     plt.clf()

#     ax1 = plt.subplot2grid((1, 5),(0, 0),colspan=4)
    
          
#     plt.imshow(pwr1,cmap='RdBu',vmax=np.max(pwr1),vmin=-np.max(pwr1),
#                extent=(min(time_windows),max(time_windows),fmin,fmax),
#                origin='lower', interpolation='none',aspect='auto')
#     plt.colorbar(fraction=0.05,pad=0.02)
#     plt.clim(mycolormap[structure][0], mycolormap[structure][1])

#     plt.axvline(x=3, color='black')
#     plt.axhline(y=3, color='black', linestyle='--', linewidth=1)
#     plt.axhline(y=8, color='black', linestyle='--', linewidth=1)
#     plt.axhline(y=13, color='black', linestyle='--', linewidth=1)
#     plt.axhline(y=25, color='black', linestyle='--', linewidth=1)
#     plt.axhline(y=50, color='black', linestyle='--', linewidth=1)

   
#     locs, labels = plt.xticks()    
#     plt.xticks(locs, ['-3', '-2', '-1', '0', '1', '2', '3'], fontsize=14)
#     plt.yticks([1.5, 6, 11, 19, 36, 65, 80], [r'$\delta$', r'$\theta$', r'$\alpha$',
#                                     r'$\beta$', r'$\gamma_{low}$', r'$\gamma_{high}$', '80'], fontsize=14)

#     #ax1.set_yscale('log')
#     ax1.set_xlabel('Time (s)', fontsize=16)
#     ax1.set_ylabel('Frequency (Hz)', fontsize=16)
#     plt.title('SRS Grand-average for {} in KO'.format(structure), fontsize=20)

    
#     if save_figs:
#         print('Saving figures for KO...')
#         plt.savefig('figs/spectrograms/{}_KO_periphery_baseline3.png'.format(structure),
#                 dpi=150, orientation='landscape')



# print('Plotting difference...')
# for structure in ['mPFC', 'NAC', 'BLA', 'vHip']:
#     print("Structure: {}".format(structure))

#     pwr1 = np.mean(grandAverage['KOs'][structure], axis=2) - np.mean(grandAverage['WTs'][structure], axis=2)
#     fmin = min(frequencies)
#     fmax = max(frequencies)

#     plt.figure(1, figsize=(10, 4))
#     plt.clf()

#     ax1 = plt.subplot2grid((1, 5),(0, 0),colspan=4)
    
          
#     plt.imshow(pwr1,cmap='RdBu',vmax=np.max(pwr1),vmin=-np.max(pwr1),
#                extent=(min(time_windows),max(time_windows),fmin,fmax),
#                origin='lower', interpolation='none',aspect='auto')
#     plt.colorbar(fraction=0.05,pad=0.02)
#     plt.clim(mycolormap[structure][0], mycolormap[structure][1])

    
#     plt.axvline(x=3, color='black')
#     plt.axhline(y=3, color='black', linestyle='--', linewidth=1)
#     plt.axhline(y=8, color='black', linestyle='--', linewidth=1)
#     plt.axhline(y=13, color='black', linestyle='--', linewidth=1)
#     plt.axhline(y=25, color='black', linestyle='--', linewidth=1)
#     plt.axhline(y=50, color='black', linestyle='--', linewidth=1)

   
#     locs, labels = plt.xticks()    
#     plt.xticks(locs, ['-3', '-2', '-1', '0', '1', '2', '3'], fontsize=14)
#     plt.yticks([1.5, 6, 11, 19, 36, 65, 80], [r'$\delta$', r'$\theta$', r'$\alpha$',
#                                     r'$\beta$', r'$\gamma_{low}$', r'$\gamma_{high}$', '80'], fontsize=14)


#     #ax1.set_yscale('log')
#     ax1.set_xlabel('Time (s)', fontsize=16)
#     ax1.set_ylabel('Frequency (Hz)', fontsize=16)
#     plt.title('SRS Grand-average KO - WT for {}'.format(structure), fontsize=20)

#     if save_figs:
#         print('Saving figures for KO...')
#         plt.savefig('/home/maspe/filer/SERT/ALL/figs/SRS/{}_KO-WT.png'.format(structure),
#                 dpi=150, orientation='landscape')






