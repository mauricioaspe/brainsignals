#!/usr/bin/env python
# coding: utf-8


### Importing modules
import os
os.chdir('/home/maspe/filer/projects/brainsignals')
print("Dir: {}".format(os.getcwd()))

import time
import pickle
import numpy as np
from matplotlib import pyplot as plt


fs = 1000.0

# Morlet parameters
dt = 1 / fs
time_windows = np.arange(0, 6, dt) ###
frequencies = np.arange(1, 80, 1) ###



def groupGenotypes(IDs, IDs_WT, IDs_KO, structures=[], this_type='baselines'):

    grandAverage_WTs = dict()
    grandAverage_KOs = dict()
    clock = time.time()
    iterator_WT = 0
    iterator_KO = 0

    ### Group by genotype and structure
    for mouse in IDs:
        print('Loading {}...'.format(mouse))

        ### Loading data
        npys_dir = '/home/maspe/filer/projects/brainsignals/DATA/MICE/' + mouse + '/npys/'

        data = pickle.load(open(npys_dir + 'morlets-{}.epochs'.format(this_type), 'rb'))

        #print("Data keys: {}".format(data.keys()))
        #print("mPFC shape: {}".format(data['mPFC'].shape))

        if mouse in IDs_WT:
            print("Mouse (WT): {}".format(mouse))
            for structure in structures: #data.keys():
                if iterator_WT == 0:
                    grandAverage_WTs[structure] = data[structure]
                else:
                    grandAverage_WTs[structure] = np.dstack((grandAverage_WTs[structure], data[structure]))

            iterator_WT += 1

        else:
            print("Mouse (KO): {}".format(mouse))
            for structure in structures: #data.keys():
                if iterator_KO == 0:
                    grandAverage_KOs[structure] = data[structure]
                else:
                    grandAverage_KOs[structure] = np.dstack((grandAverage_KOs[structure], data[structure]))

            iterator_KO += 1


    #print("Keys: {}".format(grandAverage_WTs.keys()))
    #print("mPFC shape: {}".format(grandAverage_WTs['mPFC'].shape))
    print('All mice processed in {:.2f} s.'.format(time.time() - clock))
    print('Done!')

    grandAverages = {'WT': grandAverage_WTs, 'KO': grandAverage_KOs}
    #for structure in structures.keys():
    #    grandAverage[structure] = np.mean(structures[structure], axis=2)

    print('Grand averages just created!')

    return grandAverages




def plotGenotypes(grandAverages, mycolormap, genotype='WT', this_type='baselines', color_lim=False, save_figs=False):
    print('Plotting {}...'.format(genotype))
    
    for structure in grandAverages.keys():
        print(grandAverages[genotype][structure].shape)
        pwr1 = np.mean(grandAverages[genotype][structure], axis=2)
        fmin = min(frequencies)
        fmax = max(frequencies)

        plt.figure(1, figsize=(10, 4))
        plt.clf()

        ax1 = plt.subplot2grid((1, 5),(0, 0),colspan=4)


        plt.imshow(pwr1,cmap='jet',vmax=np.max(pwr1),vmin=-np.max(pwr1),
                   extent=(min(time_windows),max(time_windows),fmin,fmax),
                   origin='lower', interpolation='none',aspect='auto')
        plt.colorbar(fraction=0.05,pad=0.02)
        if color_lim:
            plt.clim(mycolormap[structure][0], mycolormap[structure][1])

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
        plt.title('SRS Grand-average for {} in {} ({})'.format(structure, genotype, this_type), fontsize=20)

        if save_figs:
            print('Saving figures for {}...'.format(structure))
            plt.savefig('figs/spectrograms/{}_{}_{}_test.png'.format(structure, genotype, this_type), dpi=150, orientation='landscape')



def plotDifferences(grandAverages, mycolormap, structures=[], this_type='baselines', color_lim=False, save_figs=False):        
    #######################
    print('Plotting differences...')
    for structure in structures:
        print("Structure: {}".format(structure))

        pwr1 = np.mean(grandAverages['KOs'][structure], axis=2) - np.mean(grandAverages['WTs'][structure], axis=2)

        fmin = min(frequencies)
        fmax = max(frequencies)

        plt.figure(1, figsize=(10, 4))
        plt.clf()

        ax1 = plt.subplot2grid((1, 5),(0, 0),colspan=4)


        plt.imshow(pwr1,cmap='jet',vmax=np.max(pwr1),vmin=-np.max(pwr1),
                   extent=(min(time_windows),max(time_windows),fmin,fmax),
                   origin='lower', interpolation='none',aspect='auto')
        plt.colorbar(fraction=0.05,pad=0.02)
        if color_lim:
            plt.clim(mycolormap[structure][0], mycolormap[structure][1])


        plt.axvline(x=3, color='black')
        plt.axhline(y=3, color='black', linestyle='--', linewidth=1)
        plt.axhline(y=8, color='black', linestyle='--', linewidth=1)
        plt.axhline(y=13, color='black', linestyle='--', linewidth=1)
        plt.axhline(y=25, color='black', linestyle='--', linewidth=1)
        plt.axhline(y=50, color='black', linestyle='--', linewidth=1)


        locs, labels = plt.xticks()    
        plt.xticks(locs, ['-3', '-2', '-1', '0', '1', '2', '3'], fontsize=14)
        plt.yticks([1.5, 6, 11, 19, 36, 65, 80], [r'$\delta$', r'$\theta$', r'$\alpha$', r'$\beta$', r'$\gamma_{low}$', r'$\gamma_{high}$', '80'], fontsize=14)


        #ax1.set_yscale('log')
        ax1.set_xlabel('Time (s)', fontsize=16)
        ax1.set_ylabel('Frequency (Hz)', fontsize=16)
        plt.title('SRS Grand-average KO - WT for {} ({})'.format(structure, this_type), fontsize=20)

        if save_figs:
            print('Saving differences...')
            plt.savefig('figs/spectrograms/{}_KO-WT_{}_test.png'.format(structure, this_type), dpi=150, orientation='landscape')






