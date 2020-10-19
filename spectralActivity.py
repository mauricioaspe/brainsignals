#!/usr/bin/env python
# coding: utf-8

'''

ULTIMO: TOMA .EPOCHS DENTRO DE NPYS DE CADA RATON
        RESAMPLEA, MORLET TRANSFORM, EVENTUALMENTE
        RESTA BASELINE Y GUARDA EN MORLETS.EPOCHS
        DENTRO DE NPYS DIR DE CADA RATON 


Import epoched data from created .npy raw file (same as continuos_epochs) (??)

Input:
    * Epoched data stored in
    'SERT/SERTXXXX/npys/baseline_epochs.npy'

    
Output:
    np.save('/SERT/SERTXXXX/npys/morlets.epochs')


'''


### Importing required modules
import os
import sys
import time
import pickle

import numpy as np
import wavelets as wl

from scipy import signal
from matplotlib import pyplot as plt


### Configuring directories
sys.path.append('/home/maspe/filer/projects/brainsignals/modules')
os.chdir('/home/maspe/filer/projects/brainsignals')
print("Working path changed to {}".format(os.getcwd()))



### Morlet transform
def transform(mice_list=[''], structures_list=[''], n_freq=80, substract_baseline=False, epoch_length= 6, final_fs=1000.0, save_data=False):

    ### Downsampling and band-pass parameters
    ### Parameters
    # Downsampling parameters: final resolution = 1000 Hz
    # Revisar que duracion baseline vs epocas no interfiera en la transformada!
    # Morlet parameters
    dt = 1 / final_fs
    time_windows = np.arange(0, epoch_length, dt)
    frequencies = np.arange(1, n_freq, 1)
    periods = 1 / (frequencies * dt)
    scales = periods / wl.Morlet.fourierwl
    n_frequencies = frequencies.shape[0]

    
    if substract_baseline:
        this_type = 'baselines'
        conditions = ['epochs', 'baselines']
    else:
        this_type = 'no-baselines'
        conditions = ['epochs']


    IDs = {key: {} for key in mice_list}
    structures = {key: {} for key in structures_list}

    iteration = 1
    data_dict = {}
    master_clock = time.time()


    for mouse in IDs.keys():
        clock = time.time()
        print('Loading mouse {} (#{})...'.format(mouse, iteration))
        npys_dir  = '/home/maspe/filer/projects/brainsignals/DATA/MICE/' + mouse + '/npys/'

        with open(npys_dir + mouse + '.info', 'rb') as f:
            info = pickle.load(f, encoding='latin1')

        ### Loading data    
        all_data = pickle.load(open(npys_dir + mouse + '.epochs', 'rb'), encoding="latin1")


        for structure in structures.keys():
            print('\nLoading {}...'.format(structure))
            print("Keys: {}".format(all_data[structure].keys()))
            for condition in conditions: # Iterates in baselines and epochs
                if condition == 'epochs':
                    print('Processing epochs...')
                    time_points = time_windows.shape[0]
                else:
                    if substract_baseline:
                        print('Processing baselines...')
                        time_points = time_windows.shape[0] // 2

                ### Loading the data
                data = all_data[structure][condition]
                print("Datashape pre resampling: {}".format(data.shape))

                ### Downsampling
                print('Downsampling to 1 KHz...')
                data = signal.resample(x=data, num=time_points, axis=1)
                print("Datashape resampled: {}".format(data.shape))

                ### Morlet transform
                clock = time.time()
                print('Morlet transform...')
                n_channels = data.shape[0]
                n_epochs = data.shape[2]
                morlet_matrix = np.empty((n_frequencies, time_points, n_channels, n_epochs))  
                for epoch in range(n_epochs): #range(2)
                    for channel in range(n_channels):
                        morlet_matrix[:, :, channel, epoch] = wl.Morlet(data[channel, :, epoch], scales=scales).getnormpower()            


                print("Morlet matrix: {}".format(morlet_matrix.shape))

                if condition == 'epochs':
                    morlets_epochs    = morlet_matrix
                else:
                    if substract_baseline:
                        morlets_baselines = morlet_matrix

                print("Morlet transform in {} s".format(time.time() - clock))

            ### Obtaining z-scores       
            if substract_baseline:
                print('Transforming to z-score...\n')

                ### Getting average and sd on time dimension
                baseline_mean = np.mean(morlets_baselines, axis=(1,3))
                baseline_sd = np.std(morlets_baselines, axis=(1,3))
                morlets_epochs = (morlets_epochs - baseline_mean[:, None, :, None]) / baseline_sd[:, None, :, None]

            structures[structure] = np.mean(morlets_epochs, axis=3)
            print("Final shape: {}".format(structures[structure].shape))

        iteration += 1

        if save_data:
            print('Saving dictionary...')    
            pickle.dump(structures, open(npys_dir + 'morlets-{}.epochs'.format(this_type), 'wb'), protocol=2)

        print('Mouse processed in {:.2f} min.\n'.format((time.time() - clock) / 60))


    print('All mice processed in {:.2f} min.\n'.format((time.time() - master_clock) / 60))    
    print('Done!')



### Grouping genotypes
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



### Plot by genotype
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



### Plot differences between genotypes
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






