#!/usr/bin/env python
# coding: utf-8

'''

PENULTIMO: TOMA .EPOCHS DENTRO DE NPYS DE CADA RATON
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


# Import required modules
import sys
sys.path.append('/home/maspe/filer/projects/brainsignals/modules')

import numpy as np
from scipy import signal
from matplotlib import pyplot as plt
import wavelets as wl
import time
import pickle


def transform(mice_list=[''], structures_list=[''], n_freq=80, substract_baseline=False, epoch_length= 6, final_fs=1000.0, save_data= False):

    if substract_baseline:
        this_type = 'baselines'
        conditions = ['epochs', 'baselines']
    else:
        this_type = 'no-baselines'
        conditions = ['epochs']

        
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

