#!/usr/bin/env python
# coding: utf-8

""" Takes the .npy with raw data from openephys2npy, and create baseline and OF epochs 

[TODO] Review the script that detects the baseline epochs


Input: 
    'SERTXXXX/npys_dir + structure.npy' 
    
    Of .npy matrices by structure, with 
        - Channels at 30 KHz
        - All channels with mean subtracted
        - Low-passed at 300 Hz
        - Median subtracted
        
Output: [CHECKEAR]
    'SERT/SERTXXXX/npys/baselines_epochs.npy'
    
    1. Epoched data from 3s before to after (total = 6 s) each center entrance. 
    2. Epoched baselines (mouse in OF periphery) of 3s lenght (total = 3 s).
    
        
Baseline epochs are obtaines from:
    
    'SERT/ALL/npys/baselines.dict'
    
    * A dictionary with n baselines for each mouse
    * Sampled at 30 Hz 
    
    
Written by Mauricio Aspé-Sánchez. Copyleft 2020.

""" 

# Import required modules
# import glob
# import sys
import time
import pickle
import numpy as np
import pandas as pd

from scipy import signal
from matplotlib import pyplot as plt

import os
import sys

os.chdir('/home/maspe/filer/projects/brainsignals/')
#sys.path.append(r'modules/')
sys.path.insert(0, "modules/micetracker")
sys.path.insert(1, "modules")

import locoAnalysis as loco


# IDs
IDs = ['ID1597', 'ID1659']
# 'SERT1659': {}, 'SERT1665': {}, 'SERT1668': {}, 'SERT1678': {}, 'SERT1908': {}}
# 'SERT1984': {},
# 'SERT1985': {},
# 'SERT2013': {},
# 'SERT2014': {},
# 'SERT2018': {},
# 'SERT2024': {}}


# IDs = pickle.load(open('/home/maspe/filer/scripts/preprocessing/IDs.dict'))['dict']
# TODO Change to make baselines for individual mice!!
baselines_all = pickle.load(open('DATA/info/ID1597.baselines', 'rb'), encoding='latin1')

# Defining time windows as 3 seconds before plus 3 second after center entrances
fs = 1000.0
seconds_pre = 3
seconds_post = 2
window_pre = int(fs) * seconds_pre
window_post = int(fs) * seconds_post
baseline_length = 3 * int(fs)
n_baselines = 25
substracted = 'median'
plotting_epochs = True
save_data = False


epochs = {}
iterator = 1
for ID in IDs:
    clock = time.time()
    print('Processing mouse {} (# {})...'.format(ID, iterator))
    
    ### Setting working file and paths
    mypath = './RAW/MICE/' + ID + '/csv/xy.csv'
    figs_dir  = 'figs/epochs/'


    
    data = loco.getData(mypath)
    xy = loco.getXY(data, normalise=True)   
    center_vector = loco.centerEvents(xy, center=0.5, save_data=False)
    path_analysis = loco.pathAnalysis(xy, center_vector, dwell_treshold=1.5, iei_treshold = 3, center=0.5)

    # Read the entrances times
    entrances_times = path_analysis['entrances'].to_numpy() * (100 / 3) # From 30 to 30,000 Hz OJO!!
    entrances_times = entrances_times.astype(int)
    print("Entrances times: {}".format(entrances_times))

    n_epochs = len(entrances_times)
    print('Number of entrances: {}'.format(n_epochs))

    # Defining baselines as 3 s windows
    #baselines = baselines_all['SERT1597'] * 1000 # From 30 to 30,000 Hz OJO!!

    print("Pre window: {}".format(window_pre))
    print("Post window: {}".format(window_post))                

    structures = {}
    #for structure in structures.keys():
    clock = time.time()

    print("\nDATA:")
    #print("Structure: {}".format(structure))
    print("Substraction method: {}".format(substracted))
    print("Sample rate: {} Hz".format(fs))
    print("Window: [-{}, {}] s".format(seconds_pre, seconds_post))

    print("\nLoading continuous...")        
    data = np.load('DATA/continuous/{}_{}.npy'.format(ID, substracted))
    n_channels = data.shape[0]
    data_points = data.shape[1]

    print("Number of channels: {}".format(n_channels))
    print("Channels data points: {}".format(data_points))
    print("Recording length: {:.2f} min".format(data_points / fs / 60))

    print('Collecting epochs...')
    #### MAKE THEM FROM NAN INSTEAD OF ZEROS; REPLACE WITH VALUES FROM VARIABLE DWELL TIMES
    data_epochs    = np.zeros((n_channels, window_pre + window_post, n_epochs))
    print("Empty data epoch shape: {}".format(data_epochs.shape))

    #data_baselines = np.zeros((n_channels, window_pre, n_baselines))
    #print("Empty baselines epoch shape: {}".format(data_baselines.shape))        

    for channel in range(n_channels):
        # print("Channel: {}".format(channel))
        # sys.stdout.write("Channel: %#d%   \r" % (channel))

        # if channel < n_channels - 1:
        #     sys.stdout.flush()

        for epoch in range(n_epochs):
            data_epochs[channel, :, epoch] = data[channel, entrances_times[epoch] - window_pre : entrances_times[epoch] + window_post]

        # for baseline in range(n_baselines):
        #     data_baselines[channel, :, baseline] = data[channel, baselines[baseline] : baselines[baseline] + baseline_length]


    epochs[ID] = data_epochs
    #structures['baselines'] = data_baselines


    if plotting_epochs:
        print('Plotting epochs...')
        # struct_n_channels = data_epochs.shape[0]
        # struct_n_epochs = data_epochs.shape[2]
        struct_n_rows = np.int(np.ceil(n_epochs / 5.0))
        for channel in range(n_channels):
            plt.figure(figsize=(20,10))
            for epoch in range(n_epochs):
                plt.subplot(struct_n_rows, 5, epoch+1)
                plt.plot(data_epochs[channel, :, epoch])

            plt.savefig(figs_dir + ID + '_ch' + str(channel) + '.png', dpi=100)
            plt.close()


        # print('Plotting baselines...')
        # struct_n_channels = data_baselines.shape[0]
        # struct_n_epochs = data_baselines.shape[2]
        # struct_n_rows = np.int(np.ceil(struct_n_epochs / 5.0))
        # for channel in range(struct_n_channels):
        #     plt.figure(figsize=(20,10))
        #     for epoch in range(struct_n_epochs):
        #         plt.subplot(struct_n_rows, 5, epoch+1)
        #         plt.plot(data_baselines[channel, :, epoch])

        #     plt.savefig(figs_dir + '_baselines2_ch' + str(channel) + '.png', dpi=100, format='png')
        #     plt.close()


    iterator += 1


    # if save_data:
    #     print('Saving epoched dictionary...')
    #     pickle.dump(structures, open('DATA/epochs/' + ID + '.epochs', 'wb'), protocol=2)              

    # print('Epochs extracted in {:.2f} min.\n'.format((time.time() - clock) / 60))


print(epochs)
print('Done!')

