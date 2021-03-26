#!/usr/bin/env python
# coding: utf-8

""" Takes the .npy with raw data from openephys2npy, and create baseline and OF epochs 

[TODO] Review the script that detects the baseline epochs
    
Written by Mauricio Aspé-Sánchez. Copyleft 2020.

""" 

# Import required modules
import time
import pickle
import numpy as np
import pandas as pd

from scipy import signal
from matplotlib import pyplot as plt

import os
import sys

os.chdir('/home/maspe/filer/projects/brainsignals/')
sys.path.insert(0, "modules/biotracker")
sys.path.insert(1, "modules")

import locoAnalysis as loco


# IDs
IDs = ['ID1597', 'ID1659', 'ID1908'] #, 'ID1984'] #'ID1678', 'ID1985', 'ID2014', 'ID1668', 'ID1665', 'ID2018', 'ID2024', 'ID2013']

# TODO Change to make baselines for individual mice!!
#baselines_all = pickle.load(open('DATA/info/ID1597.baselines', 'rb'), encoding='latin1')

# Defining time windows as 3 seconds before plus 3 second after center entrances
fs = 1000.0
seconds_pre = 2
seconds_post = 2
window_pre = int(fs * seconds_pre)
window_post = int(fs * seconds_post)
baseline_length = 3 * int(fs)
n_baselines = 25
substracted = 'median'

extract_epochs = True
plotting_epochs = True
extract_baselines = True

save_data = True



def neuroBehavLocking(files_path, center=0.5, dwell_threshold=seconds_post, iei_threshold=seconds_pre):
    #data = loco.getData(files_path)
    xy = loco.getData(files_path, interpolate=True, upsample=False, plot=False)
    center_vector = loco.centerEvents(xy, center=center, save_data=False)
    path_analysis = loco.pathAnalysis(xy, center_vector, dwell_threshold=dwell_threshold, iei_threshold=iei_threshold, center=center, verbose=True)

    return path_analysis


def getEntrances(path_analysis):
    print("Window pre: {}".format(window_pre))
    print("Window post: {}".format(window_post))                

    # Read the entrances times
    entrances_times = path_analysis['checkin_frame'].to_numpy() # * (100 / 3) # From 30 to 30,000 Hz OJO!!
    entrances_times = entrances_times.astype(int)

    return entrances_times


def plotting(ID, data_epochs, condition='epochs'):
    print('Plotting epochs...')
    n_epochs = data_epochs.shape[2]
    struct_n_rows = np.int(np.ceil(n_epochs / 5.0))
    for channel in range(n_channels):
        plt.figure(figsize=(20,10))
        for epoch in range(n_epochs):
            plt.subplot(struct_n_rows, 5, epoch+1)
            plt.plot(data_epochs[channel, :, epoch])

        plt.savefig('figs/epochs/' + ID + '_' + condition + '_ch' + str(channel) + '.png', dpi=150)
        plt.close()
    

def extractEpochs():
    iterator = 1
    for ID in IDs:
        clock = time.time()
        print("\n\n##############################")
        print("Processing mouse {} (# {})...".format(ID, iterator))

        ### Setting working file and paths
        files_path = './RAW/MICE/' + ID + '/csv/xy.csv'
        figs_dir  = 'figs/epochs/'
        path_analysis = neuroBehavLocking(files_path, center=0.5, dwell_threshold=1.5, iei_threshold=2)
        entrances_times = getEntrances(path_analysis)
        entrances_times = entrances_times * int(1000 / 30)

        n_epochs = len(entrances_times)
        print('Number of entrances: {}'.format(n_epochs))
        print("ENTRANCES (MIN): {}".format(entrances_times / 1000 / 60))

        clock = time.time()

        print("\nDATA:")
        print("Substraction method: {}".format(substracted))
        print("Sample rate: {} Hz".format(fs))
        print("Window: [-{}, {}] s".format(seconds_pre, seconds_post))

        print("\nLoading continuous...")        
        data = np.load('DATA/continuous/{}_{}.npy'.format(ID, substracted))

        hc = data[:,:int(data.shape[1] / 2)]
        data = data[:,int(data.shape[1] / 2):]
        print("DATA SHAPE: {}".format(data.shape[1] / 1000 / 60))
        print("BASELINE SHAPE: {}".format(hc.shape[1] / 1000 / 60))    

        n_channels = data.shape[0]
        data_points = data.shape[1]

        print("Number of channels: {}".format(n_channels))
        print("Channels data points: {}".format(data_points))
        print("Recording length: {:.2f} min".format(data_points / fs / 60))

        if extract_epochs:
            print("\nCollecting epochs...")
            data_epochs = np.zeros((n_channels, window_pre + window_post, n_epochs))
            print("Empty data epoch shape: {}".format(data_epochs.shape))

            if extract_baselines:
                data_baselines = np.zeros((n_channels, window_pre + window_post, n_epochs))
                print("Empty baselines epoch shape: {}".format(data_baselines.shape))        

            for channel in range(n_channels):
                for epoch in range(n_epochs):
                    data_epochs[channel, :, epoch] = data[channel, entrances_times[epoch] - window_pre : entrances_times[epoch] + window_post]

                    if extract_baselines:
                        for this_baseline in range(n_epochs):
                            data_baselines[channel, :, this_baseline] = hc[channel, entrances_times[this_baseline] - window_pre : entrances_times[this_baseline] + window_post]

            epochs = data_epochs
            if extract_baselines:
                all_baselines = data_baselines

        if plotting_epochs:
            plotting(ID, epochs)
            plotting(ID, all_baselines, condition='baselines')        

        iterator += 1

        if save_data:
            print("Saving epoched data...")
            np.save('DATA/epochs/' + ID + '_epochs', epochs)
            np.save('DATA/epochs/' + ID + '_baselines', all_baselines)                      

        print('Epochs extracted in {:.2f} s.\n'.format(time.time() - clock))



    print('Done!')

