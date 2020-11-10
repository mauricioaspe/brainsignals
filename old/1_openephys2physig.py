#!/usr/bin/env python
# coding: utf-8

'''
The file works from terminal taking the name of the mouse as its only argument.
E.g. Run:

python3 continuous_epochs.py 'SERT1597'

to preprocess the mouse 'SERT1957'


Output: 
    1. Info files
    2. Matrix by structure, with 
        - Channels at 30 KHz
        - All channels with mean subtracted
        - Low-passed at 300 Hz
        - Median subtracted
        - Saved at SERTXXXX/npys_dir + structure.npy
'''


# Import required modules
import sys
import os
import glob
import numpy as np
import physig as ps
from scipy import signal
#from matplotlib import pyplot as plt
import pandas as pd

sys.path.insert(1, "/home/maspe/filer/projects/brainsignals/modules/")
os.chdir('/home/maspe/filer/projects/brainsignals/')


if len(sys.argv) > 1:
    ID = sys.argv[1]
else:
    print('Bad file format!')
    exit()


# Path to files
path_to_files = sorted(glob.glob('DATA/MICE/' + ID + '/'))


# ## Options
IDs = []
structures = ['mPFC', 'NAC', 'BLA', 'vHip']
structures_indexes = indexes
load_triggers = True
load_accelerometers = False
substract = 'median_all'
filter_order = 9
save_data = True


# ## Filter
# Butterpass at 300 Hz
def butter_bandpass(highcut=300.0, fs=30000.0, order=5):
    nyq  = 0.5 * fs
    high = highcut / nyq
    b, a = signal.butter(order, high)
    return b, a



def getChannIndexes(path_to_chanlocs=path_to_files + '/csv/canales.csv'):
    # Read the list of channels
    df = pd.read_csv(path_to_chanlocs)
    print (df)

    #channels_locations = np.array(df['locs'].tolist())
    #n_channels = len(channels_locations)

    # Collect the indexes for each structure
    #mPFC_indexes  = [i for i,x in enumerate(channels_locations) if x == 'mPFC_left']
    #NAC_indexes = [i for i,x in enumerate(channels_locations) if x == 'NAC_left']
    #BLA_indexes  = [i for i,x in enumerate(channels_locations) if x == 'BLA_left']
    #vHip_indexes  = [i for i,x in enumerate(channels_locations) if x == 'vHipp_left']

    return channels_indexes



# ############################
# ## Main loop
# The loop
# - Load each of the 32 .continuous Openephys files
# - Subtract the mean
# - Low-pass it at 300 Hz
# - Subtract the median by structure
# - Save it as '/home/maspe/filer/SERT/SERTXXXX/npys/mPFC', for instance.
#
# [TODO] Add downsampling option
def openephys_to_npy(path_to_continuous=path_to_files + 'continuous/*.continuous'):

    electrodes = sorted(glob.glob(path_to_continuous))
    #triggers = path_to_files + 'continuous/AUX/100_ADC1.continuous'

    # Loads and low-passes all channels of this mice
    iteration = 0
    for electrode in electrodes:

        channel = ps.loadContinuous(electrode)

        # Low-pass filter
        print('Low-pass filtering (order = {}) at {} Hz...'.format(N, highcut))
        data = channel['data'] #[start_OF - points_before : stop_OF]

        b, a    = butter_bandpass(highcut, fs, order=filter_order)

        data = signal.filtfilt(b=b, a=a, x=data - np.mean(data),
                               axis=-1, padtype='odd', padlen=None, method='pad', irlen=None)

        if iteration == 0:
            print("It 0")
            data_matrix = np.empty((len(channels_locations), len(data)))

        data_matrix[iteration, :] = data        

        iteration += 1


    print('\nCollecting all channels by structure...')    


    if substracting == 'median_all':
        print("Computing grand-median...")
        grandMedian = np.median(data_matrix)

        mPFC = data_matrix[mPFC_indexes, :] - grandMedian
        print('mPDF: Done!')
        NAC  = data_matrix[NAC_indexes, :]  - grandMedian
        print('NAC: Done!')
        BLA  = data_matrix[BLA_indexes, :]  - grandMedian
        print('BLA: Done!')
        vHip = data_matrix[vHip_indexes, :] - grandMedian
        print('vHip: Done!')


    if substracting == 'median_structure':
        mPFC = data_matrix[mPFC_indexes, :] - np.median(data_matrix[mPFC_indexes, :], axis=0)
        print('mPDF: Done!')
        NAC  = data_matrix[NAC_indexes, :]  - np.median(data_matrix[NAC_indexes, :], axis=0)
        print('NAC: Done!')
        BLA  = data_matrix[BLA_indexes, :]  - np.median(data_matrix[BLA_indexes, :], axis=0)
        print('BLA: Done!')
        vHip = data_matrix[vHip_indexes, :] - np.median(data_matrix[vHip_indexes, :], axis=0)
        print('vHip: Done!')


    del [iteration, channel, data, data_matrix]  


    # ## Loading triggers
    if load_triggers:
        print("Loading triggers file...")
        triggers_channel = ps.loadContinuous(path_to_files + 'continuous/AUX/100_ADC.continuous')

        np.save('npys/triggers', triggers_channel)


    if save_data:
        print('Saving variables...')
        np.save('npys/mPFC', mPFC)
        np.save('npys/NAC', NAC)
        np.save('npys/BLA', BLA)
        np.save('npys/vHip', vHip)


        
    print('Done!')

