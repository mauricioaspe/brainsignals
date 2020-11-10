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
from datetime import datetime
#from matplotlib import pyplot as plt
import pandas as pd

sys.path.insert(1, "/home/maspe/filer/projects/brainsignals/modules/")
os.chdir('/home/maspe/filer/projects/brainsignals/')

# ## Options
structures = ['mPFC', 'NAC', 'BLA', 'vHip']
#structures_indexes = indexes
load_triggers = True
load_accelerometers = False
substract = 'median_all'
highcut = 300.0
fs = 30000.0
filter_order = 9
save_data = True



# # ## Path to files
# def getIDs(IDs = ['SERT1597']):
#     print("Getting IDs")
#     for ID in IDs:
#         path_to_files = sorted(glob.glob('DATA/MICE/' + ID + '/'))

#         print("ID: {}".format(path_to_files))

#     return path_to_files



# ## Filter
# Butterpass at 300 Hz
def butter_bandpass(highcut=highcut, fs=fs, order=5):
    nyq  = 0.5 * fs
    high = highcut / nyq
    b, a = signal.butter(order, high)
    return b, a



def getChannIndexes(ID):
    # Read the list of channels

    df = pd.read_csv('DATA/MICE/' + ID + '/csv/canales.csv', header=None)

    channels_locations = df[0].values.tolist()

    # Collect the indexes for each structure
    mPFC = [i for i,x in enumerate(channels_locations) if x == 'mPFC_left']
    NAC = [i for i,x in enumerate(channels_locations) if x == 'NAC_left']
    BLA = [i for i,x in enumerate(channels_locations) if x == 'BLA_left']
    vHip = [i for i,x in enumerate(channels_locations) if x == 'vHipp_left']

    channels_indexes = {'channels_locations': channels_locations, 'mPFC': mPFC, 'NAC': NAC, 'BLA': BLA, 'vHip': vHip} 


    print("Channels indexes collected!")

    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    print("Process finished at {}".format(current_time))


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


def openephys_to_npy(IDs):

    for ID in IDs:
        path_to_continuous = 'DATA/MICE/' + ID + '/continuous/*.continuous'
        electrodes = sorted(glob.glob(path_to_continuous))

        channels_indexes = getChannIndexes(ID)
        n_channels = len(channels_indexes['channels_locations'])

        # Loads and low-passes all channels of this mice
        iteration = 0
        for electrode in electrodes:
            channel = ps.loadContinuous(electrode)

            # Low-pass filter
            print('Low-pass filtering (order = {}) at {} Hz...'.format(filter_order, highcut))
            b, a    = butter_bandpass(highcut=highcut, fs=fs, order=filter_order)            
            data = channel['data'] #[start_OF - points_before : stop_OF]
            data = signal.filtfilt(b=b, a=a, x=data - np.mean(data),
                                   axis=-1, padtype='odd', padlen=None, method='pad', irlen=None)


            if iteration == 0:
                data_matrix = np.empty((n_channels, len(data)))

            data_matrix[iteration, :] = data        

            iteration += 1


        print('\nCollecting all channels by structure...')    
        if save_data & (substract == 'median_all'):
            print("Computing all channels' median...")
            grandMedian = np.median(data_matrix)

            print("Substracting all channels' median and saving data...")
            for structure in structures:
                print("Processing {} channels...".format(structure))
                print("Channels numbers: {}".format(channels_indexes[structure]))
                np.save('/npys/{}_all-median'.format(structure), data_matrix[channels_indexes[structure], :] - grandMedian)
                print('{}: Done!'.format(structure))

            print('Done!')


        # if substract == 'median_structure':
        #     print("Substracting median by structure...")
        #     mPFC = data_matrix[channels_indexes['mPFC'], :] - np.median(data_matrix[channels_indexes['mPFC'], :], axis=0)
        #     print('mPDF: Done!')
        #     NAC  = data_matrix[channels_indexes['NAC'], :]  - np.median(data_matrix[channels_indexes['NAC'], :], axis=0)
        #     print('NAC: Done!')
        #     BLA  = data_matrix[channels_indexes['BLA'], :]  - np.median(data_matrix[channels_indexes['BLA'], :], axis=0)
        #     print('BLA: Done!')
        #     vHip = data_matrix[channels_indexes['vHip'], :] - np.median(data_matrix[channels_indexes['vHip'], :], axis=0)
        #     print('vHip: Done!')


        del [iteration, channel, data, data_matrix]  


        # ## Loading triggers
        if load_triggers:
            print("Loading triggers file...")
            triggers_channel = ps.loadContinuous(path_to_files + 'continuous/AUX/100_ADC.continuous')

            print("Saving triggers file...")
            np.save('npys/triggers', triggers_channel)


        print("Mouse {}: DONE!".format(ID))


    # ## Print function end time
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    print("Process finished at {}".format(current_time))

