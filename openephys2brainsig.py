#!/usr/bin/env python
# coding: utf-8

'''
[TODO] Document


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
import time

sys.path.insert(1, "/home/maspe/filer/projects/brainsignals/modules/")
os.chdir('/home/maspe/filer/projects/brainsignals/')



# ## Filter
# Butterpass at 300 Hz
highcut = 300.0
fs = 30000.0
def butter_bandpass(highcut=highcut, fs=fs, order=5):
    nyq  = 0.5 * fs
    high = highcut / nyq
    b, a = signal.butter(order, high)
    return b, a



# ## Get channels indexes
def getChannIndexes(ID):
    print("Collecting channels for mouse {}...".format(ID))
    # Read the list of channels
    df = pd.read_csv('DATA/MICE/' + ID + '/csv/canales.csv', header=None)
    channels_locations = df[0].values.tolist()

    # Collect the indexes for each structure
    mPFC = [i for i,x in enumerate(channels_locations) if x == 'mPFC_left']
    NAC = [i for i,x in enumerate(channels_locations) if x == 'NAC_left']
    BLA = [i for i,x in enumerate(channels_locations) if x == 'BLA_left']
    vHip = [i for i,x in enumerate(channels_locations) if x == 'vHipp_left']
    all_channels = [i for i,x in enumerate(channels_locations) if x in ['mPFC_left', 'NAC_left', 'BLA_left', 'vHipp_left']]
    
    channels_indexes = {'channels_locations': channels_locations, 'mPFC': mPFC, 'NAC': NAC, 'BLA': BLA, 'vHip': vHip, 'all': all_channels} 
    print("Channels indexes collected!")

    return channels_indexes



# ############################
# ## Main loop
def openephys_to_npy(IDs=['SERT1597'], structures=['mPFC'], filter_order=9, substract='median', downsampling=False, load_triggers=False, load_accelerometers=False, save_data=True):
    
    # The loop
    # - Load each of the 32 .continuous Openephys files
    # - Subtract the mean
    # - Low-pass it at 300 Hz
    # - Subtract the median by structure
    # - Save it as '/home/maspe/filer/SERT/SERTXXXX/npys/mPFC', for instance.
    #
    # [TODO] Add downsampling option


    if (substract == 'median_all') & (structures != 'all'):
        print("\n[WARNING] You can't substract=median_all when not selecting structures=all!")
        print("-> Changing to substract=median_structure !!\n")
        substract = 'median_structure'
        time.sleep(3)

        
    # ## Main loop    
    for ID in IDs:
        clock = time.time()
        path_to_continuous = 'DATA/MICE/' + ID + '/continuous/*.continuous'
        electrodes = sorted(glob.glob(path_to_continuous))

        # Get channels and structures indexes
        channels_indexes = getChannIndexes(ID)

        for structure in structures:
            if structures == 'all':
                structure_indexes = channels_indexes['all']
            else:
                structure_indexes = channels_indexes[structure]

            print("Structure indexes: {}".format(structure_indexes))

            # Selecting electrodes for structure
            print("Picking electrodes for {} #{}".format(structure, structure_indexes))
            electrodes = np.array(electrodes)[structure_indexes]
            n_channels = len(electrodes)        
            print("# electrodes: {}".format(n_channels))
            
            # Loads and low-passes all channels of this mice
            iteration = 0
            print("Loading channels...")
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

                
            #print('\nCollecting all channels by structure...')    
            if substract == 'median':
                print("Computing channels' median for {}...".format(structure))
                substractor = np.median(data_matrix)

            if substract == 'mean':
                print("Computing channels' median for {}...".format(structure))
                substractor = np.mean(data_matrix)

            print("Substracting channels' {}...".format(substract))
            print("Processing {} channels...".format(structure))
            print("Channels numbers: {}".format(channels_indexes[structure]))
            data_matrix = data_matrix - substractor

            if save_data:
                print("Saving...")
                np.save('DATA/MICE/' + ID + '/npys/{}-{}'.format(structure, substract), data_matrix)
                print('{} saved !!!'.format(structure))


            del [iteration, channel, data, data_matrix]  
            print('Done!')


        # ## Loading triggers
        if load_triggers:
            print("Loading triggers file...")
            triggers_channel = ps.loadContinuous('DATA/MICE/' + ID + '/continuous/AUX/100_ADC1.continuous')

            if save_data:
                print("Saving triggers file...")
                np.save('DATA/MICE/' + ID + '/npys/triggers', triggers_channel)


        print("Mouse {} processed in {} min".format(ID, time.time() - clock))


    # ## Print function end time
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    print("Process finished at {}".format(current_time))

