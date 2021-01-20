#!/usr/bin/env python
# coding: utf-8

'''
[TODO] To document

OUTPUT:
    - 'DATA/continuous/{}_channels.statistics'.format(ID)
    - 'DATA/continuous/{}_{}-{}'.format(ID, structure, substract), channels_matrix)
    - 'DATA/continuous/{}_triggers'


    # [PASTED FROM MAIN FUNCTION BELOW]
    # The loop
    # - Load each of the 32 .continuous Openephys files
    # - Subtract the mean
    # - Low-pass it at 300 Hz
    # - Subtract the median by structure
    # - Save it as '/home/maspe/filer/SERT/SERTXXXX/npys/mPFC', for instance.
    #
    # [TODO] Add downsampling option


'''


# Import required modules
import sys
import os
import glob
import numpy as np
import loadEphys as loadRaw

from scipy import signal
from datetime import datetime
#from matplotlib import pyplot as plt
import pandas as pd
import time
import pickle

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
    df = pd.read_csv('RAW/MICE/' + ID + '/csv/canales.csv', header=None)
    channels_locations = df[0].values.tolist()

    # Collect the indexes for each structure
    mPFC = [i for i,x in enumerate(channels_locations) if x == 'mPFC_left']
    NAC  = [i for i,x in enumerate(channels_locations) if x == 'NAC_left']
    BLA  = [i for i,x in enumerate(channels_locations) if x == 'BLA_left']
    vHip = [i for i,x in enumerate(channels_locations) if x == 'vHipp_left']
    all_channels = [i for i,x in enumerate(channels_locations) if x in ['mPFC_left', 'NAC_left', 'BLA_left', 'vHipp_left']]
    
    channels_indexes = {'channels_locations': channels_locations, 'mPFC': mPFC, 'NAC': NAC, 'BLA': BLA, 'vHip': vHip, 'all': all_channels} 
    print("Channels indexes collected!")

    return channels_indexes


# ############################
# ## Main function
def raw_to_npy(IDs=['SERT1597'], structures=['mPFC'], only=[], filter_order=9, detrend=True, substract='median', downsampling=False, load_triggers=False, load_accelerometers=False, save_data=True):

    
    if (substract == 'median_all') & (structures != 'all'):
        print("\n[WARNING] You can't substract=median_all when not selecting structures=all!")
        print("-> substract=median_all changed to substract=median_structure !!\n")
        substract = 'median_structure'
        input("Press ENTER to continue (or CTR+C to abort)...")

        
    # ## Main loop    
    for ID in IDs:
        clock = time.time()
        path_to_continuous = 'RAW/MICE/' + ID + '/continuous/*.continuous'
        electrodes = sorted(glob.glob(path_to_continuous))

        if not detrend:
            print("[Warning] Detrending set to false!")

        # Get channels and structures indexes
        #channels_indexes = 
        this_structures = list()
        this_structures.append(structures)
        print(this_structures)

                
        for structure in this_structures:
            if structures == 'all':
                channels_indexes = getChannIndexes(ID)['all']
            else:
                channels_indexes = getChannIndexes(ID)[structure]

            if len(only):
                channels_indexes = np.array(channels_indexes)[only] - 1
                
            #print("Structure indexes: {}".format(channels_indexes))

            # Selecting electrodes for structure
            print("Picking electrodes for {}".format(structure))

            
            electrodes = np.array(electrodes)[channels_indexes]
            n_channels = len(electrodes)        
            print("# electrodes: {}".format(n_channels))
            print("Electrodes indexes: {}".format(channels_indexes))
            
            # Loads and low-passes all channels of this mice
            iteration = 0
            print("\nLoading channels...")
            for electrode in electrodes:
                channel = loadRaw.loadContinuous(electrode)

                # Low-pass filter
                print('Low-pass filtering (order = {}) at {} Hz...'.format(filter_order, highcut))
                b, a    = butter_bandpass(highcut=highcut, fs=fs, order=filter_order)            
                channel = channel['data'] #[start_OF - points_before : stop_OF]
                channel_mean = np.mean(channel)
                channel_sd = np.std(channel)                

                if detrend:
                    print("Substracting channel's mean...")
                    channel = channel - channel_mean
                    
                channel = signal.filtfilt(b=b, a=a, x=channel)

                if iteration == 0:
                    channels_matrix = np.empty((n_channels, len(channel)))
                    mean_vector = np.empty((n_channels, 1))
                    sd_vector = np.empty((n_channels, 1))

                channels_matrix[iteration, :] = channel
                mean_vector[iteration] = channel_mean
                sd_vector[iteration] = channel_sd


                iteration += 1

            print("\nCalculating standard deviations...")
            corrected_sd_vector = np.empty((n_channels, 1))    
            for i in range(n_channels):
                corrected_sd_vector[i] = np.std(np.delete(channels_matrix, i, axis=0)) 

            channels_statistics = dict(mean=mean_vector, sd=sd_vector, sd_i=corrected_sd_vector)

            print("Calculating rejections...")
            rejection = np.concatenate((np.where(sd_vector < -1.5 * corrected_sd_vector)[0], np.where(sd_vector > 1.5 * corrected_sd_vector)[0]))

                
            #sd = np.std(channels_matrix)
            print("\nChannels mean:\n{}".format(mean_vector))
            print("\nVector of s.d.:\n{}".format(sd_vector))            
            print("\nVector of s.d.(j != i):\n{}".format(corrected_sd_vector))
            print("\nChannels with sd > 1.5 x s.d(j != i):\n{}".format(rejection))
            #print("\nAll channels sd:\n{}".format(np.std(channels_matrix)))
            #rejection = np.where(mean_vector < -2 * sd)

            # print("Rejection: {}".format(rejection))
            # print("S.d.: {}".format(sd))

                
            #print('\nCollecting all channels by structure...')    
            if substract == 'median':
                print("\nComputing channels' median for {}...".format(structure))
                substractor = np.median(channels_matrix)

            if substract == 'mean':
                print("\nComputing channels' median for {}...".format(structure))
                substractor = np.mean(channels_matrix)

            
            print("Substracting channels' {}...".format(substract))
            #print("Channels numbers: {}".format(channels_indexes[structure]))
            channels_matrix = channels_matrix - substractor


            if save_data:
                print("Saving...")
                pickle.dump(channels_statistics, open('DATA/continuous/{}_channels.statistics'.format(ID), 'wb'), protocol=2)
                np.save('DATA/continuous/{}_{}-{}'.format(ID, structure, substract), channels_matrix)
                print('{} saved !!!'.format(structure))


            del [channel, channels_matrix]  
            print('Done!')


        # ## Loading triggers
        if load_triggers:
            print("\nLoading triggers file...")
            triggers_channel = loadRaw.loadContinuous('RAW/MICE/' + ID + '/continuous/AUX/100_ADC1.continuous')

            if save_data:
                print("Saving triggers file...")
                np.save("DATA/continuous/{}_triggers".format(ID), triggers_channel)


        print("\nMouse {} processed in {} min".format(ID, (time.time() - clock) / 60))


    # ## Print function end time
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    print("\nProcess finished at {}".format(current_time))

