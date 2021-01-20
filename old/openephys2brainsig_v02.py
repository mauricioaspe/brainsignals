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
from matplotlib import pyplot as plt
#import pandas as pd
import time
import pickle

sys.path.insert(1, "/home/maspe/filer/projects/brainsignals/modules/")
os.chdir('/home/maspe/filer/projects/brainsignals/')

### Downsampling and band-pass parameters                                                                                                            ### Parameters                                                                                                                                       # Downsampling parameters: final resolution = 1000 Hz                                                                                                # Revisar que duracion baseline vs epocas no interfiera en la transformada!                                                                          # Morlet parameters
final_fs = 1000
epoch_length = 6
dt = 1 / final_fs
time_windows = np.arange(0, epoch_length, dt)


# ## Filter
# Butterpass at 300 Hz
highcut = 300.0
fs = 30000.0
def butter_bandpass(highcut=highcut, fs=fs, order=5):
    nyq  = 0.5 * fs
    high = highcut / nyq
    b, a = signal.butter(order, high)
    return b, a


# ############################
# ## Main function
def raw_to_npy(IDs=['SERT1597'], only=[], filter_order=9, detrend=True, substract='median', downsampling=True, load_triggers=True, load_accelerometers=False, save_data=False):

    # ## Main loop    
    for ID in IDs:
        print("Processing mouse {}".format(ID))
        clock = time.time()
        path_to_continuous = 'RAW/MICE/' + ID + '/continuous/*.continuous'

        electrodes = sorted(glob.glob(path_to_continuous))
        print("Loading channels...")        
        
        if len(only):
            print("Only {} !!".format(only))
            electrodes = [electrodes[i] for i in only]
        
        n_channels = len(electrodes)
        
        if not detrend:
            print("[Warning] Detrending set to false!")
        
        # Loads and low-passes all channels of this mice
        iteration = 0
        for electrode in electrodes:
            channel = loadRaw.loadContinuous(electrode)
            
            ### Downsampling
            time_points = int(channel['data'].shape[0] / 30)
            #print("Channel original points: {} ({:.2f} min)".format(channel['data'].shape[0], channel['data'].shape[0] / 30000 / 60))
            print('Downsampling to 1 KHz...')
            channel = signal.resample(x=channel['data'], num=time_points) #, axis=1)
            #print("Channel resampled points: {} ({:.2f} min)".format(channel.shape[0], channel.shape[0] / 1000 / 60))

            # Low-pass filter
            print('Low-pass filtering (order = {}) at {} Hz...'.format(filter_order, highcut))
            b, a    = butter_bandpass(highcut=highcut, fs=fs, order=filter_order)            

            #channel = channel['data']
            channel_mean = np.mean(channel)
            channel_sd = np.std(channel)                

            ### DETREND ANTES O DESPUES DE FILTRADO?
            if detrend:
                print("Substracting channel's mean...\n")
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

        ### Vuelta a iteracion por IDs
        print("All channels loaded !!")
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
            print("\nComputing channels' median...")
            substractor = np.median(channels_matrix)

        if substract == 'mean':
            print("\nComputing channels' median...")
            substractor = np.mean(channels_matrix)

        print("Substracting channels' {}...".format(substract))
        #print("Channels numbers: {}".format(channels_indexes[structure]))
        channels_matrix = channels_matrix - substractor


        plt.plot(channels_matrix[0, :])
        plt.savefig('DATA/continuous/{}_chan1.png'.format(ID), dpi=150)
        plt.close()

        
        if save_data:
            print("Saving...")
            pickle.dump(channels_statistics, open('DATA/continuous/{}_channels.statistics'.format(ID), 'wb'), protocol=2)
            np.save('DATA/continuous/{}_{}'.format(ID, substract), channels_matrix)
            print('{} saved !!!'.format(ID))


        del [channel, channels_matrix]  
        print('Done!')

        
        # ## Loading triggers
        if load_triggers:
            print("\nLoading triggers file...")
            triggers_channel = loadRaw.loadContinuous('RAW/MICE/' + ID + '/continuous/AUX/100_ADC1.continuous')

            if save_data:
                print("Saving triggers file...")
                np.save("DATA/continuous/{}_triggers".format(ID), triggers_channel)


        print("\nMouse {} processed in {} min\n\n".format(ID, (time.time() - clock) / 60))

        
    # ## Print function end time
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    print("\nProcess finished at {}".format(current_time))

