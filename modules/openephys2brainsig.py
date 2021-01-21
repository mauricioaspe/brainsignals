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
import time
import pickle

sys.path.insert(1, "/home/maspe/filer/projects/brainsignals/modules/")
os.chdir('/home/maspe/filer/projects/brainsignals/')

### Parameters
fs = 30000.0
final_fs = 1000.0
dt = 1 / final_fs


# ## Filter
# Butterpass at 300 Hz
def butter_bandpass(highcut=200.0, fs=30000.0, order=5):
    nyq  = 0.5 * fs
    high = highcut / nyq
    b, a = signal.butter(order, high)
    return b, a


def channels_sanity(channels_matrix):
    mean_vector = np.mean(channels_matrix, axis=0)
    sd_vector = np.std(channels_matrix, axis=0)
    
    n_channels = channels_matrix.shape[0]
    ### Standard deviation of all j != i
    print("\nCalculating standard deviations...")
    corrected_sd = np.empty((n_channels, 1))    
    for i in range(n_channels):
        corrected_sd[i] = np.std(np.delete(channels_matrix, i, axis=0)) 

    ### Calculating rejections according to criteria
    print("Calculating rejections...")
    rejection = np.concatenate((np.where(sd_vector < -1.5 * corrected_sd)[0], np.where(sd_vector > 1.5 * corrected_sd)[0]))
    statistics = dict(mean=mean_vector, sd=sd_vector, sd_i=corrected_sd, rejection=rejection)

    if False:
        print("\nChannels mean:\n{}".format(mean_vector))
        print("\nVector of s.d.:\n{}".format(sd_vector))            
        print("\nVector of s.d.(j != i):\n{}".format(corrected_sd))
        print("\nChannels with sd > 1.5 x s.d(j != i):\n{}".format(rejection))

    return statistics


def sanity_plot(ID, channels_matrix):
    ### Plotting channels
    plt.plot(channels_matrix[0, :])
    plt.savefig('DATA/continuous/{}_chan1.png'.format(ID), dpi=150)
    plt.close()


# ############################
# ## Main function
def raw_to_npy(IDs=['SERT1597'], only=[], highcut=200, filter_order=9, detrend=True, substract='median', downsampling=True, sanity=False, load_triggers=True, load_accelerometers=False, save_data=False):

    now = datetime.now()
    start_time = now.strftime("%H:%M:%S")

    ### Main loop    
    for ID in IDs:
        print("Processing mouse {}".format(ID))
        clock = time.time()
        path_to_continuous = 'RAW/MICE/' + ID + '/continuous/*.continuous'

        electrodes = sorted(glob.glob(path_to_continuous))
        n_channels = len(electrodes)
        
        if not detrend:
            print("[Warning] Detrending set to false!")
        
        if len(only):
            print("[Warning] Loading only electrodes {} !!".format(only))
            electrodes = [electrodes[i-1] for i in only]
                            
        ### ELECTRODES
        # Loads and low-passes all channels of this mice
        print("Loading channels...")        
        iteration = 0
        for electrode in electrodes:
            channel = loadRaw.loadContinuous(electrode)
            sample_rate = channel['header']['sampleRate']
            
            ### Low-pass filter
            print('Low-pass filtering (order = {}) at {} Hz...'.format(filter_order, highcut))
            b, a = butter_bandpass(highcut=highcut, fs=sample_rate, order=filter_order)            
            channel = signal.filtfilt(b=b, a=a, x=channel['data'])

            ### Downsampling
            if downsampling:
                print('Downsampling to 1 KHz...')
                ds_factor = sample_rate / final_fs
                time_points = int(channel.shape[0] / ds_factor)
                channel = signal.resample(x=channel, num=time_points)

            ### Collecting preprocessed channel in numpy array 
            if iteration == 0:
                channels_matrix = np.empty((n_channels, len(channel)))

            channels_matrix[iteration, :] = channel

            iteration += 1
            
        ### Electrodes ending
        print("All channels loaded !!")

        ### Creating sanity statistics
        statistics = channels_sanity(channels_matrix)

        if sanity:    
            sanity_plot(ID, channels_matrix)
            
        ### Detrending
        if detrend:
            print("Substracting channel's mean...")
            channel = channel - statistics['mean']            
    
        ### Substracting mean or median of all channels (INCLUIR A SI MISMO?)
        if substract == 'median':
            print("\nComputing channels' median...")
            substractor = np.median(channels_matrix, axis=0)

        if substract == 'mean':
            print("\nComputing channels' median...")
            substractor = np.mean(channels_matrix, axis=0)

        print("Substracting channels' {}...".format(substract))
        channels_matrix = channels_matrix - substractor

        ### Saving channel's statistics
        if save_data:
            print("Saving...")
            #pickle.dump(channels_statistics, open('DATA/continuous/{}_channels.statistics'.format(ID), 'wb'), protocol=2)
            np.save('DATA/continuous/{}_{}'.format(ID, substract), channels_matrix)
            print('{} saved !!!'.format(ID))

        del [channel, channels_matrix]  
        
        ### Loading triggers
        if load_triggers:
            print("\nLoading triggers file...")
            triggers_channel = loadRaw.loadContinuous('RAW/MICE/' + ID + '/continuous/AUX/100_ADC1.continuous')

            if save_data:
                print("Saving triggers file...")
                np.save("DATA/continuous/{}_triggers".format(ID), triggers_channel)

        print("\nMouse {} processed in {:.2f} min\n\n".format(ID, (time.time() - clock) / 60))

        
    ### Print function end time
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    print("\nProcess started at {} and finished at {}".format(start_time, current_time))

    print("Done!")
