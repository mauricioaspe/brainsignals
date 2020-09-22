#!/usr/bin/env python
# coding: utf-8


""" The file works from terminal taking the name of the mouse as its only argument.
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
        
Ideas:
    1. Crear una clase Mouse donde poner todo por ratón?
    2. Definitivamente crear una clase con los filtros

"""


### Import the required modules
import sys
sys.path.insert(1, '/home/diogenes/Projects/SERT/modules/')

import glob
import pickle
import time
import os
import numpy as np
import pandas as pd
import physig as ps
from scipy import signal


### Get list of mice in MICE folder
def getMice(path='/home/diogenes/Projects/SERT/DATA/MICE/'):
    return sorted(os.listdir(path))


### Get info of mice
def getInfo(mice_list, save=False):

    for mouse in mice_list:
        print('Processing mouse {}...'.format(mouse))
        print("Save info: {}\n".format(save))
        time.sleep(1)

        # Setting working file and paths
        files_dir = '/home/diogenes/Projects/SERT/DATA/MICE/' + mouse + '/continuous/'
        npys_dir  = '/home/diogenes/Projects/SERT/DATA/MICE/' + mouse + '/npys/'
        figs_dir  = '/home/diogenes/Projects/SERT/DATA/MICE/' + mouse + '/figs/'

        # Start and end of OF from xy.csv
        print("Collecting start and end of OF...")
        time.sleep(1)
        
        df = pd.read_csv(files_dir + 'xy.csv', header=None)
        start_OF = int(np.ceil(df[3][0] * 30)) # Times 30 because of sampling at 30 KHz
        stop_OF = int(np.floor(df[3][1] * 30))

        # Read the entrances times from entradas.csv
        print("Reading entrances times...")
        time.sleep(1)
        
        df = pd.read_csv(files_dir + 'entradas.csv', header=None, names=["locs"])
        entrances_times = np.array(df['locs'].tolist(), dtype='int') * 30 # Times 30 because of sampling at 30 KHz
        n_epochs = len(entrances_times)

        print("Entrances times: {}".format(entrances_times))
        print('Number of entrances: {}'.format(n_epochs))

        # Creating channels list
        print('Loading channels list...')
        channels_list = ['ch01', 'ch02', 'ch03', 'ch04', 'ch05', 'ch06', 'ch07', 'ch08', 'ch09', 'ch10',
                         'ch11', 'ch12', 'ch13', 'ch14', 'ch15', 'ch16', 'ch17', 'ch18', 'ch19', 'ch20',
                         'ch21', 'ch22', 'ch23', 'ch24', 'ch25', 'ch26', 'ch27', 'ch28', 'ch29', 'ch30',
                         'ch31', 'ch32']

        # Read the list of channels from canales.csv
        print('Reading channels locations...')
        time.sleep(1)
        
        df = pd.read_csv(files_dir + 'canales.csv', header=None, names=["locs"])
        channels_locations = np.array(df['locs'].tolist())
        n_channels = len(channels_locations)

        # Collect the indexes for each structure
        print("Collecting structure's index for each channel...")
        time.sleep(1)
        
        mPFC_indexes  = [i for i,x in enumerate(channels_locations) if x == 'mPFC_left']
        NAC_indexes = [i for i,x in enumerate(channels_locations) if x == 'NAC_left']
        BLA_indexes  = [i for i,x in enumerate(channels_locations) if x == 'BLA_left']
        vHip_indexes  = [i for i,x in enumerate(channels_locations) if x == 'vHipp_left']

        # Get the number of channels for each structure
        mPFC_nchannels = len(mPFC_indexes)
        NAC_nchannels  = len(NAC_indexes)
        BLA_nchannels  = len(BLA_indexes)
        vHip_nchannels = len(vHip_indexes)

        # List all continuous files in working folder
        channels_paths = sorted(glob.glob(files_dir + '*.continuous'))
        

        ### Create and save info file
        info = {
            'ID': mouse,
            'files_dirs': [],

            'startOF': start_OF,
            'stopOF': stop_OF,
            'entrances_times': entrances_times,
            'n_epochs': n_epochs,

            'channels_list': channels_list,
            'channels_locs': channels_locations,
            'channels_paths': channels_paths,
            'n_channels': n_channels,

            'mPFC_nchannels': mPFC_nchannels,
            'NAC_nchannels': NAC_nchannels,
            'BLA_nchannels': BLA_nchannels,
            'vHip_nchannels': vHip_nchannels,

            'mPFC_indexes': mPFC_indexes,
            'NAC_indexes': NAC_indexes,
            'BLA_indexes': BLA_indexes,
            'vHip_indexes': vHip_indexes
        }


        if save:
            print("Saving info file...")
            time.sleep(1)

            pickle.dump(info, open('/home/diogenes/Projects/SERT/results/info_files/' + mouse + '.info', 'wb'), protocol=2)
            print('Info file: saved !!\n')

        print("Done !!\n")




# Create filter
def butter_bandpass(highcut, fs, order=5):
    nyq  = 0.5 * fs
    high = highcut / nyq
    b, a = signal.butter(order, high)
    return b, a


##### Main loop #####
# Loop for loading and low-pass all channels of this mice
def openephys2physig(mice_list, highcut=300.0, fs=30000.0, N=9, reference='median', save=False):
    for mouse in mice_list:
        print('Processing mouse {}...'.format(mouse))
        print("Save info: {}\n".format(save))
        time.sleep(1)


        b, a    = butter_bandpass(highcut, fs, order=N)
        ## Main loop

        # The loop</b><br>
        #     <li>Load each of the 32 .continuous Openephys files</li>
        #     <li>Subtract the mean</li>
        #     <li>Low-pass it at 300 Hz</li>
        #     <li>Subtract the median by structure</li>
        #     <li>Save it as '/home/maspe/filer/SERT/SERTXXXX/npys/mPFC'</li>

        info = pickle.load(open('/home/diogenes/Projects/SERT/results/info_files/' + mouse + '.info', 'rb'))
        channels_paths = info['channels_paths']
        channels_locs = info['channels_locs']

        mPFC_indexes = info['mPFC_indexes']
        NAC_indexes = info['NAC_indexes']
        BLA_indexes = info['BLA_indexes']
        vHip_indexes = info['vHip_indexes']        
        
        iteration = 0
        for channel_path in channels_paths[:1]:
            channel = ps.loadContinuous(channel_path)

            print('Low-pass filtering (order = {}) at {} Hz...'.format(N, highcut))
            data = channel['data'] #[start_OF - points_before : stop_OF]
            data = signal.filtfilt(b=b, a=a, x=data - np.mean(data),
                                   axis=-1, padtype='odd', padlen=None, method='pad', irlen=None)


            if iteration == 0:
                data_matrix = np.empty((len(channels_locs), len(data)), dtype='float32')          

            data_matrix[iteration, :] = data    
            iteration += 1


        if reference == 'median':
            mPFC_reference = np.median(data_matrix[mPFC_indexes, :], axis=0)
            NAC_reference = np.median(data_matrix[NAC_indexes, :], axis=0)
            BLA_reference = np.median(data_matrix[BLA_indexes, :], axis=0)
            vHip_reference = np.median(data_matrix[vHip_indexes, :], axis=0)            

            
        print('\nCollecting all channels by structure...')    
        mPFC = data_matrix[mPFC_indexes, :] - mPFC_reference
        print('mPDF: Done!')
        NAC  = data_matrix[NAC_indexes, :]  - NAC_reference
        print('NAC: Done!')
        BLA  = data_matrix[BLA_indexes, :]  - BLA_reference
        print('BLA: Done!')
        vHip = data_matrix[vHip_indexes, :] - vHip_reference
        print('vHip: Done!')

        del [iteration, channel, data, data_matrix, reference_matrix]  


        if save:
            print('Saving variables...')
            np.save(npys_dir + 'mPFC', mPFC)
            np.save(npys_dir + 'NAC', NAC)
            np.save(npys_dir + 'BLA', BLA)
            np.save(npys_dir + 'vHip', vHip)


        print('Done!')

