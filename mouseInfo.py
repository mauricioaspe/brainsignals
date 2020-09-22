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


# #### Get the mouse name and create the file name and folder paths

#if len(sys.argv) > 1:
#    ID = sys.argv[1]
#else:
#    print('Bad file name!')
#    exit()

def getMice(path='/home/diogenes/Projects/SERT/DATA/MICE/'):
    return sorted(os.listdir(path))


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
            'vHip_nchannels': vHip_nchannels
        }


        if save:
            print("Saving info file...")
            time.sleep(1)

            pickle.dump(info, open('/home/diogenes/Projects/SERT/results/info_files/' + mouse + '.info', 'wb'), protocol=2)
            print('Info file: saved !!\n')

        print("Done !!\n")

