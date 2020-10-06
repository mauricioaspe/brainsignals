#!/usr/bin/env python
# coding: utf-8

"""
Import continuous data from .continuous Openephys raw files
 
Ideas
- Crear una clase Mouse donde poner todo por ratón?</li>
- Claramente crear una clase con los filtros</li>


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
        
""" 


### Import the required modules
import glob
import sys
import numpy as np
import pandas as pd
import physig as ps
from scipy import signal
from matplotlib import pyplot as plt


if len(sys.argv) > 1:
    ID = sys.argv[1]
else:
    print('Bad file format!')
    exit()


print('Processing mouse {}...'.format(ID))

# Setting working file and paths
files_dir = '~/SERT/DATA/mice/' + ID + '/continuous/'
npys_dir  = '~/SERT/DATA/mice/' + ID + '/npys/'
figs_dir  = '~/SERT/DATA/mice/' + ID + '/figs/'


### Lists all continuous files and read the list of channels
# List all continuous files in working folder
files = sorted(glob.glob(files_dir + '/*.continuous'))

# Read the list of channels
df = pd.read_excel(files_dir + 'canales.xlsx', sheet_name=0, header=None, names=["locs"])
channels_locations = np.array(df['locs'].tolist())
n_channels = len(channels_locations)


### Gets the start and end of the OF and the entrances times 
# Start and end of OF
df = pd.read_excel(files_dir + 'xy.xlsx', sheet_name=0, header=None)
start_OF = int(np.ceil(df[3][0] * 30)) # Times 30 because of sampling at 30 KHz
stop_OF = int(np.floor(df[3][1] * 30))

# Read the entrances times
df = pd.read_excel(files_dir + 'entradas.xlsx', sheet_name=0, header=None, names=["locs"])
entrances_times = np.array(df['locs'].tolist(), dtype='int') * 30 # Times 30 because of sampling at 30 KHz
n_epochs = len(entrances_times)
print('Number of entrances: {}'.format(n_epochs))


### Get indexes and number of channels for each structure. Write a dictionary 'info.npy'

# Output:
# info = {
#     <li>'channels_list'</li>
#     <li>'channels_locs'</li>
#     <li>'n_channels'</li>
#     <li>'mPFC_nchannels'</li>
#     <li>'NAC_nchannels'</li>
#     <li>'BLA_nchannels'</li>
#     <li>'vHip_nchannels'</li>
#     <li>'entrances_times'</li>
#     <li>'n_epochs': number of epochs,</li>
#     <li>'startOF': start of the OF, in (units),</li>
#     <li>'stopOF': stop_OF, in (units)}</li>
# </ul>
# 
# Save it to '/home/maspe/filer/SERT/SERTXXXX/npys/info.npy'


# Collect the indexes for each structure
mPFC_indexes  = [i for i,x in enumerate(channels_locations) if x == 'mPFC_left']
NAC_indexes = [i for i,x in enumerate(channels_locations) if x == 'NAC_left']
BLA_indexes  = [i for i,x in enumerate(channels_locations) if x == 'BLA_left']
vHip_indexes  = [i for i,x in enumerate(channels_locations) if x == 'vHipp_left']

# Get the number of channels for each structure
mPFC_nchannels = len(mPFC_indexes)
NAC_nchannels  = len(NAC_indexes)
BLA_nchannels  = len(BLA_indexes)
vHip_nchannels = len(vHip_indexes)

# Read header info to get sample rate
with open(files[0], 'rb') as f:
    header = ps.readHeader(f)

fs = header['sampleRate']

# Channels list
channels_list = ['ch01', 'ch02', 'ch03', 'ch04', 'ch05', 'ch06', 'ch07', 'ch08', 'ch09', 'ch10',
         'ch11', 'ch12', 'ch13', 'ch14', 'ch15', 'ch16', 'ch17', 'ch18', 'ch19', 'ch20',
         'ch21', 'ch22', 'ch23', 'ch24', 'ch25', 'ch26', 'ch27', 'ch28', 'ch29', 'ch30',
         'ch31', 'ch32']

# Create and save info file
info = {'channels_list': channels_list, 'channels_locs': channels_locations, 'n_channels': n_channels,
        'mPFC_nchannels': mPFC_nchannels, 'NAC_nchannels': NAC_nchannels, 'BLA_nchannels': BLA_nchannels,
        'vHip_nchannels': vHip_nchannels, 'entrances_times': entrances_times, 'n_epochs': n_epochs,
        'startOF': start_OF, 'stopOF': stop_OF}

np.save(npys_dir + 'info.npy', info)

print('Info file: saved\n')


### Create the filter
# Butterpass at 300 Hz, oder 9
# Create filter
def butter_bandpass(highcut, fs, order=5):
    nyq  = 0.5 * fs
    high = highcut / nyq
    b, a = signal.butter(order, high)
    return b, a

# Filter parameters
highcut = 300.0
N       = 9
b, a    = butter_bandpass(highcut, fs, order=N)


### Main loop
#
# The loop
# - Load each of the 32 .continuous Openephys files
# - Subtract the mean
# - Low-pass it at 300 Hz
# - Subtract the median by structure
# - Save it as '/home/maspe/filer/SERT/SERTXXXX/npys/mPFC'


##### Main loop #####
# Loop for loading and low-pass all channels of this mice
iteration = 0
for this_file in files:
    channel = ps.loadContinuous(this_file)

    print('Low-pass filtering (order = {}) at {} Hz...'.format(N, highcut))
    data = channel['data'] #[start_OF - points_before : stop_OF]
    data = signal.filtfilt(b=b, a=a, x=data - np.mean(data),
                           axis=-1, padtype='odd', padlen=None, method='pad', irlen=None)

    
    if iteration == 0:
        data_matrix = np.empty((len(channels_locations), len(data)))          
    
    data_matrix[iteration, :] = data    
    iteration += 1

print('\nCollecting all channels by structure...')    
mPFC = data_matrix[mPFC_indexes, :] - np.median(data_matrix[mPFC_indexes, :], axis=0)
print('mPDF: Done!')
NAC  = data_matrix[NAC_indexes, :]  - np.median(data_matrix[NAC_indexes, :], axis=0)
print('NAC: Done!')
BLA  = data_matrix[BLA_indexes, :]  - np.median(data_matrix[BLA_indexes, :], axis=0)
print('BLA: Done!')
vHip = data_matrix[vHip_indexes, :] - np.median(data_matrix[vHip_indexes, :], axis=0)
print('vHip: Done!')

del [iteration, channel, data, data_matrix]  

print('Saving variables...')
np.save(npys_dir + 'mPFC', mPFC)
np.save(npys_dir + 'NAC', NAC)
np.save(npys_dir + 'BLA', BLA)
np.save(npys_dir + 'vHip', vHip)
print('Done!')

