#!/usr/bin/env python
# coding: utf-8

# In[3]:


# Import required modules
import glob
import sys
import numpy as np
import pandas as pd
import physig as ps
from scipy import signal
from matplotlib import pyplot as plt
import wavelets as wl
import time


# In[ ]:


# ID = 'SERT1597'
if len(sys.argv) > 1:
    ID = sys.argv[1]
else:
    print('Bad file format!')
    exit()


# In[ ]:


master_clock = time.time()
print('########################\nProcessing mouse {}...'.format(ID))

# Setting working file and paths
files_dir = '/home/maspe/filer/SERT/' + ID + '/continuous/'
npys_dir  = '/home/maspe/filer/SERT/' + ID + '/npys/'
figs_dir  = '/home/maspe/filer/SERT/' + ID + '/figs/'

# List all continuous files in working folder
files = sorted(glob.glob(files_dir + '/*.continuous'))

# Read the list of channels
df = pd.read_excel(files_dir + 'canales.xlsx', sheet_name=0, header=None, names=["locs"])
channels_locations = np.array(df['locs'].tolist())
n_channels = len(channels_locations)

# Start and end of OF
df = pd.read_excel(files_dir + 'xy.xlsx', sheet_name=0, header=None)
OF_npoints = int(np.floor(df[3][1])) - int(np.ceil(df[3][0]))

start_OF = int(np.ceil(df[3][0] * 30)) # Times 30 because of sampling at 30 KHz
stop_OF = int(np.floor(df[3][1] * 30))

# Read the entrances times
df = pd.read_excel(files_dir + 'entradas.xlsx', sheet_name=0, header=None, names=["locs"])
entrances_times = np.array(df['locs'].tolist(), dtype='int') * 30 # Times 30 because of sampling at 30 KHz
n_epochs = len(entrances_times)
print('Number of entrances: {}\n'.format(n_epochs))

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

#np.save(npys_dir + 'info.npy', info)

#print('Info file: saved\n')

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


# In[ ]:


### Parameters
# Downsampling parameters: final resolution = 1000 Hz
fs = 30000
final_fs  = 1000.0
OF_points = fs * 60 * 10
#ds_factor = fs // final_fs

# Morlet parameters
dt = 1 / final_fs
time_windows = np.arange(0, 600, dt) # 600 por numero segundos en los 10 min de OF
frequencies = np.arange(1, 100, 1)
periods = 1 / (frequencies * dt)
scales = periods / wl.Morlet.fourierwl
n_frequencies = frequencies.shape[0]
time_points = time_windows.shape[0]


# In[ ]:


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
        data_matrix = np.empty((n_channels, len(data)))
              
    data_matrix[iteration, :] = data    
    
    clock = time.time()
    print('Downsampling...')
    data = signal.resample(x=data[start_OF : start_OF + OF_points], num=time_points)
    print('Downsampled in {:.2f} min.'.format((time.time() - clock) / 60))

    clock = time.time()
    print('Morlet transform...')
    if iteration == 0:
        morlet_matrix = np.empty((n_frequencies, len(data), n_channels))
    
    transformed = wl.Morlet(data, scales=scales).getnormpower()
    morlet_matrix[:, :, iteration] = transformed
    print('Transformed in {:.2f} min.\n'.format((time.time() - clock) / 60))
    
    iteration += 1
    
print('\nCollecting all channels and Morlets by structure...')    
#mPFC = data_matrix[mPFC_indexes, :] - np.median(data_matrix[mPFC_indexes, :], axis=0)
mPFC_morlet = morlet_matrix[:, :, mPFC_indexes]
print('mPFC: Done!')

#NAC  = data_matrix[NAC_indexes, :]  - np.median(data_matrix[NAC_indexes, :], axis=0)
NAC_morlet = morlet_matrix[:, :, NAC_indexes]
print('NAC: Done!')

#BLA  = data_matrix[BLA_indexes, :]  - np.median(data_matrix[BLA_indexes, :], axis=0)
BLA_morlet = morlet_matrix[:, :, BLA_indexes]
print('BLA: Done!')

#vHip = data_matrix[vHip_indexes, :] - np.median(data_matrix[vHip_indexes, :], axis=0)
vHip_morlet = morlet_matrix[:, :, vHip_indexes]
print('vHip: Done!')

del [data, transformed, data_matrix, morlet_matrix]  


print('Mouse processed in {:.2f} min.\n'.format((time.time() - master_clock) / 60))

print('Saving variables (this takes some time)...')
np.save(npys_dir + 'mPFC_morlet', mPFC_morlet)
np.save(npys_dir + 'NAC_morlet', NAC_morlet)
np.save(npys_dir + 'BLA_morlet', BLA_morlet)
np.save(npys_dir + 'vHip_morlet', vHip_morlet)
print('Done!\n\n')

