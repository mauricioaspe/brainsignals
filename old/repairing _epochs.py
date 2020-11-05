#!/usr/bin/env python
# coding: utf-8

import pickle
import time
import numpy as np
from matplotlib import pyplot as plt


ID = '2018'
mouse = 'SERT' + ID


clock = time.time()

print('Loading data...')
npys_dir  = '/home/maspe/filer/SERT/' + mouse + '/npys/'
data = pickle.load(open(npys_dir + mouse + '.epochs', 'rb'))              

print('Data loaded in {} s!'.format(time.time() - clock))


condition = 'baselines'
structure = 'vHip'

channel = 5
epoch   = 2

### Where to move
epoch2 = -1

print('Ready to modify {} channel {}, epoch {} in {} of {}'.format(condition, channel, epoch, structure, mouse))
print('Mind matrix shape {} !!!'.format(data[structure][condition].shape))
print('Original   = ({})'.format((channel, epoch)))
print('Replace by = ({})'.format((channel, epoch + epoch2)))

temp_data = data[structure][condition][channel, :, epoch-1 + epoch2]

### Plots
plt.subplot(2,1,1)
plt.plot(data[structure][condition][channel, :, epoch-1])

plt.subplot(2,1,2)
plt.plot(temp_data)
### AÑADIR AL PLOT PROMEDIO Y DESVIACION ESTANDAR !!!


accept = True
if accept:
    data[structure][condition][channel, :, epoch-1] = temp_data


save_data = False
clock = time.time()
if save_data:
    pickle.dump(data, open(npys_dir + mouse + '.epochs', 'wb'), protocol=2)
    
print('Data saved in {} s!'.format(time.time() - clock))



epoch = 8
data[structure]['epochs'][:, :, epoch-1] = data[structure]['epochs'][:, :, 5]

