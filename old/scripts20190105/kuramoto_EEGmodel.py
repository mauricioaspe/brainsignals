#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 16:59:53 2019

@author: carlosmig
Reference:
Stam, C. J., Nolte, G., & Daffertshofer, A. (2007). Phase lag index: assessment of 
functional connectivity from multi channel EEG and MEG with diminished 
bias from common sources. Human brain mapping, 28(11), 1178-1193.

"""

import numpy as np
import matplotlib.pyplot as plt
import time
import conect_measures as conn
from scipy import stats
from scipy import signal
import SWnetwork

nnodes = 64 #number of nodes/sensors
# CM = np.ones((nnodes, nnodes)) #Connectivity matrix (e.g., structural connectivity)
CM = SWnetwork.SW_network(3, 0.15, nnodes)
K = 0 #Coupling between nodes
D = 0 #Variance, not SD (for noise)
i0 = 8 #The number of shared oscillators for consecutive EEG channels

#Random frequencies from a Lorentzian distribution
#centered in 10 Hz, width equal to 0.5 Hz
# w_center, w_width = 10, 0.5 
np.random.seed(101) #This seed control the values picked from the distribution
w0 = np.random.uniform(8,12,nnodes)

#Integration step, simulation time, total points, time vector
#initial conditions, results -> theta(t)
dt = 1E-3
tmax = 10
tskip = 3
N = int(tmax / dt)
N2 = int((tmax - tskip) / dt)
Nskip = N - N2
t = np.linspace(0, tmax - tskip, N2)
ic = np.random.uniform(-np.pi, np.pi, size = nnodes)
theta_res = np.zeros((N, nnodes))
GT = np.zeros((N, nnodes)) #Ground Truth (no noise, no filtering, no volumetric cond.) 

sqDt = np.sqrt(2 * D / dt) 

#Volumetric conduction matrix
VCM = np.zeros((nnodes, nnodes))
for i in range(0, nnodes):
    j = np.arange(i - i0, i + i0 + 1, 1)
    idx = (j * np.heaviside(-j + nnodes, 0) + (j - nnodes) * np.heaviside(j - nnodes, 1)).astype(int)
    VCM[idx.tolist(),i] = 1
    
    
#Init
t0 = time.time()

np.random.seed(222)
theta = ic
for i in range(0, N):
    theta_res[i] = theta
    theta += dt * (2 * np.pi * w0 + K * np.sum(CM * np.sin(theta - theta[:,None]), 1) + \
            + sqDt * np.random.normal(0,1,nnodes))

#Resolve for GT 
theta = ic
for i in range(0, N):
    GT[i] = theta
    theta += dt * (2 * np.pi * w0 + K * np.sum(CM * np.sin(theta - theta[:,None]), 1))

theta_res = theta_res[Nskip:,:]
GT = GT[Nskip:,:]
    
print(time.time() - t0)
#End

# Amplitude = 1 #No modulation of the amplitude
# Below: modulation of the power -> length = max amplitude
# 1 / freq = oscillation cycles relative to the oscillatory frequency 
length, freq, shape = 1.5, np.repeat(w0, N2).reshape(nnodes, N2).T, 1
tm = np.repeat(t, nnodes).reshape(N2, nnodes)
Amplitude = length * np.exp(-shape * np.cos(w0 * tm)**2)

#Signals without volumetric conduction
signals = Amplitude * np.sin(theta_res) 

#Ground Truth (no noise, no filtering, no volumetric cond.)
GT = Amplitude * np.sin(GT)

#EEG signals (with volumetric conduction)
EEG = np.zeros_like(signals)
for i in range(0,nnodes):
    EEG[:,i] = signals @ VCM[:,i] / (2 * i0 + 1)

#Bessel bandpass filtering of the signals
yf = np.fft.fft(EEG, axis = 0)[0:N2//2,:]
pos = np.argmax(np.abs(yf), axis = 0)
freqs = np.fft.fftfreq(N2, dt)[0:N2//2]
freq = freqs[[int(x) for x in pos]]
Mfreq = np.mean(freq)
# Mfreq = w_center


Fmin, Fmax = Mfreq - 3, Mfreq + 3
a, b = signal.bessel(4, [2 * dt * Fmin, 2 * dt * Fmax], btype = 'bandpass')
Vfilt = signal.filtfilt(a, b, signals, axis = 0)
EEGfilt = signal.filtfilt(a, b, EEG, axis = 0) #EEG filtered signal


plt.figure(1)
plt.clf()
plt.subplot(1,2,1)
plt.plot(t, GT)
plt.xlabel('Time (sec)')
plt.ylabel('f(t)')
plt.subplot(1,2,2)
plt.plot(t, EEG)
plt.xlabel('Time (sec)')
plt.ylabel('EEG')
plt.suptitle('Raw signals')
plt.tight_layout()


plt.figure(2)
plt.clf()
plt.subplot(1,2,1)
plt.plot(t, Vfilt)
plt.xlabel('Time (sec)')
plt.ylabel('f(t)')
plt.subplot(1,2,2)
plt.plot(t, EEGfilt)
plt.xlabel('Time (sec)')
plt.ylabel('EEG')
plt.suptitle('Filtered signals')
plt.tight_layout()


# print(conn.PLV(np.angle(signal.hilbert(Vfilt.T, axis = 1)), average = True))
# print(conn.PLV(np.angle(signal.hilbert(EEGfilt.T, axis = 1)), average = True))
# # print(conn.PLI(np.angle(signal.hilbert(Vfilt.T, axis = 1)), average = True))
# # print(conn.PLI(np.angle(signal.hilbert(EEGfilt.T, axis = 1)), average = True))





