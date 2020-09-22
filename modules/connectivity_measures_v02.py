#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  2 17:23:54 2019

@author: carlosmig
Includes some of the common measures of connectivity used in M/EEG. For
more details read: Multi-Dimensional Dynamics of Human Electromagnetic
Brain Activity (Kida et al., 2016)
"""

import numpy as np
from scipy import signal
from scipy import stats


def correlation(data, average = False):
    '''Calculates the pearson correlation between signals
    Parameters
    ----------
    data : Matrix with the time series
           Rows -> oscillators 
           Columns -> time
    average: if 'True' returns the average value across all pairs
    Returns
    -------
    PhaseSynch: matrix/value with correlations
    '''
    corr_mat = np.corrcoef(data)
    
    N = data.shape[0]
    
    if average == True:
        corr_mat = np.mean(corr_mat[np.tril_indices(N, k = -1)])
    
    return(corr_mat)   
    

def coherence(data, average = False):
    '''Calculates the coherence
    Parameters
    ----------
    data : Matrix with the time series
           Rows -> oscillators 
           Columns -> time
    average: if 'True' returns the average value across all pairs
    Returns
    -------
    PhaseSynch: matrix/value with coherences
    '''
    
    N = data.shape[0]
    S = data.shape[1]
    data_mean = np.mean(data, axis = 1)
    data = (data.T - data_mean).T
    
    spectrum = np.fft.fft(data, axis = 1)[:,0:S//2]
    auto_spec = spectrum * np.conj(spectrum)
    Gxx = np.mean(np.repeat(auto_spec.T, N).reshape((S//2,N,N)).T, 2)
    Gyy = Gxx.T
    Gxy = np.mean(spectrum[None,:,:] * np.conj(spectrum[:,None,:]), axis = 2)
    Cxy = np.abs(Gxy)**2 / (Gxx * Gyy)
    Cxy = np.real(Cxy)
    
    if average == True:
        Cxy = np.mean(Cxy[np.tril_indices(N, k = -1)])
    
    return(Cxy)   
    
    
def icoh(data, average = False):
    '''Calculates the imaginary part of the coherence.
    This metric is 'blind' to zero-lag coherence
    Parameters
    ----------
    data : Matrix with the time series
           Rows -> oscillators 
           Columns -> time
    average: if 'True' returns the average value across all pairs
    Returns
    -------
    PhaseSynch: matrix/value with imag coherences
    '''
    
    N = data.shape[0]
    S = data.shape[1]
    analytic_signal = signal.hilbert(data, axis = 1)
    envelopes = np.abs(analytic_signal)
    phases = np.angle(analytic_signal)
    phase_dif = phases[:,None,:] - phases[None,:,:]
    A1 = np.repeat(envelopes.T, N).reshape((S,N,N)).T
    A2 = np.transpose(A1, axes = [1,0,2])
    iCoh = np.mean(A1 * A2 * np.sin(phase_dif), -1) / np.sqrt(np.mean(A1**2, -1) * np.mean(A2**2, -1))
    iCoh = np.abs(iCoh)
    
    if average == True:
        iCoh = np.mean(iCoh[np.tril_indices(N, k = -1)])
    
    return(iCoh)   
    

def order_parameter(data, average = False):  
    '''Calculates the Kuramoto Order Parameter
    Parameters
    ----------
    data : Matrix with the time series
           Rows -> oscillators 
           Columns -> PHASES
    average: if 'True' returns the average value across all pairs
    Returns
    -------
    PhaseSynch: matrix/value with the Kuramoto Order Parameter
    '''
    
    N = data.shape[0]
    phases = data
         
    PhaseSynch = np.mean(np.abs(np.exp(1j * phases[:,None,:]) / 2 + np.exp(1j * phases[None,:,:]) / 2), -1)
    
    if average == True:
        PhaseSynch = np.mean(PhaseSynch[np.tril_indices(N, k = -1)])
    
    return(PhaseSynch)   


def PLV(data, average = False):    
    '''Calculates the Phase Locking Value (PLV)
    ----------
    data : Matrix with the time series
           Rows -> oscillators 
           Columns -> PHASES
    average: if 'True' returns the average value across all pairs
    Returns
    -------
    PLV: matrix/value with the Phase Locking Value
    '''

    N = data.shape[0]
    phases = data
        
    PLV = np.abs(np.mean(np.exp(1j * (phases[:,None,:] - phases[None,:,:])), -1))
    if average == True:
        PLV = np.mean(PLV[np.tril_indices(N, k = -1)])
    
    return(PLV)


def envelope_corr(data, orthogonalize = False, average = False):  
    '''Calculates the Power Envelope Correlation
    Parameters
    ----------
    data : Matrix with the time series
           Rows -> oscillators 
           Columns -> time
    orthogonalize: if 'True' removes signal components sharing the same phases
                  (the linear component) for avoiding zero-lag correlations
    average: if 'True' returns the average value across all pairs
    Returns
    -------
    PhaseSynch: matrix/value with the envelope correlations
    '''
    
    N = data.shape[0]
    data_mean = np.mean(data, axis = 1)
    data = (data.T - data_mean).T
    
    if orthogonalize == True:
        envelope_corr = np.zeros((N, N))
        for i in range(0, N - 1):
            for j in range(i + 1, N):
                data_pair = data[[i,j],:]
                spec_pair = np.fft.fft(data_pair, axis = 1)
                fx, fy = spec_pair[0,:], spec_pair[1,:]
                spec_yx_lin =  np.real(np.sum(fx * np.conj(fy)) / np.sum(fx * np.conj(fx))) * fx
                y_ortho = np.real(np.fft.ifft(fy - spec_yx_lin))
                data_pair = np.row_stack((data[0,:], y_ortho))
                envelopes = np.abs(signal.hilbert(data_pair, axis = 1))
                envelope_corr[i,j] = stats.pearsonr(envelopes[0,:], envelopes[1,:])[0]
        envelope_corr = np.nan_to_num(envelope_corr)
        envelope_corr = envelope_corr + envelope_corr.T
        np.fill_diagonal(envelope_corr, 1)
    
    else:
        envelopes = np.abs(signal.hilbert(data, axis = 1))
        envelope_corr = np.corrcoef(envelopes) 
    
    if average == True:
        envelope_corr = np.mean(envelope_corr[np.tril_indices(N, k = -1)])
    
    return(envelope_corr)   
       
    
def PLI(data, average = False):    
    '''Calculates the Phase Lag Index (PLI) using the phase differences
    ----------
    data : Matrix with the time series
           Rows -> oscillators 
           Columns -> phases
    average: if 'True' returns the average value across all pairs

    Returns
    -------
    PLI: matrix/value with the Phase Lag Index
    '''
    
    N = data.shape[0]
    phases = data
    phases = (phases + 2 * np.pi) % (2 * np.pi)
    
    PLI = np.abs(np.mean(np.sign(np.sin(phases[:,None,:] - phases[None,:,:])), -1))    
   
    np.fill_diagonal(PLI, 0)
    if average == True:
        PLI = np.mean(PLI[np.tril_indices(N, k = -1)])

    return(PLI)
    
    
def PLI2(data1, data2, average = False, method = 'pli'):    
    '''Calculates the Phase Lag Index (PLI) using the cross-spectrum
    An improved index of phase-synchronization for electrophysiological data 
    in the presence of volume-conduction, noise and sample-size bias
    (Vinck et al., 2011)
    ----------
    data : Matrix with the time series
           Rows -> oscillators 
           Columns -> phases
    average: if 'True' returns the average value across all pairs
    method: 'pli' for normal PLI, and 'wpli' for weighted PLI

    Returns
    -------
    PLI: matrix/value with the Phase Lag Index
    
    '''
    
    N1 = data1.shape[0]
    S1 = data1.shape[1]
    #data_mean = np.mean(data, axis = 1)
    #data = (data.T - data_mean).T
    
    spectrum1 = np.fft.fft(data1, axis = 1)[:,0:S//2]
    spectrum2 = np.fft.fft(data2, axis = 1)[:,0:S//2]
    cross_spec = spectrum1[None,:,:] * np.conj(spectrum2[:,None,:])
    if method == 'pli':
        PLI = np.abs(np.mean(np.sign(np.imag(cross_spec)), -1))
    elif method == 'wpli':
        PLI = np.abs(np.mean(np.abs(np.imag(cross_spec)) * np.sign(np.imag(cross_spec)), -1)) \
            / np.mean(np.abs(np.imag(cross_spec)), -1)
    else:
        print('invalid method')
    
    PLI = np.nan_to_num(PLI)
    np.fill_diagonal(PLI, 0)
    if average == True:
        PLI = np.mean(PLI[np.tril_indices(N, k = -1)])

    return(PLI)   
    
    
# tmax = 6
# dt = 1E-4
# N = int(tmax / dt)
# t = np.linspace(0, tmax, N)
# y1 = np.sin(2 * np.pi * 1 * t)
# y2 = np.sin(2 * np.pi * 1 * t + np.pi)
# y3 = np.sin(2 * np.pi * 1 * t)
# y4 = np.sin(2 * np.pi * 4 * t)

# y = np.column_stack((y1,y2,y3,y4))
# phases = np.angle(signal.hilbert(y, axis = 0))


# print(PLI(phases.T, average = False))
# # print(PLI2(y.T, average = False, method = 'pli'))
# # print(PLI2(y.T, average = False, method = 'wpli'))
# print(coherence(y.T, average = False))
# print(icoh(y.T, average = False))








    
    
    








