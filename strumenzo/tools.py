# -*- coding: utf-8 -*-
"""
Created on Wed Oct 16 12:07:18 2019
@author: enzomarra
This is a collection of useful functions and objects I have used in my analyses over the years.
The Trace class is a 1D numpy array with a few extra ephys related attributes and methods.
The Recording class is an experiment/recording session containing one or more traces and additional
information on the acquisition and experimental condition. Recording is intended as a base to create
more experiment specific derived classes.
"""

import numpy as np
import scipy.signal as sig

# Classic filters
def butter_lowpass(cutoff, fs, order=4):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = sig.butter(order, normal_cutoff, btype='low', analog=True)
    return b, a

def butter_lowpass_filter(array, cutoff, fs, order=4):
    b, a = butter_lowpass(cutoff, fs, order=order)
    y = sig.lfilter(b, a, data)
    return y

def bessel_lowpass(cutoff, fs, order=4):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = sig.bessel(order, normal_cutoff, btype='low', analog=True)
    return b, a

def bessel_lowpass_filter(array, cutoff, fs, order=4):
    b, a = bessel_lowpass(cutoff, fs, order=order)
    y = sig.lfilter(b, a, data)
    return y

def butter_highpass(cutoff, fs, order=2):
    nyq=0.5*fs
    normal_cutoff = cutoff / nyq
    b, a = sig.butter(order, normal_cutoff, btype='high', analog=False)
    return b, a

def butter_highpass_filter(array, cutoff, fs, order=2):
    b, a = butter_highpass(cutoff, fs, order=order)
    y = sig.filtfilt(b, a, data)
    return y

def butter_bandpass(lowcut, highcut, fs, order=5):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = sig.butter(order, [low, high], btype='band')
    return b, a

def notch_filter(array, lowcut, highcut, fs, order=4):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = sig.lfilter(b, a, data)
    return y

# Downsampling functions

def downsample(array, factor):
    """ downsample of a certain factor a 100 samples array downsampled by 2
    will be a 50 samples array of the same time duration"""
    run_mean=running_mean(channel,factor)
    out=[]
    out_size=int(len(channel)/factor)    
    sample_size= float(len(channel))/out_size
    for i in range(out_size):
        out.append(channel[int(floor(i*sample_size))])
    try :
        out=neo.AnalogSignal(signal=np.array(out),units=channel.units, sampling_rate=channel.sampling_rate/factor)
    except: 
        pass
    return out

def downsample_to(array, out_size):
    """keeps the shape of the signal but squeezes it into a fixed-size array"""
    return downsample(channel, int(floor(len(channel)/out_size)))
 

# Finding thresholds without interpolation 

def find_nearest_ind(array,value):
    """finds index of nearest value in a array"""
    array = np.asarray(array)
    ind = (np.abs(array-value)).argmin()
    return ind
  
def find_nearest_sample(array,value):
    """finds sampled value of nearest value in a array"""
    return array[find_nearest_ind(array,value)]

# Finding events

def find_incipit(array): return np.argmax(np.abs(np.diff(trace)))+1



# Classes
  
  
