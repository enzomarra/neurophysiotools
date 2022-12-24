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
def butter_lowpass_filter(array, cutoff, sf, order=4):
    """Implements scipy Butterworth lowpass filter at given cutoff
    on an array sampled at sf sampling frequency. Default is 4th order.
    Returns the filtered array."""
    sos = sig.butter(order, cutoff, btype='lowpass', output='sos', analog=True)
    filtered = sig.sosfilter(sos, array)
    return filtered

def bessel_lowpass_filter(array, cutoff, sf, order=4):
    """Implements scipy Bessel lowpass filter at given cutoff
    on an array sampled at sf sampling frequency. Default is 4th order.
    Returns the filtered array."""       
    sos = sig.bessel(order, cutoff, btype='lowpass', output='sos', analog=True)
    filtered = sig.sosfilter(sos, array)
    return filtered

def butter_highpass_filter(array, cutoff, sf, order=2):
    """Implements scipy Butterworth highpass filter at given cutoff
    on an array sampled at sf sampling frequency. Default is 4th order.
    Returns the filtered array."""
    sos = sig.butter(order, cutoff, btype='highpass', output='sos', analog=True)
    filtered = sig.sosfilter(sos, array)
    return filtered

def bessel_highpass_filter(array, cutoff, sf, order=4):
    """Implements scipy Bessel highpass filter at given cutoff
    on an array sampled at sf sampling frequency. Default is 4th order.
    Returns the filtered array."""       
    sos = sig.bessel(order, cutoff, btype='highpass', output='sos', analog=True)
    filtered = sig.sosfilter(sos, array)
    return filtered

def bandpass_filter(array, lowcut, highcut, sf, order=4):
    """Implements scipy Butterworth bandpass filter at given a
    lowcut as the minimum frequency allowed and highpass as the maximum
    frequency allowed on an array sampled at sf sampling frequency. 
    Default is 4th order. Returns the filtered array."""    
    sos = sig.butter(order, [lowcut, highcut], btype='bandpass', output='sos', analog=True)
    filtered = sig.sosfilter(sos, data)
    return filtered

def notch_filter(array, notch=50.0, window=1.0, sf, order=4):
    """Implements a notch filter using scipy Butterworth bandstop 
    cutting at a given frequency and a window around it on an array 
    sampled at sf sampling frequency. Default is 4th order. 
    Returns the filtered array."""
    lowcut= notch - (window/2.0)
    highcut= notch +(window/2.0)
    sos = sig.butter(order, [lowcut, highcut], btype='bandstop', output='sos', analog=True)
    filtered = sig.sosfilter(sos, data)
    return filtered

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

# General measures
def coastline(channel):
    """returns the coastline using the formula in Niknazar et al.2013  
    only the array part of the neo analog signals is used"""
    return np.sum(np.absolute(np.diff(channel.as_array()[:,0])))

# Classes

class Trace(np.array):

class Recording(object):
    """
    Test
    """
    def __init__(self, traces, rec_id='', date_time='',description='', **kwargs):
  
  
