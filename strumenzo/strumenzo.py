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
from numbers import Number


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
    filtered = sig.sosfilter(sos, array)
    return filtered

def notch_filter(array, notch=50.0, window=1.0, sf, order=4):
    """Implements a notch filter using scipy Butterworth bandstop 
    cutting at a given frequency and a window around it on an array 
    sampled at sf sampling frequency. Default is 4th order. 
    Returns the filtered array."""
    lowcut= notch - (window/2.0)
    highcut= notch +(window/2.0)
    sos = sig.butter(order, [lowcut, highcut], btype='bandstop', output='sos', analog=True)
    filtered = sig.sosfilter(sos, array)
    return filtered

def running_mean(array,window):
    """Gets a moving average over a window of int points
    similar to Matlab smooth function but with may
    have more artefacts at start and end. Returns the filtered array."""    
    avg_mask=np.ones(window) / window
    run_mean=np.convolve(array, avg_mask, 'same')
    return run_mean

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

def find_perc(trace,ind,percentage, t_to_bsl= 30,onset='onset'):
    """finds a certain percentage of change for a given event, part_value 
    should be between 0 and 1. For example find_part(trace,70,0.2) 
    will find the index of the 20% value of the event in position 70 of trace. 
    If onset='onset', the analysis is performed on the reversed array to make sure the peak comes 
    before the baseline. t_to_bsl is the max number of points between-baseline and peak, 
    too large may interfere with previous points, must be greater than 9"""
    
    if onset=='onset':event=trace[ind:ind-t_to_bsl:-1] #select the period before the peak but reversed limiting the search
    if onset=='offset':event=trace[ind-t_to_bsl:ind] #select the period after the peak 
    incipit= find_incipit(event) #incipit point to approximate baseline
    bsl=np.mean(event[incipit+int(t_to_bsl/5.):]) #mean baseline
    peak=trace[ind] #value of the peak for readability 
    if bsl>peak: perc_value=bsl+(val_dist(bsl,peak)*percentage) #calculate the value of the trace at the desired percentage
    if bsl<peak: perc_value=bsl-(val_dist(bsl,peak)*percentage)

def find_art(array):
    """find large events that go above threshold on more than one channel at the time
    and returns their position. This only works with 16Channels array end files 
    exported with  16 channels""" 
    #pick channels to be used
    chs=[rec.analogsignals[0].as_array(),rec.analogsignals[15].as_array()]#for now limited to two channels
    art_threshold=0.8  #intended to be relative to the max value
    
    #take both positive and negative for both channels
    #for channel in chs:
    
    art_inds=list(set(supra_thres[0])&set(supra_thres[1]))
    
    
    return art_inds

def find_events(trace, event=''):
    """finds events in trace and returns index of the event max amplitude (peak or trough)."""
    event_inds=[]
    
    event_types=['ap','epsc','epsp']
    peak=np.greater #set the search for positive events
    trough=np.less  #set the search for negative events
    try: #checks the event type is supported
        event_types.index(event.lower())
    except ValueError:
        print("This type of event is not supported. Events currently supported are:\n")
        print("event_types")
    reject=[]#testing only
    #assign polarity (direction of event) and number of points describing the event based on common sampling rates
    if event.lower()=='ap':
        event_points=4 #number of points describing the event, high numbers will reduce false positive but increase false negatives
        polarity= peak #search for a max or a min
        #find all possible peaks 
        event_inds=(sig.argrelextrema(trace,polarity, axis=0, order=(event_points), mode='clip'))
        event_inds=event_inds[0]
        #threshold filter
        event_inds=event_inds[trace[event_inds]>(np.mean(trace)+np.std(trace))]#event_inds[trace[event_inds]>(np.mean(trace)+np.std(trace))]
        
    elif event.lower()=='epsc':
        event_points=20 #number of points describing the event, high numbers will reduce false positive but increase false negatives
        polarity= trough #search for a max or a min
        #find all possible peaks 
        event_inds=(sig.argrelextrema(trace,polarity, axis=0, order=(event_points), mode='clip'))
        event_inds=event_inds[0]
        #threshold filter
        #event_inds=event_inds[trace[event_inds]<(np.mean(trace)-(1*np.std(trace)))]
        x=np.arange(1000)
        trise=20.0;tdecay=200.0
        y=-1.0*(1-np.exp(-x/trise))*np.exp(-x/tdecay)#the argmin is 48
        reject=[]
        for i in event_inds:
            if abs(trace[i]-trace[i-int(event_points*0.4)])<=abs(trace[i]-trace[i+int(event_points*0.4)]): #find_perc(trace,i,0.5, onset='onset')
                reject.append(i)
        return event_inds, reject
       
    elif event.lower()=='epsp':
        event_points=80 #number of points describing the event, high numbers will reduce false positive but increase false negatives
        polarity= peak #search for a max or a min
        #find all possible peaks 
        event_inds=(sig.argrelextrema(trace,polarity, axis=0, order=(event_points), mode='wrap'))
        event_inds=event_inds[0]
        #threshold filter
    return event_inds


# General measures
def coastline(channel):
    """returns the coastline using the formula in Niknazar et al.2013  
    only the array part of the neo analog signals is used"""
    return np.sum(np.absolute(np.diff(channel.as_array()[:,0])))

# General utilities

def batch_open(folder_name, extension='.'):
    try: 
        all_files=listdir(folder_name)
    except:
        from os import listdir
        all_files=listdir(folder_name)

    rec_list=[]
    for file in all_files:
        if file.find(extension)!=-1:
            rec_list.append(file)
    print(rec_list)
    return rec_list

def plot_spikes(array,spk_ind=None,freq=200):
    """plots spikes as red dots on the black trace"""

    time=np.linspace(0.0, len(array)/freq,num=len(array))
    plt.figure;
    plt.plot(time,array,'k')
    plt.plot(time[spk_ind],array[spk_ind],'ro')

# Classes
class Trace(np.ndarray):
    def __new__(cls, input_array, sampling_frequency=1., signal_units=str, channel_id=None, pre_filtered=None):
        # Create the ndarray instance
        obj = np.asarray(input_array).view(cls)
        # Add the new attribute to the created instance
        obj.sampling_frequency = sampling_frequency
        obj.signal_units= signal_units
        obj.channel_id=channel_id
        obj.pre_filtered=pre_filtered
        
        # Return the newly created object
        return obj

    def __array_finalize__(self, obj):
        if obj is None: return
        # Set the default value for the sampling_frequency attribute
        self.sampling_frequency = getattr(obj, 'sampling_frequency', float)
        self.signal_units = getattr(obj, 'signal_units', str)
        self.channel_id = getattr(obj, 'channel_id', None)
        self.pre_filtered = getattr(obj, 'pre_filtered', None)

    def t_axis(self, start=0.0):
        return np.linspace(start,self.size/self.sampling_frequency,self.size)


class BaseRecording(object):
    """
    To be described
    signals must be float in a list, a tuple, a numpy.ndarray, a Trace object or a list of Trace Objects

    """
    infos={}
    def __init__(self, signals, sampling_frequency=None, 
    signal_units=None,channel_id=None, pre_filtered=None,rec_id='', 
    date_time='',description='', **kwargs):
        
        if isinstance(signals,Trace):
            self.signals=signals
        elif isinstance(signals,list) or isinstance(signals,tuple) or isinstance(signals,np.ndarray):
            if all(isinstance(x, Number) for x in signals):
                self.signals=Trace(signals)
                if sampling_frequency!=None: self.signals.sampling_frequency=sampling_frequency
                if signal_units!=None: self.signals.signal_units=signal_units
                if channel_id!=None: self.signals.channel_id=channel_id
                if pre_filtered!=None: self.signals.pre_filtered=pre_filtered
            elif all(isinstance(y, Trace) for y in signals):
                self.signals={}
                assigned_channel=0
                for trace in signals:
                    if trace.channel_id==None:
                        self.signals["_"+str(assigned_channel)]=trace
                    else: self.signals[str(trace.channel_id)]=trace
        else:
            print("The signals argument must be a list, a tuple, a numpy.ndarray, a Trace object or a list of Trace Objects")
        
        self.date_time=date_time
        self.description=description
        self.infos = kwargs




    def get_info(self, key):
        if key in self.infos.keys(): return self.infos.get(key)
        else: return False
  
  
