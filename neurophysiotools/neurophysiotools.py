# -*- coding: utf-8 -*-
"""
Created on Wed Oct 16 12:07:18 2019
@author: enzomarra
This is a collection of useful functions and objects I have used in my analyses over the years.
The Trace class is a 1D numpy array with a few extra ephys related attributes and methods.
The Recording class is an experiment/recording session containing one or more traces and additional
information on the acquisition and experimental condition. Recording is intended as a base to create
more experiment specific derived classes. Recommended import neurophysiotools as npt 
"""


import numpy as np
import pandas as pd
import scipy.signal as sig
from numbers import Number


# Collection of filters
def butter_lowpass_filter(array, cutoff, sf, order=4):
    """This function applies a lowpass Butterworth filter to a given array of data.
    Inputs:
        array: the array of data to be filtered. This should be a 1D numpy array.
        cutoff: the desired cutoff frequency of the filter, in Hz. This should be a float.
        sf: the sample frequency of the data, in Hz. This should be a float.
        order: the order of the filter. This should be an integer. The default value is 4.
    Output:
        filtered: the array of data after being filtered with the lowpass Butterworth filter. This will have the same shape as the input array."""
    sos = sig.butter(order, cutoff, btype='lowpass', sf, output='sos', analog=False)
    filtered = sig.sosfilter(sos, array)
    return filtered

def bessel_lowpass_filter(array, cutoff, sf, order=4):
    """This function applies a lowpass Bessel filter to a given array of data.
    Inputs:
        array: the array of data to be filtered. This should be a 1D numpy array.
        cutoff: the desired cutoff frequency of the filter, in Hz. This should be a float.
        sf: the sample frequency of the data, in Hz. This should be a float.
        order: the order of the filter. This should be an integer. The default value is 4.
    Output:
        filtered: the array of data after being filtered with the lowpass Bessel filter. This will have the same shape as the input array."""   
    sos = sig.bessel(order, cutoff, btype='lowpass', sf, output='sos', analog=False)
    filtered = sig.sosfilter(sos, array)
    return filtered

def butter_highpass_filter(array, cutoff, sf, order=2):
    """This function applies a highpass Butterworth filter to a given array of data.
    Inputs:
        array: the array of data to be filtered. This should be a 1D numpy array.
        cutoff: the desired cutoff frequency of the filter, in Hz. This should be a float.
        sf: the sample frequency of the data, in Hz. This should be a float.
        order: the order of the filter. This should be an integer. The default value is 4.
    Output:
        filtered: the array of data after being filtered with the highpass Butterworth filter. This will have the same shape as the input array."""
    sos = sig.butter(order, cutoff, btype='highpass', sf, output='sos', analog=False)
    filtered = sig.sosfilter(sos, array)
    return filtered

def bessel_highpass_filter(array, cutoff, sf, order=4):
    """This function applies a highpass Bessel filter to a given array of data.
    Inputs:
        array: the array of data to be filtered. This should be a 1D numpy array.
        cutoff: the desired cutoff frequency of the filter, in Hz. This should be a float.
        sf: the sample frequency of the data, in Hz. This should be a float.
        order: the order of the filter. This should be an integer. The default value is 4.
    Output:
        filtered: the array of data after being filtered with the highpass Bessel filter. This will have the same shape as the input array."""   
    sos = sig.bessel(order, cutoff, btype='highpass', sf, output='sos', analog=False)
    filtered = sig.sosfilter(sos, array)
    return filtered

def bandpass_filter(array, lowcut, highcut, sf, order=4):
    """This function applies a bandpass Butterworth filter to a given array of data.
    Inputs:
        array: the array of data to be filtered. This should be a 1D numpy array.
        lowcut: the lower bound of the desired frequency range for the filter, in Hz. This should be a float.
        highcut: the upper bound of the desired frequency range for the filter, in Hz. This should be a float.
        sf: the sample frequency of the data, in Hz. This should be a float.
        order: the order of the filter. This should be an integer. The default value is 4.
    Output:
        filtered: the array of data after being filtered with the bandpass Butterworth filter. This will have the same shape as the input array."""    
    sos = sig.butter(order, [lowcut, highcut], btype='bandpass', sf, output='sos', analog=False)
    filtered = sig.sosfilter(sos, array)
    return filtered

def notch_filter(array, notch=50.0, window=1.0, sf, order=4):
    """This function applies a notch filter to a given array of data.
    Inputs:
        array: the array of data to be filtered. This should be a 1D numpy array.
        notch: the center frequency of the notch filter, in Hz. This should be a float. The default value is 50.0.
        window: the width of the notch filter, in Hz. This should be a float. The default value is 1.0.
        sf: the sample frequency of the data, in Hz. This should be a float.
        order: the order of the filter. This should be an integer. The default value is 4.
    Output:
        filtered: the array of data after being filtered with the notch filter. This will have the same shape as the input array."""
    lowcut= notch - (window/2.0)
    highcut= notch +(window/2.0)
    sos = sig.butter(order, [lowcut, highcut], btype='bandstop', sf, output='sos', analog=False)
    filtered = sig.sosfilter(sos, array)
    return filtered

def running_mean(array,window):
    """This function applies a running mean filter to a given array of data.
    Inputs:
        array: the array of data to be filtered. This should be a 1D numpy array.
        window: the size of the window for the running mean filter. This should be an integer.
    Output:
        run_mean: the array of data after being filtered with the running mean filter. This will have the same shape as the input array."""    
    
    avg_mask=np.ones(window) / window
    run_mean=np.convolve(array, avg_mask, 'same')
    return run_mean

# Downsampling functions
def downsample(array, factor):
    """Downsample an array by a given factor.
    Parameters:
    array (array-like): The array to downsample.
    factor (int): The factor by which to downsample the array.
    Returns:
    array: The downsampled array."""
    return sig.decimate(array, factor)

def downsample_to(array, out_size):
    """Downsample an array to a certain length.
    Parameters:
    array (array-like): The array to downsample.
    out_size (int): The length of the array after downsampling.
    Returns:
    array: The downsampled array."""
    return downsample(array, int(np.floor(len(array)/out_size)))
 

# Finding thresholds without interpolation 
def find_threshold_crossings(signal, threshold):
    """Returns the indices of the signal where it crosses the given threshold.
    Parameters:
    signal: 1D numpy array of the signal
    threshold: float, the threshold value to compare against
    Returns:
    crossings: list of integers, the indices of the signal where it crosses the threshold
    """
    # Find the indices where the signal is above the threshold
    above_threshold = np.where(signal > threshold)[0]

    # Initialize a list to store the threshold crossings
    crossings = []

    # Iterate over the above_threshold indices
    for i in above_threshold:
        # If the signal at the current index is above the threshold and the signal at the previous index was below the threshold,
        # add the current index to the crossings list
        if signal[i] > threshold and signal[i-1] < threshold:
            crossings.append(i)
            
    return crossings

def find_closest_index(array, value):
    """Returns the index of the value in the array that is closest to the given value.
    Parameters:
        array: 1D numpy array, the array to search in
        value: float, the value to compare against
    Returns:
        index: integer, the index of the value in the array that is closest to the given value
    """
    # Calculate the absolute difference between the value and each element in the array
    diff = np.abs(array - value)

    # Find the index of the minimum difference
    index = np.argmin(diff)

    return index

def find_closest_value(array, value):
    """Returns the value in the array that is closest to the given value.
    Parameters:
        array: 1D numpy array of values
        value: float, the value to compare against
    Returns:
        closest_value: float, the value in the array that is closest to the given value
    """
    # Return the value at the minimum difference index
    return array[find_closest_index(array,value)]


# Finding events
"""scipy signal has lots of useful functions already"""

def find_incipit(signal, threshold = None):
    """Returns the indices of the start of a rapid change in the signal. A rapid change is defined as a change in the signal that
    exceeds the given threshold. If no threshold is given 2 * signal standard deviation will be used
    Parameters:
    signal: 1D numpy array of the signal
    threshold: float, the threshold value to compare against, if no threshold is given 2 * signal standard deviation will be used
    Returns:
    change_starts: list of integers, the indices of the start of a rapid change in the signal
    """
    if threshold==None: threshold= 2* np.std(signal)
    # Calculate the absolute difference between successive elements in the signal
    diff = np.abs(np.diff(signal))

    # Find the indices where the difference exceeds the threshold
    above_threshold = np.where(diff > threshold)[0]

    # Initialize a list to store the change starts
    change_starts = []

    # Iterate over the above_threshold indices
    for i in above_threshold:
        # If the difference at the current index is above the threshold and the difference at the previous index was below the threshold,
        # add the current index to the change_starts list
        if diff[i] > threshold and diff[i-1] < threshold:
            change_starts.append(i)
            
    return change_starts


# General measures
def coastline(array):
    """This function calculates the coastline of a given neural recording channel.
    Inputs:
        channel: the neural recording channel to be analyzed. This should be a neo.AnalogSignal object.
    Output:
        coastline: the calculated coastline of the channel. This will be a float.
    This function uses the numpy library to calculate the coastline using the formula from Niknazar et al. (2013). 
    The absolute value of the difference between successive samples is then taken, 
    and the sum of these values is returned as the output of the function."""
    return np.sum(np.absolute(np.diff(array)))

# General utilities

def batch_open(folder_name, extension='.'):
    """This function generates a list of files in a given folder with a given extension.
    Inputs:
        folder_name: the name of the folder to be searched for files. This should be a string.
        extension: the desired file extension. This should be a string. The default value is '.'.
    Output:
        rec_list: a list of files in the given folder with the given extension. This will be a list of strings.
    This function first uses the listdir() function to generate a list of all the files in the given folder. 
    It then filters this list to only include files with the desired extension, and returns this filtered list as the output of the function. 
    If the listdir() function is not found in the current environment, it is imported from the os library."""

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


# Wrappers and decorators
"""The section below is a collection of decorators to be used within this module.
It's not ideal to use decorators outside the module so use it as a regular function 
(example 2) if you import it somewhere else
Use smooth_by_2 example are as decorator inside this module:
@smooth_by_2
def myfunction(): return usmoothed_array

You can also apply multiple decorators to a single function by stacking them using the @ symbol. For example:

@decorator_1
@decorator_2
@decorator_3
def some_function(x):
return x

This will apply decorator_1, decorator_2, and decorator_3 to some_function() in that order. 
The output of each decorator will be passed as the input to the next decorator, and the final output will be the result of the decorated function.

However, decorators can become messy when applying across modules so to apply a smoothing to an existing method outside this module
it's advisable to use the target function as argument of the wrapper function. For example:

arr_smoother= smooth_by_2(array.tolist)
arr_smoother() #will return the list of arr.tolist but with running mean"""


def smooth_by_2(func):
    """This function is a decorator that applies a running mean filter with a window size of 2 to the output of another function.
    Inputs:
        func: the function whose output will be filtered with the running mean. This should be a function object.
    Output:
        inner: a new function that wraps the original function and applies the running mean filter to its output.
    The inner function of the decorator takes an arbitrary number of arguments using the *args syntax and passes them to the original function 
    using the func(*args) syntax. It then applies the running mean filter with a window size of 2 to the output of the original function using 
    the running_mean() function, and returns the filtered result."""
    def inner(*args):
        return running_mean(func(*args),2)
    return inner


# Classes
class Trace(np.ndarray):
    """The Trace class is a subclass of numpy.ndarray, meaning it is a modified version of a numpy array that has additional attributes and methods.

    The __new__ method is called when the object is first created and is responsible for creating the object and adding additional attributes to it. 
    The __new__ method takes in the following arguments:
        cls: a reference to the class itself
        input_array: the data that will be stored in the Trace object, which is passed to the np.asarray function to create a numpy array
        sampling_frequency: the sampling frequency of the data (default value is 1)
        signal_units: the units of the data (default value is an empty string)
        channel_id: an identifier for the channel (default value is None)
        pre_filtered: a boolean value indicating whether the data has been filtered (default value is None)

    The __array_finalize__ method is called when the object is created as a view of another object, such as when slicing. 
    It sets default values for the attributes that may not have been specified in the original object.

    The t_axis method returns a time axis for the data based on the sampling frequency and the size of the Trace object. 
    The optional start parameter allows the user to specify a starting time for the axis (default value is 0).

    The to_dict method returns a dictionary representation of the Trace object, with keys for the time axis and the signal data. 
    If the channel_id attribute is not None, it is used as the key for the signal data. 
    If the signal_units attribute is not an empty string, it is used as the key for the signal data. 
    If neither of these conditions are met, the key for the signal data is set to 'signal'.

    The to_df method returns a pandas DataFrame representation of the Trace object, created from the dictionary returned by the to_dict method."""

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

    def to_dict(self):
        if self.channel_id!=None:
            return {'time':self.t_axis(), self.channel_id:self.tolist()}
        elif self.signal_units!=str:
            return {'time':self.t_axis(), self.signal_units:self.tolist()}
        else: 
            return {'time':self.t_axis(),'signal':self.tolist()}

    def to_df(self):
        return pd.DataFrame.from_dict(self.to_dict())

    


class BaseRecording(object):
    """
    A base class for recording objects.

    This class stores signals (which can be a single `Trace` object or a list of `Trace` objects), 
    as well as metadata such as a date and time, a description, and additional information stored in a dictionary. 
    The class provides a method for accessing the values in the information dictionary.

    :param signals: The signals to be stored in the recording object. Must be float in a list, tuple, numpy.ndarray, Trace object, or list of Trace objects.
    :param sampling_frequency: (optional) The sampling frequency for the signals, if provided as a list, tuple, or numpy.ndarray.
    :param signal_units: (optional) The units of the signals, if provided as a list, tuple, or numpy.ndarray.
    :param channel_id: (optional) The channel id for the signals, if provided as a list, tuple, or numpy.ndarray.
    :param pre_filtered: (optional) A boolean indicating whether the signals have been pre-filtered, if provided as a list, tuple, or numpy.ndarray.
    :param rec_id: (optional) A string identifier for the recording.
    :param date_time: (optional) A string representing the date and time of the recording.
    :param description: (optional) A string description of the recording.
    :param **kwargs: (optional) Additional information to be stored in the recording's infos dictionary.
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
  
  
