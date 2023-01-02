# NeurophysioTools

A collection of Python functions and classes for the analysis of neurophysiological data.

This is a simple module intended for prototyping your neurophysiological analyses or in case you need to start your own project from scratch. The module requires only a limited number of packages, namely:  *numpy* and *scipy*; however, to plot spectrograms and convert **Trace** to a dataframe *pandas* and *plotly* are also required. There are many packages and modules out there that are far more sophisticated and better maintained but often have a steeper learning curve. Here, I am providing  simple classes and a handful of functions to make standard analyses slightly easier to write with only rudimentary knowledge of _numpy_ and _scipy_. However, a theoretical knowledge of signal analysis is strongly reccomended.   

To install NeurophysioTools, use pip:

    pip install neurophysiotools

Usage

    import neurophysiotools as nt

## Use the provided functions and classes for analyzing neurophysiological data

Here is a list of classes and functions, use help(FUNCTION_NAME) for more detailed documentation for each function, the intended use should be clear enough. Here, I will provide a little detail on the intended use of the classes and decorators. 


### Classes
These are the most useful elemts of the **neurophysiotools** module and one can read about the practical aspects of their use with help(CLASS_NAME). Their intended use, with examples, can be found below.
- **Trace**: Subclass of numpy.ndarray for storing and analyzing 1D neural signals.
- **BaseRecording**: Base Class for storing traces and additional information about a recording. 

The **Trace** object is a subclass of **numpy.ndarray**, so it inherits its attributes and methods, as well as additional ones usefult to analyse a physiological signal. Here, a physiological signal can represent a variety of  things: intracellular or extracellular electrophysiology channel, grey values of a time lapse microscopy acquisition, accelerometry channel or even pose estimation outputs. The common element is that the signal can be stored in a 1D array and, to make full use of the methods, has been sampled at evenly spaced intervals. The complete list of **Trace** attributes and methods can be obtained using the *help* function. **Trace** intended use is as follow:
    
    import neurophysiotools as nt
    import numpy as np
    
    sf= 100.0 #sampling frequency in Hz
    units='Volts' #label for signal units SI
    id= 'Ch1' #label for the channel
    filters=[0.1,500] #range of frequencies in the recording normally set in hardware
    import numpy as np; signal_array= np.random.random(10000) 
    trace = nt.Trace(signal_array, sampling_rate=sf, signal_units=units, channel_id=id, pre_filtered=filters)
    type(trace)
    out: <class 'neurophysiotools.Trace'>

Given the variety of physiological recordings possible, here a very simple **BaseRecording** object is provided, the complete list of its attributes and methods can be obtained using the *help* function. The important point to keep in mind is that **BaseRecording** is not intended to be used in its *base* format but to be the parent of a customised object. Here is an example of how to use **BaseRecording** to write an object that better suits your data. In this examples, two arrays of value need to be manipulated to generate a biologically meaningful signal.

    def standard_norm(array1,array2):
        #do something

    class Photometry_Recording(nt.BaseRecording):
        
        infos={}
        def __init__(self, ch490,ch405, normalization="standard", sampling_rate=1017.25, 
        rec_id='Photometry Exp',date_time='',description='', **kwargs):
            
            self.normalization=normalization
            photosignals=[nt.Trace(ch490, sampling_rate=sampling_rate,signal_units='AU',channel_id='ch490')
                                , nt.Trace(ch405, sampling_rate=sampling_rate, signal_units='AU',channel_id='ch405')]

            if self.normalization=='standard':
                photosignals.append(nt.Trace(standard_norm(ch490,ch405,sf=sampling_rate, settling=10),signal_units='DeltaF',channel_id='DeltaF'))
            elif self.normalization=='':
                            photosignals.append(nt.Trace(np.ones(len(ch490)), sampling_rate=sampling_rate,signal_units='DeltaF',channel_id='DeltaF'))
            
            super().__init__(photosignals)
        

        def technique_specific_method(self):
        """Details"""


### Decodators
These decorators are intended to be used when defining methods and not outside the module, at least not as decorators. Please see intended use of **array2trace** below, as an example.
- array2trace : changes the output *numpy.ndarray* of a function to a *Trace* object, the input of the changed function must be a **Trace**.
- smooth_by_2 : applies a running mean with window 2 to an array.

To avoid confusion in your module namespace, the classic decorator notation *@array2trace*  is not recommended outside the neurophysiotools module. However, *array2trace* can be imported tidily as follow:
    
    #imports and variables defined above
    # standard numpy functions can take **Trace** objects as input but will return *numpy.ndarray* as outputs
    numpy_convolved = np.convolve(trace, [1,1])
    type(numpy_convolved)
    <class 'neurophysiotools.Trace'>
    
    # using *nt.array2trace* one can avoid namespace issues and still make numpy functions return *Trace* objects
    convolve_trace=nt.array2trace(np.convolve)
    convolved = convolve_trace(trace, [1,1])
    type(convolved)
    <class 'neurophysiotools.Trace'>

**NB:** Since **Trace** is a subclass of *numpy.ndarray* the use of *nt.array2trace* is not required for numpy's universal funnctions *(ufunc)*. Please read numpy's documentation: [numpy ufunc](https://numpy.org/doc/stable/reference/ufuncs.html)  
    

### Functions
Functions divided by use:

#### Filter functions
Detailed descriptions are provided in the **help**, a working understanding of filters is expected.Please read around the difference between Butterworth and Bessel filters, briefly, the former has a steeper cutoff than the latter and neither causes ripples in standard conditions. The filters provided assume even sampling, which is a fair assumption with most high-end scientific instruments but may not be the case with hobbyist or consumer electronics. 
- bandpass_filter : a 4th Butterworth bandpass filter, frequencies outside the indicated values will be attenuated.  
- notch_filter : a 4th Butterworth 'notch' filter, meaning that all passes except a small window of frequncies. 
- butter_lowpass_filter : a 4th Butterworth lowpass filter, frequencies higher than the indicated values will be attenuated.
- butter_highpass_filter : a 4th Butterworth highpass filter, frequencies lower than the indicated values will be attenuated.
- bessel_lowpass_filter : a 4th Bessel lowpass filter, frequencies higher than the indicated values will be attenuated.
- bessel_highpass_filter : a 4th Bessel highpass filter, frequencies lower than the indicated values will be attenuated.
- running_mean: a running mean of the signal with a given window, similar to Matlab's popular smooth function.

#### Downsample functions
Simply reduce the size of the array by averaging the sampling points, detailed descriptions are provided in the **help**. If using the **Trace** object equivalent methods are available which will also adjust the sampling rate and time axis accordingly. The downsampling functions assume even sampling, which is a fair assumption with most high-end scientific instruments but may not be the case with hobbyist or consumer electronics. 
- downsample : downsamples the array of a given factor.
- downsample_to : downsamples the array to a given length. 

#### Functions for general measures 
Detailed descriptions are provided in the **help**, also here even sampling is assumed.
- coastline : the sum of the absolute value of the difference between point, you can think of it as the amount of 'ink' needed to plot the trace 
- get_freq_band : takes a frequency value and returns the band it belongs in a classic EEG frequency band classification, allows custom band definitions
- welch_an : !requires *pandas*! returns a pandas dataframe of Welch's spectral analysis.
- plot_welch : ! requires *pandas* and *plotly*! returns a plotly figure of a spectrogram

#### Functions for finding events
Detailed descriptions are provided in the **help**, also here even sampling is assumed.
- find_closest_value : provides the sampled value closest to an arbitrary value provided.
- find_closest_index : provides the index of a sampled value closest to an arbitrary value provided.
- find_threshold_crossings : list of indeces where the signal crosses a given threshold.
- find_incipit : provides the value before the beginning of an above thresholf event. NB: this is not the last point below threshold but where trend starts.


#### Other functions
- batch_open: finds all the files in a folder with a given extension.



### Contributions

Pull requests and suggestions for new features are welcome. Please feel free to get in touch via github or email; contribution guidelines will be uploaded soon.

### License

NeurophysioTools is licensed under the MIT License. See LICENSE for more information.
