# Perform FRB burst profile using nested sampling.

import re
import bilby
import csv

import math
import sys,os
import numpy as np
import pandas as pd

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import gridspec

from scipy import stats
import scipy.signal as ss

from lmfit.models import ExponentialGaussianModel, ExponentialModel, GaussianModel

from astropy import modeling
from astropy.modeling import models, fitting

from scipy.signal import convolve

from bokeh.models import ColumnDataSource, Div
from scipy.optimize import curve_fit

MIN_FLOAT = sys.float_info[3]

#Formatting
font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 14}

matplotlib.rc('font', **font)

def boxcar_kernel(width):
    """Returns the boxcar kernel of given width normalized by
    sqrt(width) for S/N reasons.
    Parameters
    ----------
    width : int
        Width of the boxcar.
    Returns
    -------
    boxcar : array_like
        Boxcar of width `width` normalized by sqrt(width).
    """
    width = int(round(width, 0))
    return np.ones(width, dtype="float32") / np.sqrt(width)


def find_burst(ts, width_factor=4, min_width=1, max_width=2048):
    """Find burst peak and width using boxcar convolution.

    Parameters
    ----------
    ts : array_like
        Time-series.
    width_factor : int, optional
        Windowing factor for on and off-pulse determination.
    min_width : int, optional
        Minimum width to search from, in number of time samples.
        1 by default.
    max_width : int, optional
        Maximum width to search up to, in number of time samples.
        128 by default.

    Returns
    -------
    peak : int
        Index of the peak of the burst in the time-series.
    width : int
        Width of the burst in number of samples.
    snr : float
        S/N of the burst.

    """
    min_width = int(min_width)
    max_width = int(max_width)

    # do not search widths bigger than timeseries
    widths = list(range(min_width,
                        min(max_width + 1, int((len(ts) - 50) // 6))))

    # envelope finding
    snrs = np.empty_like(widths, dtype=float)
    peaks = np.empty_like(widths, dtype=int)

    # borders for on and off-pulse determination
    outer = 3 * width_factor // 2
    inner = width_factor // 2

    for i in range(len(widths)):
        convolved = ss.convolve(ts, boxcar_kernel(widths[i]))
        peaks[i] = np.nanargmax(convolved)
        # peak should not be on the edge of time-series
        if (peaks[i] > 0.999 * ts.shape[0]) or (peaks[i] < 0.001 * ts.shape[0]):
            snrs[i] = np.nan
        else:
            # get RMS for S/N weighting, as in PRESTO's single_pulse_search.py
            baseline = np.concatenate(
                [
                    convolved[0 : max(0, peaks[i] - 3 * widths[i])],
                    convolved[peaks[i] + 3 * widths[i] :],
                ]
            )

            # cutoff of at least 50 samples is a bit arbitrary, but seems
            # reasonable
            if baseline.shape[0] > 50:
                rms = np.std(baseline)
            else:
                rms = np.nan

            snrs[i] = convolved[peaks[i]] / rms

    best_idx = np.nanargmax(snrs)

    return peaks[best_idx]-widths[best_idx]//2, widths[best_idx], snrs[best_idx]

def Gaussian1D(x,sig,x0):
    """
    Returns 1D Gaussian curve.
    """
    return np.exp(-(x-x0)*(x-x0)/(2*sig*sig + MIN_FLOAT))

def exp_decay(x,tau,x0):
    """
    Returns 1D one-sided exponential curve.
    """
    res = np.zeros(len(x)) + MIN_FLOAT
    res[x > x0] = np.exp(-(x[x>x0]-x0)/(tau+MIN_FLOAT))
    return res

def exp_gauss(x,x0,amp,sig,tau):
    """
    Returns Gaussian convolved with a one-sided exponential.
    """
    gx0 = np.mean(x)
    g = Gaussian1D(x,sig,gx0)
    ex = exp_decay(x,tau,x0)
    conv = convolve(g,ex,"same")
    conv /= np.max(conv) + MIN_FLOAT
    return amp*conv

  
def exp_gauss_1(x,x1,amp1,sig1,tau1):
    g1 = exp_gauss(x,x1,amp1,sig1,tau1)
    return g1

def exp_gauss_2(x,x1,amp1,sig1,tau1,
              x2,amp2,sig2,tau2):
    g1 = exp_gauss(x,x1,amp1,sig1,tau1)
    g2 = exp_gauss(x,x2,amp2,sig2,tau2)
    return g1 + g2

def exp_gauss_3(x,x1,amp1,sig1,tau1,
              x2,amp2,sig2,tau2,
              x3,amp3,sig3,tau3):
    g1 = exp_gauss(x,x1,amp1,sig1,tau1)
    g2 = exp_gauss(x,x2,amp2,sig2,tau2)
    g3 = exp_gauss(x,x3,amp3,sig3,tau3)
    return g1 + g2 + g3

def exp_gauss_4(x,x1,amp1,sig1,tau1,
              x2,amp2,sig2,tau2,
              x3,amp3,sig3,tau3,
              x4,amp4,sig4,tau4):
    g1 = exp_gauss(x,x1,amp1,sig1,tau1)
    g2 = exp_gauss(x,x2,amp2,sig2,tau2)
    g3 = exp_gauss(x,x3,amp3,sig3,tau3)
    g4 = exp_gauss(x,x4,amp4,sig4,tau4)
    return g1 + g2 + g3 + g4


#class fit_function:
#    def __init__(self,
#                 x=None,
#                 x1=None,
#                 amp1=None,
#                 sig1=None,
#                 tau1=None,
#                 x2=None,
#                 amp2=None,
#                 sig2=None,
#                 tau2=None,
#                 x3=None,
#                 amp3=None,
#                 sig3=None,
#                 tau3=None,
#                 x4=None,
#                 amp4=None,
#                 sig4=None,
#                 tau4=None,
#                 comp_num=1):
#        self.x=x
#        self.x1=x1
#        self.amp1=amp1
#        self.sig1=sig1
#        self.tau1=tau1
#        self.x2=x2
#        self.amp2=amp2
#        self.sig2=sig2
#        self.tau2=tau2
#        self.x3=x3
#        self.amp3=amp3
#        self.sig3=sig3
#        self.tau3=tau3
#        self.x4=x4
#        self.amp4=amp4
#        self.sig4=sig4
#        self.tau4=tau4
#        self.comp_num = comp_num
#        self.curve = None
#        if self.comp_num == 1:
#            self.curve = self.exp_gauss_1()
#        if self.comp_num == 2:
#            self.curve = self.exp_gauss_2()
#        if self.comp_num == 3:
#            self.curve = self.exp_gauss_3()
#        if self.comp_num == 4:
#            self.curve = self.exp_gauss_4()
            
#def exp_gauss_1(self):
#    g1 = exp_gauss(self.x,self.x1,self.amp1,self.sig1,self.tau1)
#    return g1
#
#def exp_gauss_2(self):
#    g1 = exp_gauss(self.x,self.x1,self.amp1,self.sig1,self.tau1)
#    g2 = exp_gauss(self.x,self.x2,self.amp2,self.sig2,self.tau2)
#    return g1 + g2
#
#def exp_gauss_3(self):
#    g1 = exp_gauss(self.x,self.x1,self.amp1,self.sig1,self.tau1)
#    g2 = exp_gauss(self.x,self.x2,self.amp2,self.sig2,self.tau2)
#    g3 = exp_gauss(self.x,self.x3,self.amp3,self.sig3,self.tau3)
#    return g1 + g2 + g3
#
#def exp_gauss_4(self):
#    g1 = exp_gauss(self.x,self.x1,self.amp1,self.sig1,self.tau1)
#    g2 = exp_gauss(self.x,self.x2,self.amp2,self.sig2,self.tau2)
#    g3 = exp_gauss(self.x,self.x3,self.amp3,self.sig3,self.tau3)
#    g4 = exp_gauss(self.x,self.x4,self.amp4,self.sig4,self.tau4)
#    return g1 + g2 + g3 + g4

def sub_npy(npy_fil, subfactor):
    """
    Returns original numpy array, a sub-banded numpy array scrunched by a subfactor,
    and the array timeseries.
    """
    npy = np.load(npy_fil)
    #npy = npy[30:, :]
    npy_fsub = np.flipud(np.nanmean(npy.reshape(-1, int(subfactor), int(npy.shape[1])), axis=1))
    #npy_tsub = np.nanmean(npy_fsub.reshape(-1, int(npy_fsub.shape[0]), int(npy_fsub.shape[1]//time_subfactor)))
    timeseries = npy_fsub.sum(0)
    return npy, npy_fsub, timeseries

def nested_sampling(npy_fil, p0, comp_num, nlive=1000, bandwidth=400., center_frequency=800., file_duration=81.39, subfactor=1., debug=False):
    """
    Perform nested sampling profile fitting with bilby.
    """
    npy, npy_sub, timeseries = sub_npy(npy_fil, subfactor)
    peaks, widths, snrs = find_burst(timeseries)

    # Calculate the time and frequency resolutions of the array
    time_res = file_duration / npy.shape[1]
    print('Raw Time Resolution (microsec): ', time_res*1e3)
    num_chan = npy.shape[0]
    freq_res = bandwidth / num_chan
    print('Raw Frequency Resolution (kHz): ', freq_res*1e3)

    #Define windowing depending on where burst sits in dynspec
    window_left = int(peaks - 1*widths)
    window_right = int(peaks + 1*widths)
    
    if window_right - window_left <= 100:
        window_right = window_right + 20
        window_left = window_left - 20
    if window_left <= 0:
        window_left = 0
    if window_right >= npy.shape[1]:
        window_right = npy.shape[1]

    print('Window Left: ', window_left)
    print('Window Right: ', window_right)
    sub_factor_time = 1
    y_data = (npy[:, :].sum(0)/np.max(npy[:, :].sum(0)))[window_left:window_right]
    print('Lenght of Data: ', len(y_data))
    sampling_time = (file_duration / npy.shape[1]) * sub_factor_time
    print('Sampling Time (ms): ', sampling_time)
    time = np.arange(len(y_data)) * sampling_time
    sigma = np.repeat(sampling_time, len(time))
    print('Time: ', time)
    print('Sigma: ', sigma)

    if debug:
        # Random Exponential Gaussian to test fit
        fig = plt.figure()
        plt.errorbar(time, y_data, yerr=sigma)
        plt.show()

    print('Fitting Initiated')
    
    label = str(npy_file.split('/')[-1].split('_3')[0])
    outdir = str(npy_file.split('/')[-1].split('_3')[0]) + '_profilefit'
    
    if comp_num == '1':
        #Upper and lower bounds for prior
        lower_bounds = [(i - i) for i in p0[0:4]]
        upper_bounds = [(i + i) for i in p0[0:4]]
        lower_bounds = [round(i, 2) for i in lower_bounds]
        upper_bounds = [round(i, 2) + 1. for i in upper_bounds]
        #upper_bounds = list(map(lambda i: 1 if i==0 else i, upper_bounds))

        print('Lower: ', lower_bounds)
        print('Upper: ', upper_bounds)

        likeli = bilby.core.likelihood.GaussianLikelihood(time, y_data, exp_gauss_1, sigma = sigma)
        prior = dict(x1 = bilby.core.prior.Uniform(lower_bounds[0], upper_bounds[0],'x1'),
                 amp1 = bilby.core.prior.Uniform(lower_bounds[1], upper_bounds[1],'amp1'),
                 sig1 = bilby.core.prior.Uniform(lower_bounds[2], upper_bounds[2],'sig1'),
                 tau1 = bilby.core.prior.Uniform(lower_bounds[3], upper_bounds[3],'tau1'))
        injection_params = dict(x1=p0[0],
                                    amp1=p0[1],
                                    sig1=p0[2],
                                    tau1=p0[3])
        print('Sampler Running')

        # Do the nested sampling
        result = bilby.run_sampler(likelihood=likeli, priors=prior, injection_parameters=injection_params, \
                                sampler='dynesty', nlive=nlive, outdir=outdir, label=label)

        result.plot_corner()

        print('Fit Complete!')

    elif comp_num == '2':
        #Upper and lower bounds for prior
        lower_bounds = [(i - i/2) for i in p0[0:8]]
        upper_bounds = [(i + i/2) for i in p0[0:8]]
        lower_bounds = [round(i, 2) for i in lower_bounds]
        upper_bounds = [round(i, 2) + 1. for i in upper_bounds]
        #upper_bounds = list(map(lambda i: 1 if i==0 else i, upper_bounds))

        print('Lower: ', lower_bounds)
        print('Upper: ', upper_bounds)
        likeli = bilby.core.likelihood.GaussianLikelihood(time, y_data, exp_gauss_2, sigma = sigma)
        prior = dict(x1 = bilby.core.prior.Uniform(lower_bounds[0], upper_bounds[0],'x1'),
                 amp1 = bilby.core.prior.Uniform(lower_bounds[1], upper_bounds[1],'amp1'),
                 sig1 = bilby.core.prior.Uniform(lower_bounds[2], upper_bounds[2],'sig1'),
                 tau1 = bilby.core.prior.Uniform(lower_bounds[3], upper_bounds[3],'tau1'),
                   x2 = bilby.core.prior.Uniform(lower_bounds[4], upper_bounds[4],'x2'),
                 amp2 = bilby.core.prior.Uniform(lower_bounds[5], upper_bounds[5],'amp2'),
                 sig2 = bilby.core.prior.Uniform(lower_bounds[6], upper_bounds[6],'sig2'),
                 tau2 = bilby.core.prior.Uniform(lower_bounds[7], upper_bounds[7],'tau2'))
        injection_params = dict(x1=p0[0],
                                    amp1=p0[1],
                                    sig1=p0[2],
                                    tau1=p0[3],                     
                                      x2=p0[4],
                                    amp2=p0[5],
                                    sig2=p0[6],
                                    tau2=p0[7])

        print('Sampler Running')

        # Do the nested sampling
        result = bilby.run_sampler(likelihood=likeli, priors=prior, injection_parameters=injection_params, \
                                sampler='dynesty', nlive=nlive, outdir=outdir, label=label)

        result.plot_corner()

        print('Fit Complete!')

    elif comp_num == '3':
        #Upper and lower bounds for prior
        lower_bounds = [(i - i) for i in p0[0:12]]
        upper_bounds = [(i + i) for i in p0[0:12]]
        lower_bounds = [round(i, 2) for i in lower_bounds]
        upper_bounds = [round(i, 2) + 1. for i in upper_bounds]
        #upper_bounds = list(map(lambda i: 1 if i==0 else i, upper_bounds))

        print('Lower: ', lower_bounds)
        print('Upper: ', upper_bounds)
        likeli = bilby.core.likelihood.GaussianLikelihood(time, y_data, exp_gauss_3, sigma = sigma)
        prior = dict(x1 = bilby.core.prior.Uniform(lower_bounds[0], upper_bounds[0],'x1'),
                 amp1 = bilby.core.prior.Uniform(lower_bounds[1], upper_bounds[1],'amp1'),
                 sig1 = bilby.core.prior.Uniform(lower_bounds[2], upper_bounds[2],'sig1'),
                 tau1 = bilby.core.prior.Uniform(lower_bounds[3], upper_bounds[3],'tau1'),
                   x2 = bilby.core.prior.Uniform(lower_bounds[4], upper_bounds[4],'x2'),
                 amp2 = bilby.core.prior.Uniform(lower_bounds[5], upper_bounds[5],'amp2'),
                 sig2 = bilby.core.prior.Uniform(lower_bounds[6], upper_bounds[6],'sig2'),
                 tau2 = bilby.core.prior.Uniform(lower_bounds[7], upper_bounds[7],'tau2'),
                   x3 = bilby.core.prior.Uniform(lower_bounds[8], upper_bounds[8], 'x3'),
                 amp3 = bilby.core.prior.Uniform(lower_bounds[9], upper_bounds[9],'amp3'),
                 sig3 = bilby.core.prior.Uniform(lower_bounds[10], upper_bounds[10],'sig3'),
                 tau3 = bilby.core.prior.Uniform(lower_bounds[11], upper_bounds[11],'tau3'))
        injection_params = dict(x1=p0[0],
                                    amp1=p0[1],
                                    sig1=p0[2],
                                    tau1=p0[3],                     
                                      x2=p0[4],
                                    amp2=p0[5],
                                    sig2=p0[6],
                                    tau2=p0[7],                        
                                      x3=p0[8],
                                    amp3=p0[9],
                                    sig3=p0[10],
                                    tau3=p0[11])

        print('Sampler Running')

        # Do the nested sampling
        result = bilby.run_sampler(likelihood=likeli, priors=prior, injection_parameters=injection_params, \
                                sampler='dynesty', nlive=nlive, outdir=outdir, label=label)

        result.plot_corner()

        print('Fit Complete!')

    elif comp_num == '4':
        #Upper and lower bounds for prior
        lower_bounds = [(i - i/2) for i in p0[0:16]]
        upper_bounds = [(i + i/2) for i in p0[0:16]]
        lower_bounds = [round(i, 2) for i in lower_bounds]
        upper_bounds = [round(i, 2) + 1. for i in upper_bounds]
        #upper_bounds = list(map(lambda i: 1 if i==0 else i, upper_bounds))

        print('Lower: ', lower_bounds)
        print('Upper: ', upper_bounds)
        likeli = bilby.core.likelihood.GaussianLikelihood(time, y_data, exp_gauss_4, sigma = sigma)
        prior = dict(x1 = bilby.core.prior.Uniform(lower_bounds[0], upper_bounds[0],'x1'),
                 amp1 = bilby.core.prior.Uniform(lower_bounds[1], upper_bounds[1],'amp1'),
                 sig1 = bilby.core.prior.Uniform(lower_bounds[2], upper_bounds[2],'sig1'),
                 tau1 = bilby.core.prior.Uniform(lower_bounds[3], upper_bounds[3],'tau1'),
                   x2 = bilby.core.prior.Uniform(lower_bounds[4], upper_bounds[4],'x2'),
                 amp2 = bilby.core.prior.Uniform(lower_bounds[5], upper_bounds[5],'amp2'),
                 sig2 = bilby.core.prior.Uniform(lower_bounds[6], upper_bounds[6],'sig2'),
                 tau2 = bilby.core.prior.Uniform(lower_bounds[7], upper_bounds[7],'tau2'),
                   x3 = bilby.core.prior.Uniform(lower_bounds[8], upper_bounds[8], 'x3'),
                 amp3 = bilby.core.prior.Uniform(lower_bounds[9], upper_bounds[9],'amp3'),
                 sig3 = bilby.core.prior.Uniform(lower_bounds[10], upper_bounds[10],'sig3'),
                 tau3 = bilby.core.prior.Uniform(lower_bounds[11], upper_bounds[11],'tau3'),
                   x4 = bilby.core.prior.Uniform(lower_bounds[12], upper_bounds[12],'x4'),
                 amp4 = bilby.core.prior.Uniform(lower_bounds[13], upper_bounds[13],'amp4'),
                 sig4 = bilby.core.prior.Uniform(lower_bounds[14], upper_bounds[14],'sig4'),
                 tau4 = bilby.core.prior.Uniform(lower_bounds[15], upper_bounds[15],'tau4'))
        injection_params = dict(x1=p0[0],
                                    amp1=p0[1],
                                    sig1=p0[2],
                                    tau1=p0[3],                     
                                      x2=p0[4],
                                    amp2=p0[5],
                                    sig2=p0[6],
                                    tau2=p0[7],                        
                                      x3=p0[8],
                                    amp3=p0[9],
                                    sig3=p0[10],
                                    tau3=p0[11],
                                      x4=p0[12],
                                    amp4=p0[13],
                                    sig4=p0[14],
                                    tau4=p0[15])

        print('Sampler Running')

        # Do the nested sampling
        result = bilby.run_sampler(likelihood=likeli, priors=prior, injection_parameters=injection_params, \
                                sampler='dynesty', nlive=nlive, outdir=outdir, label=label)

        result.plot_corner()

        print('Fit Complete!')

    else:
        print('No Burst Component Number Specified -- Please Identify and Define Number of Sub-Bursts')
    

    print('Injection Parameters: ', injection_params)
        
    # Adjust initial conditions and priors depending the number of burst components
    #if comp_num == 1:
    #    injection_params = dict(list(injection_params.items())[0:4]) 
    #    prior = dict(list(prior.items())[0:4])  
    #elif comp_num == 2:
    #    injection_params = dict(list(injection_params.items())[0:8]) 
    #    prior = dict(list(prior.items())[0:8])  
    #elif comp_num == 3: 
    #    injection_params = dict(list(injection_params.items())[0:12]) 
    #    prior = dict(list(prior.items())[0:12]) 
    #elif comp_num == 4:
    #    injection_params = dict(list(injection_params.items())[0:16]) 
    #    prior = dict(list(prior.items())[0:16])

#if __name__ == '__main__':

npy_fil = sys.argv[1]
injection_csv = sys.argv[2]
comp_num = sys.argv[3]
file_dur = sys.argv[4]

cwd = os.getcwd() #get current working directory

npy_file = str(cwd) + '/' + str(npy_fil) #get npy array file name
print('Numpy File Path: ', npy_file)
 #get csv of initial conditions

# Read in csv containing initial conditions for each profile
with open(str(cwd) + '/' + str(injection_csv)) as csv_file:
    reader = csv.reader(csv_file)
    injection_dict = dict(reader)
    
p0_str = injection_dict[npy_file.split('/')[-1].split('_3')[0]]
p0 = [float(i) for i in p0_str.split('[')[1].split(']')[0].split(',')]
print('P0: ', p0)
print('P0 Shape: ', len(p0))

nested_sampling(npy_file, list(p0), str(comp_num), file_duration=float(file_dur))

        
