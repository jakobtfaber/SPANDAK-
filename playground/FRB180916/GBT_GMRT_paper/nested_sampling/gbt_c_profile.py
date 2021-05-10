import numpy as np
import matplotlib.pyplot as plt
import csv
from matplotlib import gridspec
from scipy import stats
import pandas as pd
import scipy.signal as ss
import re
from burst_utils import find_burst, boxcar_kernel
from astropy import modeling
from astropy.modeling import models, fitting
import matplotlib
from lmfit.models import ExponentialGaussianModel, ExponentialModel, GaussianModel
import math
from scipy.signal import convolve
import sys,os
from bokeh.models import ColumnDataSource, Div
MIN_FLOAT = sys.float_info[3]
from scipy.optimize import curve_fit
import bilby

#Formatting
font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 14}

matplotlib.rc('font', **font)

def Gaussian1D(x,sig,x0):
    return np.exp(-(x-x0)*(x-x0)/(2*sig*sig + MIN_FLOAT))

def exp_decay(x,tau,x0):
    res = np.zeros(len(x)) + MIN_FLOAT
    res[x > x0] = np.exp(-(x[x>x0]-x0)/(tau+MIN_FLOAT))
    return res

def exp_gauss(x,x0,amp,sig,tau):
    gx0 = np.mean(x)
    g = Gaussian1D(x,sig,gx0)
    ex = exp_decay(x,tau,x0)
    conv = convolve(g,ex,"same")
    conv /= np.max(conv) + MIN_FLOAT
    return amp*conv

def sub_npy(npy_fil, subfactor, file_duration, bandwidth, center_frequency):
    npy = np.load(npy_fil)
    npy_sub = np.flipud(np.nanmean(npy.reshape(-1, subfactor, npy.shape[1]), axis=1))
    timeseries = npy_sub.sum(0)
    return npy, npy_sub, timeseries

npy_fil = 'C_1164_dm348.8_lores.npy'
bandwidth = 400. #MHz
center_frequency = 800. #MHz
file_duration = 333.33 #ms
subfactor = 1

npy, npy_sub, timeseries = sub_npy(npy_fil, subfactor, file_duration, bandwidth, center_frequency)
peaks, widths, snrs = find_burst(timeseries)
nchan = npy.shape[0]
freq_res = bandwidth / nchan


#Resolutions
tres = file_duration / npy.shape[1]
print('Raw Time Resolution (microsec): ', tres*1e3)
nchan = npy.shape[0]
fres = bandwidth / nchan
print('Raw Frequency Resolution (kHz): ', fres*1e3)

#Define windowing depending on where burst sits in dynspec
window_left = int(peaks - 1*widths)
window_right = int(peaks + 1*widths)

if window_right - window_left <= 100:
    window_right = window_right + 300
    window_left = window_left - 300

sub_factor_time = 1
y_data = (npy[:].sum(0)/np.max(npy[:].sum(0)))[window_left:window_right]
#y_data = y_data.reshape(-1, sub_factor_time).mean(axis=1)
sampling_time = (file_duration / npy.shape[1]) * sub_factor_time
print('Sampling Time (ms): ', sampling_time)
time = np.arange(len(y_data)) * sampling_time
sigma = np.repeat(sampling_time, len(time))

#Initial Guesses (fit with Wael's slider)
#GBTA
p0 = [60.8, 0.48, 2.2, 1.6]

#Upper and lower bounds for prior
lower_bounds = [(i - i/2) for i in p0]
upper_bounds = [(i + i/2) for i in p0]
lower_bounds = [round(i, 2) for i in lower_bounds]
upper_bounds = [round(i, 2) for i in upper_bounds]

print('Lower: ', lower_bounds)
print('Upper: ', upper_bounds)

#Define data to fit to
injection_params = dict(x0=p0[0],
                        amp=p0[1],
                        sig=p0[2],
                        tau=p0[3])

#Define data to fit to
#sampling_time = 0.08192
#time = np.arange(window_right-window_left) * sampling_time
#injection_params = dict(x0=2., amp=10., sig=0.8, tau=4.)
#exg = exp_gauss4(time, **injection_params)
#sigma = np.repeat(sampling_time, len(exg))

#Random Exponential Gaussian to test fit
fig = plt.figure()
plt.errorbar(time, y_data, yerr=sigma)
plt.show()

print('Fitting Initiated')

label = 'GBT_C'
outdir = 'frb_profiles'

likeli = bilby.core.likelihood.GaussianLikelihood(time, y_data, exp_gauss, sigma = sigma)

prior =   dict(x0 = bilby.core.prior.Uniform(lower_bounds[0], upper_bounds[0],'x1'),
             amp = bilby.core.prior.Uniform(lower_bounds[1], upper_bounds[1],'amp1'),
             sig = bilby.core.prior.Uniform(lower_bounds[2], upper_bounds[2],'sig1'),
             tau = bilby.core.prior.Uniform(lower_bounds[3], upper_bounds[3],'tau1'))
     
print('Sampler Running')

result = bilby.run_sampler(
    likelihood=likeli, priors=prior, injection_parameters=injection_params, sampler='dynesty', nlive=500, outdir=outdir, label=label)

result.plot_corner()

