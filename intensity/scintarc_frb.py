import itertools
import numpy as np
import operator
import sys

import matplotlib.pyplot as plt
plt.switch_backend('agg')
import scipy.signal as ss
import scipy.stats as sst
from numpy.fft import fftshift, fft2, ifft2
from matplotlib.colors import LogNorm


def subband(intensity_array, sub_factor):
    #Load in .npz or .npy file to generate waterfall plot
    sub_factor = np.int(sub_factor)

    #Sub-band frequency channels
    arr = np.nanmean(intensity_array.reshape(
        -1, sub_factor, intensity_array.shape[1]), axis=1)
    return arr


def bandpasscorr_first(initrow):
    #Effective bandpass correction
    row = (initrow - np.mean(initrow)) / np.std(initrow)
    return row


def bandpasscorr_sec(arr):
    arr = [bandpasscorr_first(row) for row in arr]
    arr = np.asarray(arr)
    arr = np.nan_to_num(arr)
    return arr


def nonoise(arr, axis):
    #Take 1/4 of the time series on either end and calculate an approximate snr
    summ = arr.sum(axis=axis)
    maxloc = np.where(summ == np.max(summ))[0]
    if len(maxloc) > 1:
        maxloc = 128
    else:
        maxloc = maxloc
    lenn = len(summ)/6
    xcisef = summ[:int(int(maxloc) - lenn)]
    xcisel = summ[int(int(maxloc) + lenn):]
    no_signal = np.concatenate((xcisef, xcisel))
    std_noise = np.std(no_signal)
    signal = summ[int(int(maxloc) - (lenn+1)):int(int(maxloc) + (lenn - 1))]
    snr = np.log(np.max(summ) / std_noise)

    #Set all values below the calculated snr to 0
    for n, i in enumerate(summ):
        if i < (np.max(summ) / snr):
            summ[n] = 0.0

    #Create a masked array to isolate the signal
    summask = []

    for i in summ:
        if i > 0:
            i = 1
            summask.append(i)
        else:
            summask.append(i)

    #Group "components" together based on divisions in the mask
    component_class = [[i for i, value in it] for key, it in \
        itertools.groupby(enumerate(summask), \
        key=operator.itemgetter(1)) if key != 0]

    #Enhance signal with an alternate mask
    summask_2 = []
    summask_2.append(summ)

    masked_arr = arr * summask_2

    return summ, summask, masked_arr, component_class


def gen_sec(intensity_array, eventid, sub_factor=32, plot_sec=False):

    fig = plt.figure(figsize = (35, 20))

    arr = subband(intensity_array, sub_factor)
    
    arr = bandpasscorr_sec(arr)
    
    summ, summask, masked_arr, component_class = nonoise(arr, 0)
    
    #Make a secondary spectrum
    hamm = np.hamming(len(arr[0]))
    arr_tap = arr * hamm
    ss = fftshift(fft2(arr_tap))
    ss = ss*np.conjugate(ss)
    ss_full = ss.real
    ss_crop = ss_full[256:306, 60:195]
    #ss = ss.real

    ax1 = fig.add_subplot(231)

    plt.imshow(arr, aspect = 'auto')
    plt.gca().invert_yaxis()
    plt.title(str(eventid) + " Dynamic Spectrum")
    plt.xlabel("Time (ms)")
    plt.ylabel("Frequency (MHz)")
    y = np.arange(0, len(arr.sum(1)), 64)
    x = np.arange(0, len(arr.sum(0)), 32)
    plt.yticks(y, np.arange(400, 800, 50))
    plt.xticks(x, np.arange(0, 256, 32))

    ax2 = fig.add_subplot(232)

    plt.plot(arr.sum(0))
    plt.title(str(eventid) + " Timeseries")
    plt.xlabel("Time (ms)")
    x = np.arange(0, len(arr.sum(0)), 32)
    plt.xticks(x, np.arange(0, 256, 32))

    ax3 = fig.add_subplot(233)

    plt.plot(arr.sum(1))
    plt.title(str(eventid) + " Spectrum")
    plt.xlabel("Frequency (MHz)")
    x = np.arange(0, len(arr.sum(1)), 64)
    plt.xticks(x, np.arange(400, 800, 50))

    ax4 = fig.add_subplot(234)

    #plt.plot(ss.sum(0))
    plt.imshow(ss_full, norm = LogNorm(), aspect = 'auto')
    plt.gca().invert_yaxis()
    plt.title(str(eventid)+ " " + "Secondary Spectrum")
    plt.ylabel("Time Delay " + r"($\mu$s)")
    plt.xlabel("Fringe Frequency " + "(Hz)")
    y = np.arange(0, len(ss.sum(1)), 64)
    x = np.arange(0, len(ss.sum(0)), 32)
    plt.yticks(y, np.arange(-2.5, 2.5, 0.625))
    plt.xticks(x, np.arange(-4, 4, 1))

    ax5 = fig.add_subplot(235)

    plt.imshow(ss_crop, norm = LogNorm(), aspect = 'auto')
    plt.gca().invert_yaxis()
    plt.title(str(eventid)+ " " + "Secondary Spectrum Cropped")
    plt.ylabel("Time Delay " + r"($\mu$s)")
    plt.xlabel("Fringe Frequency " + "(Hz)")
    y = np.arange(0, len(ss_crop.sum(1)), 6.2)
    x = np.arange(0, len(ss_crop.sum(0)), 11)
    plt.xticks(x, np.arange(-0.8, 0.8, 0.032))
    plt.yticks(y, np.arange(0, 0.128, 0.016))

    ax6 = fig.add_subplot(236)

    plt.imshow(ss_crop, aspect = 'auto')
    plt.gca().invert_yaxis()
    plt.title(str(eventid)+ " " + "Secondary Spectrum Cropped (No LogNorm)")
    plt.ylabel("Time Delay " + r"($\mu$s)")
    plt.xlabel("Fringe Frequency " + "(Hz)")
    y = np.arange(0, len(ss_crop.sum(1)), 6.2)
    x = np.arange(0, len(ss_crop.sum(0)), 11)
    plt.xticks(x, np.arange(-0.8, 0.8, 0.032))
    plt.yticks(y, np.arange(0, 0.128, 0.016))
    
    plt.tight_layout()
    fig.savefig(str(eventid) + " Secondary Spectrum")

    return
    
if __name__ == "__main__":
	
	#npz = input("Enter npz file :")
	#waterfall = np.load(npz)
	#intensity_array = waterfall['intensity']
	#eventid = npz.split("_")[0]
	gen_sec(intensity_array, eventid, plot_sec = True)
