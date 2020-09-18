import itertools
import numpy as np
import operator
import sys
from iautils.conversion import chime_intensity
from frb_common import common_utils
from iautils import fitburst_utils

import matplotlib.pyplot as plt

plt.switch_backend("agg")
import scipy.signal as ss
import scipy.stats as sst
from scipy.optimize import curve_fit
from astropy.modeling import models, fitting, polynomial

#sys.path.append('/Users/Documents/chime/')

def waterfall(snapshot_data):

    # load in snapshot data and bad-channel mas.
    #snapshot_data = np.load(filenames)
    #bad_channel_mask = np.loadtxt(filenames[1], dtype=np.int)

    # extract dynamic spectrum and important time/frequency/DM info.
    intensity_array = snapshot_data["intensity"]
    nfreq, ntime = intensity_array.shape
    fpga0s = snapshot_data["fpga0s"]
    binning = snapshot_data["binning"]
    event_time_fpga = int(
        (snapshot_data["fpgastart"] + snapshot_data["fpgaend"]) / 2
    )
    # `binning // 2`: correct for odd bining, slightly off for even.
    time_fpga = np.arange(ntime, dtype=np.int64) + binning // 2
    time_fpga *= (
        common_utils.l0_num_frames_sample * common_utils.l0_upchan_factor * binning
    )
    time_fpga += snapshot_data["fpgastart"]

    freq = np.arange(nfreq, dtype=np.float64)
    freq *= -common_utils.channel_bandwidth_mhz
    freq += common_utils.fpga_freq0_mhz
    freq += 400.0 / 1024 / 2  # half the FPGA channel bandwidth
    freq = freq[::-1].copy()  # unpack_datafiles reverses the order.

    # compute weights of input spectrum (i.e. channels that should be used in fit)
    #weights = []

    #for ii in range(nfreq):
        # weights.append([current_weight] * ntime)
        #if ii in bad_channel_mask:
            #weights.append([False] * ntime)

        #else:
            #weights.append([True] * ntime)

    #weights = np.array(weights)

    return intensity_array, freq, time_fpga, event_time_fpga #,weights

def subband(intensity_array, sub_factor):
    # Load in .npz or .npy file to generate waterfall plot
    sub_factor = np.int(sub_factor)

    # Sub-band frequency channels
    arr = np.nanmean(
        intensity_array.reshape(-1, sub_factor, intensity_array.shape[1]), axis=1
    )
    return arr


def bandpasscorr_first(initrow):
    # Effective bandpass correction
    row = (initrow - np.mean(initrow)) / np.std(initrow)
    return row


def bandpasscorr_sec(arr):
    arr = [bandpasscorr_first(row) for row in arr]
    arr = np.asarray(arr)
    arr = np.nan_to_num(arr)
    return arr


def nonoise(arr, axis):
    # Take 1/4 of the time series on either end and calculate an approximate snr
    summ = arr.sum(axis=axis)
    maxloc = np.where(summ == np.max(summ))[0]
    if len(maxloc) > 1:
        maxloc = 128
    else:
        maxloc = maxloc
    lenn = len(summ) / 6
    xcisef = summ[: int(int(maxloc) - lenn)]
    xcisel = summ[int(int(maxloc) + lenn) :]
    no_signal = np.concatenate((xcisef, xcisel))
    std_noise = np.std(no_signal)
    signal = summ[int(int(maxloc) - (lenn + 1)) : int(int(maxloc) + (lenn - 1))]
    snr = np.log(np.max(summ) / std_noise)

    # Set all values below the calculated snr to 0
    for n, i in enumerate(summ):
        if i < (np.max(summ) / snr):
            summ[n] = 0.0

    # Create a masked array to isolate the signal
    summask = []

    for i in summ:
        if i > 0:
            i = 1
            summask.append(i)
        else:
            summask.append(i)

    # Group "components" together based on divisions in the mask
    component_class = [
        [i for i, value in it]
        for key, it in itertools.groupby(enumerate(summask), key=operator.itemgetter(1))
        if key != 0
    ]

    # Enhance signal with an alternate mask
    summask_2 = []
    summask_2.append(summ)

    masked_arr = arr * summask_2

    return summ, summask, masked_arr, component_class


def tseriesdict(arr, axis):
    secs = arr.shape[np.int(axis)]

    # Make dictionary of 1D timeseries' or spectra depending on axis
    # 1 key per sample.
    tsdict = {}

    # Find maximum values for each array in the dictionary.
    maxvals = []

    # Find the array element location of these maximum values.
    yvals = []

    # Make a linspace to plot the yvals against.
    xtent = arr.shape[np.int(axis)]
    xvals = np.linspace(0, np.int(xtent), num=np.int(xtent))
    slices = np.split(arr, secs, axis=np.int(axis))

    for i in range(secs):
        tsdict[i] = np.ravel(slices[i])
        maxvals.append(np.max(tsdict[i]))
        y = list(np.where(tsdict[i] == tsdict[i].max()))

        # Replace RFI/already ommitted channels that, as expected,
        # contain multiple maximum vals with 0
        for n, i in enumerate(y):
            if len(i) > 1:
                y[n] = 0
        yvals.append(np.int(y[0]))

    return tsdict, maxvals, yvals, xvals


def detect_peaks(arr, axis, xval, yval, freq, sub_factor=32):
    #Ensure axis dependent binning for 1D arrays of max's

    if axis == 1:
        value = yval
        binnum = int(np.max(value))
    if axis == 0:
        value = yval
        binnum = int(np.max(value))

    peak_ts = arr.sum(0)
    peak_spec = arr.sum(1)

    #max_bin_hist = np.histogram(value, bins=binnum)

    #Determine number of peaks found in histogram

    ts_peakloc = ss.find_peaks(peak_ts)[0]
    # if len(ts_peakloc) == 0:
    # sys.exit("No Hint of Burst Detected")
    ts_peak_prom = ss.peak_prominences(peak_ts, ts_peakloc)[0]
    ts_peak_widths = ss.peak_widths(peak_ts, ts_peakloc)[0] * 0.96 #ms

    # Normalize prominences or "grades"
    ts_peak_grade = ts_peak_prom / np.max(ts_peak_prom)
    #peaks = np.where(max_bin_hist[0] > 0)
    #exclude peak in the first bin that is due to excised rfi channels
    #peaks = peaks[0][1:]

    freq_width = []
    spec_peaks = []
    for i in range(len(ts_peakloc)):
        #Fit a gaussian to determine spectral widths and other properties
        half_width_peak = np.int(np.int(ts_peak_widths[i]+1) / 2)
        low = np.int(ts_peakloc[i] - half_width_peak)
        upp = np.int(ts_peakloc[i] + half_width_peak)
        y = arr[:, low:upp].sum(1)
        x = np.linspace(0,len(y),len(y))
        a = np.max(y)
        m = np.mean(y)
        s = np.std(y)
        g_init = models.Gaussian1D(amplitude = a, mean = m, stddev = s)
        fit_g = fitting.LevMarLSQFitter()
        g = fit_g(g_init, x, y)
        gauss_fit = g(x)
        gauss_max = (np.max(g(x)))
        gauss_std = (np.std(g(x)))
        idx_max = (np.abs(g(x) - gauss_max)).argmin()
        freq_sub = freq.reshape(-1, sub_factor).mean(axis=1)
        #Append gaussian peaks for determining centroids
        y_values = []
        #If the fit is unsatisfactory, simply use averaged maxima in time bins to determine centroids
        if idx_max == 0:
            ts_peakloc_andwidth = list(range(low, upp+1))
            if len(ts_peakloc_andwidth) == 0:
                spec_peaks.append(yval[int(ts_peakloc[i])])
            else:
                for j in ts_peakloc_andwidth:
                    y_values.append(yval[int(j)])
                y_values = [y for y in y_values if y > 0]
                if len(y_values) > 0:
                    spec_peaks.append(int(np.mean(y_values)))
                else:
                    spec_peaks.append(0)
        else:
            spec_peaks.append(idx_max)
        print('spec_peaks', spec_peaks)

        freq_sub = freq.reshape(-1, sub_factor).mean(axis=1)

        #Again, if the fit is unsatisfactory, use maxima-determined centroid to obtain frequency widths

        if idx_max == 0:
            idx_std = (np.abs(g(x) - gauss_std)).argmin()
            new_max_idx = spec_peaks[i]
            idx_std = (np.abs(g(x) - gauss_std)).argmin()
            freq_width_l = np.abs(freq_sub[idx_std] - freq[new_max_idx])
            #freq_width.append(2 * freq_width_l)
            freq_width_r = []
            if new_max_idx >= (idx_std - new_max_idx):
                freq_width_r.append(freq_width_l)
            else:
                freq_width_r.append(freq_sub[new_max_idx] - 400)
            if (freq_width_l + freq_width_r[0]) > 400: #MHz
                freq_width.append(400)
            else:
                freq_width.append(freq_width_l + freq_width_r[0])
        else:
            idx_std = (np.abs(g(x) - gauss_std)).argmin()
            freq_width_l = np.abs(freq_sub[idx_std] - freq[idx_max])
            #freq_width.append(2 * freq_width_l)
            freq_width_r = []
            if idx_max >= (idx_std - idx_max):
                freq_width_r.append(freq_width_l)
            else:
                freq_width_r.append(freq_sub[idx_max] - 400)
            if (freq_width_l + freq_width_r[0]) > 400: #MHz
                freq_width.append(400)
            else:
                freq_width.append(freq_width_l + freq_width_r[0])

    #SIDELOBE DETECTION CHECKER

    # ACF Analysis
    peak_normal = peak_spec / np.max(peak_spec)
    peak_norm = peak_spec / np.max(peak_spec)
    peak_corr = ss.correlate(peak_spec, peak_spec)
    peak_corr_norm = np.asarray(peak_corr / np.max(peak_corr))[len(peak_norm) - 1 :]

    # Find peaks in ACF and their prominences
    peaks_in_corr = ss.find_peaks(peak_corr_norm)[0]
    peaks_prom = ss.peak_prominences(peak_corr_norm, peaks_in_corr)[0]
    peaks_prom_sorted = np.sort(peaks_prom)[::-1]
    peaks_prom_sorted_norm = peaks_prom_sorted / np.max(peaks_prom_sorted)
    #('peaks prom sorted norm', peaks_prom_sorted_norm)

    # Find the amplitude threshold and appropriate lag time attributed to
    # sidelobe spectral periodicity
    maxindex = np.where(peaks_prom == peaks_prom_sorted[0])[0]
    secmaxindex = np.where(peaks_prom == peaks_prom_sorted[2])[0]
    peak_thresh = peaks_prom[np.int(secmaxindex[0])] / peaks_prom[np.int(maxindex[0])]

    # Remove values that fall below the amplitude threshold--determined
    # by the ACF--from the normalized spectrum
    for n, i in enumerate(peak_norm):
        if i < peak_thresh:
            peak_norm[n] = 0

    # Conduct a second peak search, this time to isolate the peaks of
    # highest prominence so we derive a periodicity, characteristic of a
    # sidelobe detection
    peaks_in_spec = ss.find_peaks(peak_norm)[0]
    peaks_mask = ss.peak_prominences(peak_norm, peaks_in_spec)[0]

    # Determine the error in the period based on the widths of the
    # highest prominence peaks --the ones we're using
    width_errors = np.max(ss.peak_widths(peak_corr_norm, peaks_in_corr)[0]) / 2

    # Create a mask to isolate the peaks of highest prominence
    for n, i in enumerate(peaks_mask):
        if i < peak_thresh:
            peaks_mask[n] = 0
        else:
            peaks_mask[n] = 1

    # Apply the mask and remove masked values
    peaks_in_spec = peaks_in_spec * peaks_mask
    peaks_in_spec = [value for value in peaks_in_spec if value != 0]

    # Calculate the period -- if an outlier peak -- > 1 sigma -- is found
    # by calculating the zscore it's removed
    peak_diff = np.diff(peaks_in_spec)


    if len(peak_diff) > 1:
        peak_diff = [i for i in peak_diff if i >= 16] #MHz ripple
        zscore = np.abs(sst.zscore(peak_diff))
        z_greatersigma = (list(np.where(zscore > 1)[0]))[::-1]
        [peak_diff.pop(np.int(i)) for i in z_greatersigma]
        peak_diff = list(peak_diff)
        peak_period = np.mean(peak_diff)
    elif len(peak_diff) == 1:
        peak_period = peak_diff[0]
    else:
        peak_period = 0

    # Array index for peak of highest prominence in ACF
    spec_corr_maxindex = peaks_in_corr[np.int(maxindex[0])]

    #Calculate S/N for ACF to determine periodic validity
    signal_corr = peaks_prom[np.int(maxindex[0])]
    noise_base = np.std(peak_corr_norm[np.int(-(len(peak_corr_norm)/4)):])
    sidelobe_snr = np.log((signal_corr)/(noise_base))

    return ts_peakloc, ts_peak_grade, ts_peak_widths, peak_corr_norm, peak_normal, \
        peak_period, spec_corr_maxindex, width_errors, sidelobe_snr, freq_width, spec_peaks

def find_centroids(yvals_ax1, peakloc, freq, time_fpga, spec_peaks, sub_factor=32):

    centroid_coord_x_units = []
    centroid_coord_x_array = []
    for i in peakloc:
        centroid_coord_x_units.append(time_fpga[i])
        centroid_coord_x_array.append(i)
    #centroid_coord_x = [int(i) for i in centroid_coord_x]

    #centroid_coord_y_arr = []

    #Find correspoinding maxima between the axis 1 sweep and axis 0
    #sweep of original 2D array
    #for i in range(len(peakloc)):
        #centroid_coord_y_arr.append(yvals_ax1[int(peakloc[i])])

    centroid_coord_y_units = []
    centroid_coord_y_array = []
    freq_sub = freq.reshape(-1, sub_factor).mean(axis=1)
    for i in spec_peaks:
        centroid_coord_y_units.append(freq_sub[i])
        centroid_coord_y_array.append(i)

    return centroid_coord_x_units, centroid_coord_x_array, centroid_coord_y_units, centroid_coord_y_array


def find_subbursts(snapshot_data, eventid, sub_factor=32, plot_dyn=False,
                   plot_mask=False, src="./"):
    """Locate and classify sub-bursts in FRBs and identify sidelobe
    detections.

    Sub-bursts in multi-component FRBs are identified in part by
    time-series (ts) and frequency spectrum (fs) characteristics,
    including peaks, peak prominences, and an effective signal-to-noise
    threshold. Sub-bursts are also identified by the points in the
    intensity array at which 'sample maxima' (maxima calculated per
    sample in frequency and time) cluster or overlap. These are points
    of greatest intensity. (Disclaimer: no unsupervised ML clustering
    algorithm is currently being used). These two indicators combine to
    locate and classify sub-bursts. It should be mentioned that
    burst-components that do not exhibit a 'nulling' or distinct
    separation from their neighboring sub-bursts may be assigned
    multiple centroids and classified as a single sub-burst.

    Sidelobe detections are identified using a cautious spectral
    autocorrelation method to detect meaningful periodicities in a
    noisey frequency spectrum containing a signal of varying amplitudes,
    which are characteristic of sidelobe detections.

    Parameters
    ----------
    intensity_array : array_like
        2D array of intensity data.
    sub_factor : int, optional
        The factor by which the array of intensity data will be
        sub-banded in frequency. 32 by default.
    plot_dyn : bool, optional
        Plots dynamic spectrum with overlayed centroid scatter plot.
        False by default.
    plot_mask: bool, optional
        Plots masked dynamic spectrum. False by default.

    Returns
    -------
    centroid_idx : list of tuples
        x and y indices of the centroids in the dynamic spectrum.
    centroid_fpga : list of tuples
        FPGA timestamps and frequencies of the centroids in the dynamic
        spectrum.
    centroid_prominence : list
        List of normalized centroid prominences, or 'grades' (see
        scipy.signal.peak_prominences).
    centroid_widths : list
        List of pulse widths attributed to each centroid (see
        scipy.signal.peak_widths).
    centroid_bandwidths : list
        List of pulse bandwidths attributed to each centroid (see
        scipy.signal.peak_widhts).
    component_number : int
        Number of separate components, or sub-bursts, in the signal.
    drift_rate : float, bool
        Slope of linear least squares fit to centroids. False if only one
        centroid.
    sidelobe_detection : bool, float
        If dynamic spectrum contains a 'knotty' structure attributable to
        a slidelobe detection, returns True alongside S/N likelihood value.

    """
    intensity_array, freq, time_fpga, event_time_fpga = waterfall(
        snapshot_data)

    arr = subband(intensity_array, sub_factor)

    arr = bandpasscorr_sec(arr)

    summ, summask, masked_arr, component_class = nonoise(arr, 0)

    tsdict1, maxvals1, yvals1, xvals1=tseriesdict(masked_arr, 1)
    ts_peakloc1, ts_peak_grade1, ts_widths1, peak_corr_norm1, peak_normal1, \
        peak_period1, spec_corr_maxindex1, width_errors1, sidelobe_snr1, \
        freq_width1, spec_peaks1 = detect_peaks(masked_arr, 1, xvals1, yvals1, \
                                                freq)

    tsdict0, maxvals0, yvals0, xvals0=tseriesdict(masked_arr, 0)
    ts_peakloc0, ts_peak_grade0, ts_widths0, peak_corr_norm0, peak_normal0, \
        peak_period0, spec_corr_maxindex0, width_errors0, sidelobe_snr0, \
        freq_width0, spec_peaks0 = detect_peaks(masked_arr, 0, xvals0, yvals0, \
                                                freq)

    centroid_coord_x_units, centroid_coord_x_array, centroid_coord_y_units, \
        centroid_coord_y_array = find_centroids(yvals1, ts_peakloc0, freq, \
                                                time_fpga, spec_peaks1)

    centroid_location_time = list(ts_peakloc0)
    centroid_prominence = list(ts_peak_grade0)
    centroid_widths = list(ts_widths0)
    component_number = len(component_class)
    centroid_bandwidths = freq_width1

    #Sidelobe Detection Monitor
    if ((spec_corr_maxindex0 + width_errors0) >= peak_period0 \
        >= (spec_corr_maxindex0 - width_errors0)):
        sidelobe_detection = (True, sidelobe_snr1)
    else:
        sidelobe_detection = False

    # Calculate Drift Rate
    if len(centroid_coord_x_array) > 1:
        p1 = np.polyfit(centroid_coord_x_array, centroid_coord_y_array, 1)
        p3 = np.poly1d(p1)
        # xp = np.linspace((np.min(x)-2), (np.max(x)+2), 100)
        # plt.plot(xp, p3(xp), c='k', linewidth=4)
        centroid_drift_rate = p1[0]
    else:
        centroid_drift_rate = False

    if plot_dyn:
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        plt.imshow(arr, aspect="auto")
        plt.title("Dynamic Spectrum [Containing Centroids]")
        plt.scatter(centroid_coord_x_array, centroid_coord_y_array, c="r")
        plt.ylabel("Frequency (MHz)")
        plt.xlabel("Time (ms)")
        y = np.arange(0, len(masked_arr.sum(1)), 40)
        plt.yticks(y, np.arange(400, 800, np.int(400 / len(y))))
        plt.gca().invert_yaxis()
        #parse_filename= (filepath.split('_')[1]).split('/')[0]
        fig.savefig(src + "/DS_Centroids_{}.png".format(str(eventid)), dpi=150,
                    bbox_inches="tight")

    if plot_mask:
        fig2 = plt.figure()
        ax2 = fig2.add_subplot(111)
        plt.imshow(masked_arr, aspect="auto")
        plt.title("Masked Dynamic Spectrum")
        plt.ylabel("Frequency (MHz)")
        plt.xlabel("Time (ms)")
        y = np.arange(0, len(masked_arr.sum(1)), 40)
        plt.yticks(y, np.arange(400, 800, np.int(400 / len(y))))
        plt.gca().invert_yaxis()
        # parse_filename = (filepath.split('_')[1]).split('/')[0]
        fig2.savefig(src + "/Masked_DS_{}.png".format(str(eventid)), dpi=150,
                     bbox_inches="tight")

    # grab the x and y coordinates together
    centroid_idx = np.vstack((centroid_coord_x_array,
                              centroid_coord_y_array)).T
    centroid_fpga = np.vstack((centroid_coord_x_units,
                               centroid_coord_y_units)).T

    return centroid_idx, centroid_fpga, centroid_prominence, centroid_widths, \
        centroid_bandwidths, component_number, drift_rate, sidelobe_detection


if __name__ == "__main__":
    #filenames = input("Enter npz file :")
    #waterfall = np.load(npz)
    #intensity_array = waterfall['intensity']
    intensity_array, freq, time_fpga,\
    event_time_fpga, centroid_coord_x_array, centroid_coord_x_units, \
    centroid_coord_y_array, centroid_coord_y_units, \
    centroid_prominence, centroid_widths, component_number,\
    centroid_drift_rate, sidelobe_detection, centroid_bandwidths = find_subbursts(snapshot_data, eventid)
    #eventid = npz.split("_")[0]
    print(freq, time_fpga,\
        event_time_fpga, centroid_coord_x_array, centroid_coord_x_units, \
        centroid_coord_y_array, centroid_coord_y_units, \
        centroid_prominence, centroid_widths, component_number,\
        centroid_drift_rate, sidelobe_detection, centroid_bandwidths)
    #with open('sbclass_results.txt', 'a') as f:
        #f.write("%s\n" % results)
