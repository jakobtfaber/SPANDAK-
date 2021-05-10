import numpy as np
from scipy.signal import convolve
from burst_utils import find_burst, boxcar_kernel
import sys
MIN_FLOAT = sys.float_info[3]

#Fitting Functions

def get_mad(ts):
    return np.median(np.abs(ts - np.median(ts)))

def normalise(ts):
    return ts/(1.4826*get_mad(ts))


def Gaussian1D(x,sig,x0):
    return np.exp(-(x-x0)*(x-x0)/(2*sig*sig + MIN_FLOAT))

def linear(x,a,b):
    return a*x + b

def exp_decay(x,tau,x0):
    res = np.zeros(len(x)) + MIN_FLOAT
    #res[x <= x0] = MIN_FLOAT
    res[x > x0] = np.exp(-(x[x>x0]-x0)/(tau+MIN_FLOAT))
    return res

def exp_gauss(x,x0,amp,sig,tau,eps):
    gx0 = np.mean(x)
    g = Gaussian1D(x,sig,gx0)
    ex = exp_decay(x,tau,x0)
    conv = convolve(g,ex,"same")
    conv /= np.max(conv) + MIN_FLOAT
    return amp*conv + eps

def exp_gauss_4(x,x1,amp1,sig1,tau1,
                  x2,amp2,sig2,tau2,
                  x3,amp3,sig3,tau3,
                  x4,amp4,sig4,tau4):
    g1 = exp_gauss(x,x1,amp1,sig1,tau1,0)
    g2 = exp_gauss(x,x2,amp2,sig2,tau2,0)
    g3 = exp_gauss(x,x3,amp3,sig3,tau3,0)
    g4 = exp_gauss(x,x4,amp4,sig4,tau4,0)

    fwhm = 2.355*sig1 + 2.355*sig2 + 2.355*sig3 + 2.355*sig4
    
    if sig2 > 0.001 and sig3 > 0.001 and sig4 > 0.001:
        fwhm2 = (x4 + (2.355*sig4)/2) - (x1 - (2.355*sig1)/2)
    elif sig2 > 0.001 and sig3 > 0.001:
        fwhm2 = (x3 + (2.355*sig3)/2) - (x1 - (2.355*sig1)/2)
    elif sig2 > 0.001:
        fwhm2 = (x2 + (2.355*sig2)/2) - (x1 - (2.355*sig1)/2)
    else:
        fwhm2 = 2.355 * sig1
    
    print('Sigmas: ', sig1, sig2, sig3, sig4)
    g = g1+g2+g3+g4
    
    return g, fwhm, fwhm2

def lnlike(theta, x, y):
    p0 = theta
    model = exp_gauss_4(x,*theta)
#    inv_sig = 1./(model**2)
    #print('Sigma2: ', sigma2)
    chisqr = -0.5*np.sum((y-model)**2)
    return chisqr


def burst_info(idx):
    npy_fils = ["propagation_effects/A_117_dm348.8.fits.npy",
                "propagation_effects/B_686_dm348.8.fits.npy",
                "propagation_effects/C_1164_dm348.8_lores.npy",
                "propagation_effects/D_267_dm348.8_lores.npy",
                "propagation_effects/E_579_dm348.8.fits.npy",
                "propagation_effects/F_639_dm348.8.fits.npy",
                "propagation_effects/G_1549_dm348.8.fits.npy",
                "propagation_effects/GMRT_A.dynamicspec_348.8.npy",
                "propagation_effects/GMRT_B.dynamicspec_349.19.npy",
                "propagation_effects/GMRT_C.dynamicspec_350_scrunch20.npy",
                "propagation_effects/GMRT_D.dynamicspec_349.19_scrunch4.npy"]

    npy = np.load(npy_fils[idx])
    timeseries = npy.sum(0)

    keys = ['gbta', 'gbtb', 'gbtc', 'gbtd', 'gbte', 'gbtf', 'gbtg', 'gmrta', 'gmrtb', 'gmrtc', 'gmrtd']
    pfits_key = keys[idx]
    print('Burst: ', pfits_key)
    bw = [400., 400., 400., 400., 400., 400., 400., 200., 200., 200., 200.] #MHz
    cf = [800., 800., 800., 800., 800., 800., 800., 400., 400., 400., 400.] #MHz
    fd = [83.33, 83.33, 333.33, 166.68, 83.33, 83.33, 83.33, 163.84, 122.88, 400.0, 300.] #ms
    subfactor = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
    bandwidth = bw[idx]
    print('Bandwidth (MHz): ', bw[idx])
    center_frequency = cf[idx]
    print('Center Frequency (MHz): ', cf[idx])
    file_duration = fd[idx]
    print('File Duration (ms): ', fd[idx])
    sampling_time = (file_duration / npy.shape[1])
    print('Sampling Time (ms): ', sampling_time)
    print('Number of Samples: ', npy.shape[1])
    nchan = npy.shape[0]
    print('Number of Channels: ', nchan)
    
    #def sub_npy(npy_fil, file_duration, bandwidth, center_frequency, subfactor = 1):
    #    npy = np.load(npy_fil)
    #    npy_sub = np.flipud(np.nanmean(npy.reshape(-1, subfactor, npy.shape[1]), axis=1))
    #    timeseries = npy_sub.sum(0)
        #return npy, npy_sub, timeseries
    #npy, npy_sub, timeseries = sub_npy(npy, file_duration, bandwidth, center_frequency)
    
    #Find burst peaks, widths & snrs
    peaks, widths, snrs = find_burst(timeseries)
    print('Peak Location (grid): ', peaks)
    print('Peak Width (grid): ', widths)
    
    #Resolutions
    tres = file_duration / npy.shape[1]
    print('Raw Time Resolution (microsec): ', tres*1e3)
    fres = bandwidth / nchan
    print('Raw Frequency Resolution (kHz): ', fres*1e3)
    
    #Define windowing depending on where burst sits in dynspec
    window_left = int(peaks - 1*widths)
    window_right = int(peaks + 1*widths)

    if window_right - window_left <= 100:
        window_right = window_right + 300
        window_left = window_left - 300
    if window_left <= 0:
        window_left = 0
    if window_right >= npy.shape[1]:
        window_right = npy.shape[1]
    print('Window (left): ', window_left)
    print('Window (right): ', window_right)
    
    return npy_fils

def pfits_dict():
    pfits = {  
    'gbta': [1.56, 0.39, 0.17, 0.99,
             0.001, 0.001, 0.001, 0.001,
             0.001, 0.001, 0.001, 0.001,
             0.001, 0.001, 0.001, 0.001],
      'gbtb': [4.22, 0.16, 0.18, 0.04,
               4.59, 0.18, 0.001, 2.03,
               5.08, 0.51, 1.14, 0.49,
               5.11, 0.23, 0.1, 0.02],
     #'gbtb(fail)': [6.66-2, 0.55, 0.78, 0.26,
     #         7.06-2, 0.31, 0.14, 0.27,
     #         7.76-2, 0.34, 0.18, 0.45,
     #         8.35-2, 0.26, 0.21, 1.35],
     #'gbtb(fail)': [5.00, 0.69, 0.83, 0.20, 
     #        5.93, 0.23, 0.11, 0.20, 
     #        6.46, 0.26, 0.14, 0.34, 
     #        7.05, 0.20, 0.16, 1.01],
     #'gbtb(fail)': [5.20, 0.48, 1.17, 0.38,
     #        5.06, 0.22, 0.07, 0.14, 
     #        4.10, 0.18, 0.17, 0.25,
     #        4.64, 0.19, 0.10, 1.95],
     'gbtc': [60.80, 0.48, 2.21, 1.61,
              0.001, 0.001, 0.001, 0.001,
              0.001, 0.001, 0.001, 0.001,
              0.001, 0.001, 0.001, 0.001],
    'gbtd': [25.84, 0.66, 0.83, 1.42,
            0.001, 0.001, 0.001, 0.001,
            0.001, 0.001, 0.001, 0.001,
            0.001, 0.001, 0.001, 0.001],
    'gbte': [2.54, 0.76, 0.67, 0.1,
             3.57, 0.70, 0.39, 0.38,
             4.29, 0.46, 0.18, 1.33,
             0.001, 0.001, 0.001, 0.001],
    'gbtf': [8.07, 0.40, 0.41, 0.06,
             9.00, 0.40, 0.29, 3.29,
             11.54, 0.21, 0.12, 1.27,
             13.95, 0.15, 0.04, 0.66],
    'gbtg': [7.46, 0.38, 1.02, 2.46,
             7.58, 0.27, 0.14, 0.03,
             10.001, 0.27, 0.47, 2.43,
             13.68, 0.24, 0.05, 0.68],
    'gmrta': [11.76, 0.78, 2.15, 0.1,
              16.89, 0.80, 1.68, 0.46, 
              20.18, 0.35, 0.83, 4.66,
              0.001, 0.001, 0.001, 0.001],
    'gmrtb': [3.80, 0.96, 0.70, 1.92,
             0.001, 0.001, 0.001, 0.001,
             0.001, 0.001, 0.001, 0.001,
             0.001, 0.001, 0.001, 0.001],
    'gmrtc': [109.41, 0.99, 3.33, 0.09,
             0.001, 0.001, 0.001, 0.001,
             0.001, 0.001, 0.001, 0.001,
             0.001, 0.001, 0.001, 0.001],
    'gmrtd': [28.50, 0.4, 1.16, 1.92,
             29.82, 0.64, 0.08, 0.82,
             0.001, 0.001, 0.001, 0.001,
             0.001, 0.001, 0.001, 0.001]
    }
    pfits_errors = {  
    'gbta': [[0.01, 0.00, 0.01, 0.02]],
    'gbtb': [[0.03, 0.01, 0.01, 0.04],
             [0.01, 0.01, 0.02, 0.12],
             [0.03, 0.01, 0.01, 0.04],
             [0.02, 0.01, 0.01, 0.02]],
     'gbtc': [[0.48, 0.04, 0.39, 0.55]],
     'gbtd': [[0.07, 0.02, 0.07, 0.13]],
     'gbte': [[[0.02,0.01, 0.01, 0.01]], 
              [0.02, 0.02, 0.02, 0.04], 
              [0.01, 0.02, 0.00, 0.04]],
    'gbtf': [[0.01, 00.00, 0.01, 0.01],
             [0.02, 0.0, 0.0, 0.06], 
             [0.01, 0.0, 0.01, 0.07], 
             [0.01, 0.00, 0.00, 0.01]],
    'gbtg': [[0.01, 0.01, 0.0, 0.0],
             [0.02, 0.0, 0.0, 0.04], 
             [0.02, 0.0, 0.01, 0.05], 
             [0.01, 0.01, 0.0, 0.03]],
    'gmrta': [[0.11, 0.02, 0.09, 0.01], 
              [0.1, 0.02, 0.13, 0.08], 
              [0.24, 0.02, 0.07, 0.41]],
    'gmrtb': [[0.03, 0.02, 0.04, 0.09]],
    'gmrtc': [[9.41, 0.08, 0.30, 0.01]],
    'gmrtd': [[1.35, 0.04, 0.09, 0.17], 
              [0.95, 0.06, 0.01, 0.07]]
    }
    
    #Width Errors:
    fwhm_width_errors = {'gbta': 0.01,
    'gbtb': np.sqrt((0.01)**2+(0.01)**2+(0.01)**2+(0.02)**2),
    'gbtc': 0.39,
    'gbtd': 0.07,
    'gbte': np.sqrt((0.01)**2+(0.02)**2+(0.00)**2+(0.01)**2),
    'gbtf': np.sqrt((0.01)**2+(0.01)**2+(0.01)**2+(0.0)**2),
    'gbtg': np.sqrt((0.00)**2+(0.00)**2+(0.01)**2+(0.00)**2),
    'gmrta': np.sqrt((0.09)**2+(0.13)**2+(0.07)**2),
    'gmrtb': 0.04,
    'gmrtc': 0.30,
    'gmrtd': np.sqrt((0.09)**2+(0.01)**2)
    }

    fwhm2_width_errors = {'gbta': 0.01,
    'gbtb': np.sqrt((0.01)**2+(0.02)**2),
    'gbtc': 0.39,
    'gbtd': 0.07,
    'gbte': np.sqrt((0.01)**2+(0.01)**2),
    'gbtf': np.sqrt((0.01)**2+(0.0)**2),
    'gbtg': np.sqrt((0.01)**2+(0.0)**2),
    'gmrta': np.sqrt((0.09)**2+(0.07)**2),
    'gmrtb': 0.04,
    'gmrtc': 0.30,
    'gmrtd': np.sqrt((0.09)**2+(0.01)**2)
    }
    
    return pfits#, pfits_errors, fwhm_width_errors, fwhm2_width_errors

    
def burst_info_all(npy_fils, idx):
    
    pfits, pfits_errors, fwhm_width_errors, fwhm2_width_errors = pfits_dict()
    
    keys = ['gbta', 'gbtb', 'gbtc', 'gbtd', 'gbte', 'gbtf', 'gbtg', 'gmrta', 'gmrtb', 'gmrtc', 'gmrtd']
        
    npy = np.load(npy_fils[idx])
    timeseries = npy.sum(0)
        
    if keys[idx] == 'gbta':
        npy = npy[1600:, :]
    
    pfits_key = keys[idx]
    print('Burst: ', pfits_key)
    bw = [400., 400., 400., 400., 400., 400., 400., 200., 200., 200., 200.] #MHz
    cf = [800., 800., 800., 800., 800., 800., 800., 400., 400., 400., 400.] #MHz
    fd = [83.33, 83.33, 333.33, 166.68, 83.33, 83.33, 83.33, 163.84, 122.88, 400.0, 300.] #ms
    subfactor = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
    print('Archive File Duration (ms): ', fd[idx])
    bandwidth = bw[idx]
    print('Bandwidth (MHz): ', bw[idx])
    center_frequency = cf[idx]
    file_duration = fd[idx]
    sampling_time = (file_duration / npy.shape[1])
    nchan = npy.shape[0]
    freq_res = bandwidth / nchan
    
    #Resolutions
    tres = file_duration / npy.shape[1]
    fres = bandwidth / nchan
    
    peaks, widths, snrs = find_burst(timeseries)
    
    #Define windowing depending on where burst sits in dynspec
    window_left = int(peaks - 1*widths)
    window_right = int(peaks + 1*widths)
    
    if keys[idx] == 'gmrtc' or keys[idx] == 'gmrtd':
        wrr = 50
        wll = 50
    else:
        wrr = 300
        wll = 300
        
    if window_right - window_left <= 100:
        window_right = window_right + wrr
        window_left = window_left - wll
    if window_left <= 0:
        window_left = 0
    if window_right >= npy.shape[1]:
        window_right = npy.shape[1]
    
    y_data = (npy[:].sum(0)/np.max(npy[:].sum(0)))[window_left:window_right]
    x = np.arange(len(y_data)) * sampling_time
    
    y, burst_width, alt_burst_width = exp_gauss_4(x, *pfits[str(pfits_key)])
    
    print(str(keys[idx]) + ' width (ms): ' + str(round(burst_width, 2)) + \
          ' +/- ' + str(round(fwhm_width_errors[str(pfits_key)], 2)))
    
    if keys[idx] == 'gbtb':
        print('No Alternative Burst Width')
        
    else:
        print(str(keys[idx]) + ' width (ms): ' + str(round(alt_burst_width, 2)) + \
              ' +/- ' + str(round(fwhm2_width_errors[str(pfits_key)], 2)))
        
    return x, y, y_data, pfits_key
