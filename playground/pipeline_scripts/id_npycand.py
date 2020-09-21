import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import scipy.signal as ss
from scipy import stats
import sys
import os
#sys.path.insert(0, os.path.abspath('../fits2npy_test_files'))

#Path to fits directory
#f2ndir = '/Users/jakobfaber/Documents/spandak_extended/SPANDAK_extension/pipeline_playground/fits2npy_test_files'
f2ndir = '/datax/scratch/jfaber/SPANDAK_extension/pipeline_playground/R3/0002/16.384175488_20.223649024_3.8_9/npys'#FRB121102_npy/1703.15379733_1708.74620267_npy'

def id_cand(f2ndir=f2ndir):
    
    npy_fils = [i for i in os.listdir(f2ndir) if i.endswith('fits.npy')]#[1:]

    for fil in npy_fils:    

        #Load npy file    
        ar = np.load(f2ndir + '/' + fil)

        #Sub-band npy array
        sub_factor = 32
        ar_sb = np.nanmean(ar.reshape(-1, sub_factor, ar.shape[1]), axis=1)

        #Integrate to get absolute-valued and normalized timeseries and calculate 'snr'
        ar_ts = np.abs(ar_sb.sum(0) / np.max(ar_sb.sum(0)))
        fig = plt.figure(figsize = (20, 10))
        ax1 = fig.add_subplot(121)
        #plt.plot(ar_ts)

        #Smooth timeseries with Savitzky Golay filter
        ts_sg = ss.savgol_filter(ar_ts, 115, 9)[100:1900]
        ts_sg_snr = 10 * np.log10(np.max(ts_sg) / np.min(ts_sg))
        print('SNR: ', ts_sg_snr)
	print('Pulse File: ', fil)

        #Signal search for peaks, and normalized peak prominence
        ar_pks = ss.find_peaks(ts_sg)
        #print('Peaks: ', ar_pks[0])
        ar_pk_prom = ss.peak_prominences(ts_sg, ar_pks[0])[0]
        norm_pk_prom = ar_pk_prom / np.max(ar_pk_prom)
        peak_prom_snr = 10 * np.log10(np.max(norm_pk_prom) / np.min(norm_pk_prom))
        print('Prominence SNR: ', peak_prom_snr)
        plt.title('Time Series | Peak Prominence SNR: ' + str(peak_prom_snr))
	plt.plot(ts_sg)
	#print('Peak Prominences: ', norm_pk_prom)


        #Plot dynamic spectrum
        if ts_sg_snr: 
	#ar_corr = ss.correlate(ar_sb, ar_sb, mode = 'full')
        	ax2 = fig.add_subplot(122)
        	plt.title('Dynamic Spectrum | Candidate Detected! | SNR: ' + str(ts_sg_snr))
		plt.imshow(ar_sb, aspect = 'auto')
        	plt.gca().invert_yaxis()
		plt.savefig(fil + '_' +  str(ts_sg_snr) + '.png')
		#plt.show()
	else:
		ax2 = fig.add_subplot(122)
                plt.title('Dynamic Spectrum | No Candidate Found | SNR: ' + str(ts_sg_snr))
                plt.imshow(ar_sb, aspect = 'auto')
                #plt.gca().invert_yaxis()
                #plt.savefig(fil + '.png')

id_cand()
