import sys
import os

os.system('source /home/vgajjar/spandakenv/bin/activate')

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import scipy.signal as ss
from scipy import stats
import numpy as np
import psrchive
import pylab

def fits2numpy(fitsdir):
	npydir = str(fitsdir) + '/npys'
	if not os.path.exists(npydir):
		os.mkdir(npydir)
	os.chdir(npydir)
	for fits in os.listdir(fitsdir):
		#print(fits)
		#npar = 'pulse_120390656' + '_secondtry.npy'
		if fits.endswith('.fits'):
			npar = str(fits) + '.npy'
			with open(npar, 'wb') as npar_file:		
				#arch = psrchive.Archive_load('/datax/scratch/jfaber/SPANDAK_extension/pipeline_playground/61.4627973333_67.0552026667_fits/pulse_120390656.fits')
				#arch = psrchive.Archive_load(directory + '/' + fits)
				arch = psrchive.Archive_load(fitsdir + '/' + fits)
				#os.system('psrplot -p F -jD' + directory + '/' + fits)
				arch.dedisperse()
				arch.remove_baseline()
				arch.convert_state('Stokes')
				data = arch.get_data()

				#Write Numpy Arrays to npys directory in fits directory
				np.save(npar_file, data[:, 0, :, :].mean(0))
				print('Numpy Array Written...')
	return npydir

def id_cand(npydir):

	npy_fils = [i for i in os.listdir(npydir) if i.endswith('.npy')]#[1:]

	png_directory = str(npydir) + '/pngs'
	if not os.path.exists(png_directory):	
		os.mkdir(png_directory)
	os.chdir(png_directory)	

	pulse_files = []
	npy_snrs = []

	for fil in npy_fils:    

		pulse_files.append(fil)

		#Load npy file    
		ar = np.load(npydir + '/' + fil)

		#Sub-band npy array
		sub_factor = 32
		ar_sb = np.nanmean(ar.reshape(-1, sub_factor, ar.shape[1]), axis=1)

		#Integrate to get absolute-valued and normalized timeseries to be able to calculate 'snr'
		ar_ts = np.abs(ar_sb.sum(0) / np.max(ar_sb.sum(0)))

		fig = plt.figure(figsize = (20, 10))
		ax1 = fig.add_subplot(121)

		#Smooth timeseries with Savitzky Golay filter
		ts_sg = ss.savgol_filter(ar_ts, 115, 9)[100:-100]
		ts_sg_snr = 10 * np.log10(np.max(ts_sg) / np.mean(ts_sg))
		#print('SNR: ', ts_sg_snr)
		#print('Pulse File: ', fil)

		#Signal search for peaks, and normalized peak prominence
		ar_pks = ss.find_peaks(ts_sg)
		ar_pk_prom = ss.peak_prominences(ts_sg, ar_pks[0])[0]
		norm_pk_prom = ar_pk_prom / np.max(ar_pk_prom)
		peak_prom_snr = 10 * np.log10(np.max(norm_pk_prom) / np.min(norm_pk_prom))
		npy_snrs.append(ts_sg_snr)
		#print('Prominence SNR: ', peak_prom_snr)
		
		plt.title('Time Series | Peak Prominence SNR: ' + str(peak_prom_snr))
		plt.plot(ts_sg)
		#print('Peak Prominences: ', norm_pk_prom)


		#Plot diagnostic dynamic spectrum
		if ts_sg_snr > 5: 
			ax2 = fig.add_subplot(122)
			plt.title('Dynamic Spectrum | Candidate Likely | SNR: ' + str(ts_sg_snr))
			plt.imshow(ar_sb, aspect = 'auto')
			plt.gca().invert_yaxis()
			plt.savefig(fil + '_' +  str(ts_sg_snr) + '.png')
			#plt.show()
		else:
			ax2 = fig.add_subplot(122)
			plt.title('Dynamic Spectrum | Candidate Unlikely | SNR: ' + str(ts_sg_snr))
			plt.imshow(ar_sb, aspect = 'auto')
			plt.gca().invert_yaxis()
			plt.savefig(fil + '_' +  str(ts_sg_snr) + '.png')
			#plt.savefig(fil + '.png')
	
	pulse_fits = pulse_files[np.where(npy_snrs == np.max(npy_snrs[1:-1]))[0][0]][:-4]
	
	return pulse_fits


def polfluxcal(pulse_fits, calib_file_path):
	
	new_calibdir = str(fitsdir) + '/calib'
	if not os.path.exists(new_calibdir):	
		os.mkdir(new_calibdir)
	os.chdir(new_calibdir)
	os.system('cp ' + str(calib_file_path) + '/* .')
	os.system('cp ' + str(fitsdir) + '/' + str(pulse_fits) +  ' .')

	database = 'pac -wp . -u fits'
	os.system(database)
	print('Database Command', database)
	print('Database created')
	fluxcal = 'fluxcal -i 15 -d database.txt -c fluxcal.cfg'
	os.system(fluxcal)
	print('Fluxcal Command', fluxcal)
	print('Flux & Pol calibration initiated')
	cfreq_adjust = 'psredit -c "freq=6407.714800" -m ' + str(pulse_fits)
	os.system(cfreq_adjust)
	print('Center Frequency Adjusted to Exactly 6407.414800 MHz')
	calib = 'pac -x -d database.txt ' + str(pulse_fits)
	os.system(calib)
	print('Calibration Command', calib)
	print('Calibration complete')
	return new_calibdir

def rmfit(pulse_fits, new_calibdir):

	os.chdir(new_calibdir)

	RM_fit_command = 'python ' + '/datax/scratch/jfaber/SPANDAK_extension/pipeline_playground/RMfit_curve.py ' + str(pulse_fits)[:-4] + '.calib'
	os.system(RM_fit_command)
	print('RM_fit Command', RM_fit_command)
	return



if __name__ == '__main__':

	fitsdir = sys.argv[1]
	calib_file_path = sys.argv[2]
	npydir = str(fitsdir) + '/npys'
	npydir = fits2numpy(fitsdir)
	print('Numpy Arrays Written')
	pulse_fits = id_cand(npydir)
	print('Pulse Fits File Identified (Remember to Check Bookend PNGs)', pulse_fits)
	#pulse_fits = 'pulse_120303971.fits'	
	new_calibdir = polfluxcal(pulse_fits, calib_file_path)
	print('Polarization and Flux Calibration Complete')
	#rmfit(pulse_fits, new_calibdir)
	#print('Rotation Measure Fitting Complete')









