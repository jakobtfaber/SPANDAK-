import numpy as np
import matplotlib.pyplot as plt
import fnmatch
import sys
import os
import pandas as pd
import csv
import itertools

def _cdd_auto(cdd=True):

	'''Performs coherent dedispersion with DSPSR on spliced raw voltages'''

	#Start times
	start_times = [117, 686, 1164, 1549, 639, 267, 579]
	par_fil_paths = ['A_117/A_117.par', 'B_686/B_686.par', 'C_1164/C_1164.par', 'D_1549/D_1549.par', 'E_639/E_639.par', 'F_267/F_267.par', 'G_579/G_579.par']

	#Specify CDD parameters
	cepoch = 58178
	polar = 4
	phasebin = 2048
	p = 0
	chan = 1024
	samples = 1024 #number of MB

	#Parse and Form Coherent Dispersion Commands
	cdd_run_commands = []

	for B in np.arange(len(start_times)):
		cdd_run = 'dspsr ' + '-U ' + str(samples) + ' -F ' + str(chan) + ':D ' \
		+ ' -K ' + ' -d ' + str(polar) + ' -b  ' + str(phasebin) + ' -E ' \
		+ '/datax/scratch/jfaber/SPANDAK++/pipeline_playground/R3_GBT_1/GBT_GMRT_paper/bursts_bpc/' + par_fil_paths[B] + ' -s -a psrfits -e fits ' \
		+ '/datax/scratch/jfaber/SPANDAK++/pipeline_playground/R3_GBT_1/GBT_GMRT_paper/bursts_bpc/' + str(par_fil_paths[B].split('/')[0])  \
		+ '/spliced_' \
		+ str(start_times[B]) + '.raw'
		cdd_run_commands.append(cdd_run)

	#Coherently Dedisperse Spliced Raw Voltages
	if cdd == True:

		for cdd in cdd_run_commands:
			
			#Make directies for raw voltages (if necessary -- probably won't need it), but also for fits file storage
			time_dir = '/datax/scratch/jfaber/SPANDAK++/pipeline_playground/R3_GBT_1/GBT_GMRT_paper/bursts_bpc/' + str(cdd.split('/')[9])
			fits_dir = '/datax/scratch/jfaber/SPANDAK++/pipeline_playground/R3_GBT_1/GBT_GMRT_paper/bursts_bpc/' + str(cdd.split('/')[9]) + '/fits'
			#print('Time', time_dir)
			#print('Fits', fits_dir)

			if not os.path.exists(time_dir):
				os.mkdir(time_dir)
			if not os.path.exists(fits_dir):
				os.mkdir(fits_dir) 
#
			##Funnel dspsr output into fits file directory made previously
			os.chdir(fits_dir)
			print('DSPSR Output Funnelling Into: ' + os.getcwd())
			print(cdd)
			os.system(cdd)
			print('Coherent Dedispersion Complete')



	return cdd_run_commands

if __name__ == '__main__':

	cdd_run_commands = _cdd_auto(cdd=True)
