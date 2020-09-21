import matplotlib
matplotlib.use('Agg')
import numpy as np
import psrchive
import pylab
import sys
import os

#fits_dir = os.fsencode('/datax/scratch/jfaber/SPANDAK_extension/pipeline_playgroune/61.4627973333_67.0552026667_fits')

directory = r'/datax/scratch/jfaber/SPANDAK_extension/pipeline_playground/R3/0002/16.384175488_20.223649024_3.8_9/fits'#/FRB121102_fits/1703.15379733_1708.74620267_fits'
def fits2numpy():
	for fits in os.listdir(directory):
		#print(fits)
		#npar = 'pulse_120390656' + '_secondtry.npy'
		if fits.endswith('.fits'):
			npar = str(fits) + '.npy'
			with open(npar, 'wb') as npar_file:		
				#arch = psrchive.Archive_load('/datax/scratch/jfaber/SPANDAK_extension/pipeline_playground/61.4627973333_67.0552026667_fits/pulse_120390656.fits')
				arch = psrchive.Archive_load(directory + '/' + fits)
				#os.system('psrplot -p F -jD' + directory + '/' + fits)
				arch.dedisperse()
				arch.remove_baseline()
				arch.convert_state('Stokes')
				data = arch.get_data()
				np.save(npar_file, data[:, 0, :, :].mean(0))
				print('Array Written...')

fits2numpy()

