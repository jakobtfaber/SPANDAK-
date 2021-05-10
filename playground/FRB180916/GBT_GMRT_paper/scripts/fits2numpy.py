import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import psrchive
import pylab
import sys
import os

#fits_dir = os.fsencode('/datax/scratch/jfaber/SPANDAK_extension/pipeline_playgroune/61.4627973333_67.0552026667_fits')

#directory = r'/datax/scratch/jfaber/SPANDAK_extension/pipeline_playground/R3_GBT_1/R3_analysis/D_1164/fits'

directory = sys.argv[1]

def fits2numpy(directory):
	for fits in os.listdir(directory):
		#print(fits)
		#npar = 'pulse_120390656' + '_secondtry.npy'
		if fits.endswith('.fits') or fits.endswith('.norm'):
			npar = str(fits) + '.npy'
			with open(npar, 'wb') as npar_file:		
				#arch = psrchive.Archive_load('/datax/scratch/jfaber/SPANDAK_extension/pipeline_playground/61.4627973333_67.0552026667_fits/pulse_120390656.fits')
				arch = psrchive.Archive_load(directory + '/' + fits)
				#os.system('psrplot -p F -jD' + directory + '/' + fits)
				arch.setnbin(256)
				arch.dedisperse()
				arch.remove_baseline()
				arch.convert_state('Stokes')
				data = arch.get_data()
				np.save(npar_file, data[:, 0, :, :].mean(0))
				print('Array Written...')


def plot(directory):
	for npy in os.listdir(directory):
		if npy.endswith('.npy'):
			ar = np.load(str(npy))
			plt.imshow(ar)
			plt.show()


if __name__ == 'main':

	fits2numpy(directory)
	#plot(directory)

