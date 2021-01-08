import matplotlib
matplotlib.use('Agg')
import numpy as np
import psrchive
import pylab
import sys
import os

#Directory containing fits files
directory = sys.argv[1]
end = sys.argv[2]
DM = sys.argv[3]
burstname = sys.argv[4]


def fits2numpy():
        for fits in os.listdir(directory):
		if fits.endswith(str(end)): #or fits.endswith('.tsch256') or fits.endswith('.fits'):
                        #npar = str(fits.split('.')[0]) + '.dm' + str(DM) + '.npy'
                        npar = str(burstname) + '.dm' + str(DM) + '.npy'
			with open(npar, 'wb') as npar_file:
                                arch = psrchive.Archive_load(directory + '/' + fits)
                                arch.pscrunch()
				arch.bscrunch_to_nbin(64)
				arch.fscrunch_to_nchan(608)
				arch.set_dispersion_measure(float(DM))
				arch.dedisperse()
                                arch.remove_baseline()
                                #arch.convert_state('Stokes')
                                data = arch.get_data()
                                np.save(npar_file, data[:, 0, :, :].mean(0))
                                print('Array Written...')

fits2numpy()
