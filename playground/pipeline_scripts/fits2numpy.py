import matplotlib
matplotlib.use('Agg')
import numpy as np
import psrchive
import pylab
import sys
import os

#Directory containing fits files
directory = argsv[1]

def fits2numpy():
        for fits in os.listdir(directory):
                if fits.endswith('.fits') or fits.endswith('.ar') or fits.endswith('.tsch'):
                        npar = str(fits.split('.')[0]) + '.npy'
                        with open(npar, 'wb') as npar_file:
                                arch = psrchive.Archive_load(directory + '/' + fits)
                                arch.dedisperse()
                                arch.remove_baseline()
                                #arch.convert_state('Stokes')
                                data = arch.get_data()
                                np.save(npar_file, data[:, 0, :, :].mean(0))
                                print('Array Written...')

fits2numpy()
