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


def polfluxcal(fitsdir, pulse_fits, calib_file_path):

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

        RM_fit_command = 'python ' + '/datax/scratch/jfaber/FLITS/playground/FRB121102/rm_pa_fit/RMfit_curve.py ' + str(pulse_fits)[:-3] + '.calib'
        os.system(RM_fit_command)
        print('RM_fit Command', RM_fit_command)
        return



if __name__ == '__main__':

        fitsdir = sys.argv[1]
        pulse_fits = sys.argv[2]
	calib_file_path = sys.argv[3]
        print('Pulse Fits File Identified (Remember to Check Bookend PNGs)', pulse_fits)
        new_calibdir = polfluxcal(fitsdir, pulse_fits, calib_file_path)
        print('Polarization and Flux Calibration Complete')
        rmfit(pulse_fits, new_calibdir)
        print('Rotation Measure Fitting Complete')
