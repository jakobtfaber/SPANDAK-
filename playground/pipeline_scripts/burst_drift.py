import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import scipy.signal as ss
from matplotlib import rc
from astropy.modeling.models import Gaussian2D
from photutils.isophote import EllipseGeometry
from photutils import EllipticalAperture
from photutils.isophote import Ellipse
from photutils.isophote import build_ellipse_model
plt.rcParams.update({'font.size': 20})
plt.rc('font', family='serif')

def burst_drift(npy, fslice = ':', tslice = ':'):
	
	#Load npy array
	burst = np.load(npy)

	#Subband in Frequency
	subfac = 16
	rb_sub = np.nanmean(R3H.reshape(-1, subfac, R3H.shape[1]), axis=1)

	rb_sub_zm = rb_sub[fslice,tslice]

	ynorm = (rb_sub_zm.sum(0)/np.max(rb_sub_zm.sum(0)))
	x = np.linspace(0, len(ynorm)-1, len(ynorm))

	#Smooth 

	sav = ss.savgol_filter(rb_sub_zm, 49, 6)

	#Calculate 2D ACF

	acf = ss.correlate(sav, sav)

	#Provide the initial Ellipse to be fitted
	#Calculate Ellipse Geometry
	geometry = EllipseGeometry(x0 = acf.shape[1]/2, y0 = acf.shape[0]/2, sma = 20, eps = 0.4, pa = 60 * np.pi/180.)
	#Show Initial Guess Ellipsee
	aper = EllipticalAperture((geometry.x0, geometry.y0), geometry.sma, geometry.sma*(1-geometry.eps),geometry.pa)

	#Plot Initial Guess Ellipse on ACF

	#plt.imshow(acf, aspect = 'auto')
	#aper.plot(color='white')

	#Fit Ellipse to 2D ACF

	ellipse = Ellipse(acf, geometry)
	isolist = ellipse.fit_image()
	model_image = build_ellipse_model(acf.shape, isolist)
	residual = acf - model_image

	fig = plt.figure(figsize = (40, 10))
	ax1 = fig.add_subplot(131)
	plt.gca().invert_yaxis()
	x = np.linspace(0, acf.shape[1] -1, acf.shape[1])
	y = lambda x: (0.809) * x + 37.65


	plt.plot(x, y(x), color = 'white', label = 'Drift Rate = -33 MHz/ms')
	plt.imshow(sav,interpolation = None)
	plt.yticks(np.arange(0, sav.shape[0], 60), [1000, 900, 800, 700, 600])
	plt.xticks(np.arange(0, sav.shape[1], 48), [0,2,4,6,8,10])
	plt.ylabel('Frequency (MHz)')
	plt.xlabel('Time (ms)')
	plt.title('Dynamic Spectrum R3H')
	plt.legend()

	ax2 = fig.add_subplot(132)

	plt.imshow(acf)#, origin='lower')
	plt.ylabel('Frequency Shift (MHz)')
	plt.xlabel('Time Shfit (ms)')
	plt.yticks(np.arange(0, acf.shape[0], 79), [1000, 600, 200, 0, -200, -600, -1000])
	plt.xticks(np.arange(0, acf.shape[1], 49), [-10, -8, -6, -4, -3, 0, 2, 4, 6, 8, 10])

	#plt.plot(x, y(x))
	plt.title('2D ACF R3H')

	smas = np.linspace(10, 100, 8)
	for sma in smas:
		iso = isolist.get_closest(sma)
		x, y, = iso.sampled_coordinates()
		plt.plot(x, y, color='white')

	#ax3.imshow(model_image, origin='lower')
	#ax3.set_title('Ellipse Model')

	#ax3.imshow(residual, origin='lower')

	ax3 = fig.add_subplot(133)

	#plt.set_aspect('auto')
	x = np.linspace(0, len(ynorm)-1, len(ynorm))
	plt.plot(x, ss.savgol_filter(ynorm, 19, 6))
	plt.title('Timeseries R3H')
	plt.ylabel('Normalized Intensity')
	plt.xlabel('Time (ms)')
	plt.xticks(np.arange(0, sav.shape[1], 48), [0,2,4,6,8,10])
	#ax3.set_title('Residual')
	fig.savefig('Drift Rate R3H')
	return

if __name__ == '__main__':

	npy = sys.argv[1]
	burst_drift(npy)




