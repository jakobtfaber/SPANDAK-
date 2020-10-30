import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import scipy.signal as ss
from matplotlib import rc
from astropy import modeling
from photutils.isophote import EllipseGeometry
from photutils import EllipticalAperture
from photutils.isophote import Ellipse
from photutils.isophote import build_ellipse_model
import sys
plt.rcParams.update({'font.size': 10})
plt.rc('font', family='serif')

def identify_pulse(burst_npy, bandwidth, centfreq, tres=0.08192, plot=True):

	#Load numpy array
	npy = np.load(str(burst_npy))

	#Sub-band data
	subfac = 16
	sub = np.nanmean(npy.reshape(-1, subfac, npy.shape[1]), axis=1)

	#Smooth Data w Savitzky Golay filter (t - time) (s - spectrum)
	twinlen = sub.shape[1] // 4
	if twinlen % 2 == 0:
		twinlen = twinlen + 1
	swinlen = sub.shape[0] // 2
	if swinlen % 2 == 0:
		swinlen = swinlen + 1
	polyo = 6
	savts_2d = ss.savgol_filter(sub, twinlen, polyo, axis=1)
	savsp_2d = ss.savgol_filter(sub, swinlen, polyo, axis=0)
	print('DYN Time Window Length', twinlen)
	print('DYN Freq Window Length', swinlen)

	#Calculate initial guess parameters timeseries - normalized by div by max
	savts = savts_2d.sum(0) / np.max(savts_2d.sum(0))
	maximts = np.max(savts)
	stdts = np.std(savts)
	meants = np.where(savts == maximts)[0][0]
	xts = np.linspace(0, len(savts), len(savts))

	#Fit 1D Gaussian to Timeseries
	fitterts = modeling.fitting.LevMarLSQFitter()
	modelts = modeling.models.Gaussian1D(amplitude=maximts, mean=meants, stddev=stdts)
	fitted_modelts = fitterts(modelts, xts, savts)

	#Crop pulse to 1/2 pulse width on either side
	pulsewidth = 3 * fitted_modelts.parameters[2] #3sigma
	rcut = int(fitted_modelts.parameters[1] + pulsewidth)
	lcut = int(fitted_modelts.parameters[1] - pulsewidth)
	if rcut >= len(savts):
		rcut = len(savts)
	if lcut <= 0:
		lcut = 0
	
	print('Crop Time Bin Right', rcut)
	print('Crop Time Bin Left', lcut)
	sav_c = sub[:, lcut:rcut]

	#Apply proper units to axes
	nchan = npy.shape[1]
	bw = 400 #MHz
	cfreq = 800 #MHz
	fres = bw / nchan
	print('Frequency Resolution', str(fres) + ' MHz')
	subchan = nchan / subfac
	subfres = bw / subchan
	print('Sub Frequency Resolution', str(subfres) + ' MHz')

	#Frequency
	flabellocs = np.arange(0, sav_c.shape[0], step=sav_c.shape[0]//((2*bw/100)))
	flabel_noend = np.arange(cfreq - bw/2, cfreq+ bw/2, step=bw//((2*bw/100)))
	flabels = np.flip(np.append(flabel_noend, [cfreq+ bw/2]))

	#Time
	tbin = sav_c.shape[1]
	tlabellocs = np.arange(0, sav_c.shape[1], step=sav_c.shape[1]//(8+1))
	tlabels_noend = np.arange(0, sav_c.shape[1]*tres, step=(sav_c.shape[1]*tres)//8)
	tlabels = np.append(tlabels_noend, [int(sav_c.shape[1]*tres)])

	#Plot
	if plot == True:
		fig = plt.figure(figsize = (10, 5))
		ax = fig.add_subplot(121)
		plt.imshow(sav_c, aspect = 'auto')
		plt.xticks(tlabellocs, tlabels)
		plt.yticks(flabellocs, flabels)
		plt.ylabel('Frequency (MHz)')
		plt.xlabel('Time (ms)')
		fig.add_subplot(122)
		plt.plot(xts[lcut:rcut], savts[lcut:rcut])
		plt.plot(xts[lcut:rcut], fitted_modelts(xts)[lcut:rcut])
		fig.savefig(str(burst_npy) + '_dyn_ts.png')

	#Crop Smoothed Dyn Spec on Both Axes
	savts_2d_c = savts_2d[:, lcut:rcut]
	savsp_2d_c = savsp_2d[:, lcut:rcut]

	return savts_2d_c, savsp_2d_c, tres, subfres

def prep_2dacf(time_smooth, freq_smooth, time_res, subband_freq_res):

	#Calculate 2d acf
	acf2d = ss.correlate(savts_2d_c, savsp_2d_c)

	#Cap spiked central values in acf
	cap = np.mean(acf2d[len(acf2d.sum(1))//2 +10:len(acf2d.sum(1))//2 +10, len(acf2d.sum(0))//2 +10:len(acf2d.sum(0))//2 +10])
	acf2d_cap = np.where(acf2d > cap, cap, acf2d)

	#Smooth Data w Savitzky Golay filter (t - time) (s - spectrum)
	twinlen = acf2d.shape[1] // 4
	if twinlen % 2 == 0:
	    twinlen = twinlen + 1
	swinlen = acf2d.shape[0] // 4
	if swinlen % 2 == 0:
	    swinlen = swinlen + 1
	polyo = 6
	savacft = ss.savgol_filter(acf2d, twinlen, polyo, axis=1)
	savacff = ss.savgol_filter(acf2d, swinlen, polyo, axis=0)
	print('ACF Time Window Length: ', twinlen)
	print('ACF Freq Window Length: ', swinlen)

	#Calculate initial guess parameters spectrum acf time
	savt = savacft.sum(0) / np.max(savacft.sum(0))
	maximt = np.max(savt)
	stdt = np.std(savt)
	meant = np.where(savt == maximt)[0][0]
	xt = np.linspace(0, len(savt), len(savt))

	#Fit 1D Gaussian to Spectrum
	fittert = modeling.fitting.LevMarLSQFitter()
	modelt = modeling.models.Gaussian1D(amplitude=maximt, mean=meant, stddev=stdt)
	fitted_modelt = fittert(modelt, xt, savt)

	#Calculate initial guess parameters spectrum acf freq
	savsp = savacff.sum(1) / np.max(savacff.sum(1))
	maximsp = np.max(savsp)
	stdsp = np.std(savsp)
	meansp = np.where(savsp == maximsp)[0][0]
	xsp = np.linspace(0, len(savsp), len(savsp))

	#Fit 1D Gaussian to Spectrum
	fittersp = modeling.fitting.LevMarLSQFitter()
	modelsp = modeling.models.Gaussian1D(amplitude=maximsp, mean=meansp, stddev=stdsp)
	fitted_modelsp = fittersp(modelsp, xsp, savsp)

	#Get Ellipse Ratio
	sigmat = fitted_modelt.stddev.value * 0.75
	sigmaf = fitted_modelsp.stddev.value * 0.75

	#Sigmas form a rectangle, get slope of the rectangle diagonal to estimate semi major axis PA
	hyp = np.sqrt(sigmat**2 + sigmaf**2)
	estpa = np.arccos(sigmat / hyp) #in radians

	#Estimate ellipticity (eps) with sigma ratio
	oppestpa = np.arccos(sigmaf / hyp)
	estsmajax = np.tan(oppestpa)*(hyp / 2)
	estsminax = hyp / 2
	esteps = 1 - (estsminax / estsmajax)

	print(estsmajax, estsminax)

	print('Estimated Ellipticity: ', esteps)
	print('Estmated Semimajor Axis: ', estsmajax)
	print('Estimated PA: ', estpa)

	print('Initial guess ellipse applied!')

	#Provide the initial ellipse to be fitted
	#Calculate ellipse geometry
	geometry = EllipseGeometry(x0 = acf2d.shape[1]/2, \
			y0 = acf2d.shape[0]/2, sma = estsmajax, eps = esteps, pa = estpa)
	#Show initial guess ellipse
	aper = EllipticalAperture((geometry.x0, geometry.y0), \
			geometry.sma, geometry.sma*(1-geometry.eps),geometry.pa)

	print('Now for the fit...')
	
	return acf2d, geometry, aper, sigmat, sigmaf

def fit_ellipse(burst_npy, acf, geometry, aper, sigma_t, sigma_f, plot=True):

	#Fit Ellipse to 2D ACF
	try:
		ellipse = Ellipse(acf, geometry)
		isolist = ellipse.fit_image()
		model_image = build_ellipse_model(acf.shape, isolist)
		residual = acf - model_image
		
		if plot == True:
			smas = np.linspace(0, int(sigma_f), 3)
			for sma in smas:
				iso = isolist.get_closest(sma)
				x, y, = iso.sampled_coordinates()
				fig = plt.figure()
				plt.imshow(acf, aspect = 'auto')
				plt.plot(x, y, color='white')
				fig.savefig(str(burst_npy) + '_ellipse_fit.png')
	except OverflowError:
	    print('Note: Overflow Error')
	    pass
	except ValueError:
	    print('Note: Value Error')
	    pass
	except IndexError:
	    print('Ellipse Fit Failed!')
	    pass

	print('Fit completed!')

	slope = np.tan(np.max(isolist.pa))
	std_sma_dr_unround = np.tan(np.std(isolist.pa))
	std_sma_dr = '%s' % float('%.2g' % std_sma_dr_unround)
	drift_rate_unround = -1 * (slope * (subfres / tres)) #MHz/ms
	drift_rate = '%s' % float('%.4g' % drift_rate_unround) #MHz/ms
	print('Drift Slope: ', str(drift_rate) + ' +/- ' + str(std_sma_dr) + ' MHz/ms') 

	return drift_rate

if __name__ == "__main__":

	burst_npy = sys.argv[1]
	savts_2d_c, savsp_2d_c, tres, subfres = identify_pulse(burst_npy, bandwidth=400, centfreq=800, tres=0.04069, plot=True)
	acf2d, geometry, aper, sigmat, sigmaf = prep_2dacf(savts_2d_c, savsp_2d_c, tres, subfres)
	drift_rate = fit_ellipse(burst_npy, acf2d, geometry, aper, sigmat, sigmaf, plot=True)



