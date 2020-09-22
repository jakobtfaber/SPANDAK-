#!/usr/bin/env python

from argparse import ArgumentParser
import os,sys,math
import timeit
import optparse
import numpy as np
import glob
from itertools import chain
import smtplib
from os.path import basename
import subprocess as sb
import numpy as np
import pypulse as pp
import matplotlib.pyplot as plt
from scipy.signal import detrend 
from peakutils.baseline import baseline
import pandas as pd
import psrchive as psr
import scipy.ndimage as nd
from scipy.optimize import curve_fit
from scipy import asarray as ar,exp

def run(gonogo,cmd):
	cmd = cmd.split(" ")
	if(gonogo): 
		#os.system(cmd)
		print cmd
		p = sb.Popen(cmd)
	else:	p = sb.Popen('echo')
	return p 

def simrun(execute,cmd):
	print cmd
	if(execute):
                os.system(cmd)

  	
def thresholding_algo(y, lag, threshold, influence):
    signals = np.zeros(len(y))
    filteredY = np.array(y)
    avgFilter = [0]*len(y)
    stdFilter = [0]*len(y)
    avgFilter[lag - 1] = np.mean(y[0:lag])
    stdFilter[lag - 1] = np.std(y[0:lag])
    for i in range(lag, len(y)):
        if abs(y[i] - avgFilter[i-1]) > threshold * stdFilter [i-1]:
            if y[i] > avgFilter[i-1]:
                signals[i] = 1
            else:
                signals[i] = -1

            filteredY[i] = influence * y[i] + (1 - influence) * filteredY[i-1]
            avgFilter[i] = np.mean(filteredY[(i-lag):i])
            stdFilter[i] = np.std(filteredY[(i-lag):i])
        else:
            signals[i] = 0
            filteredY[i] = y[i]
            avgFilter[i] = np.mean(filteredY[(i-lag):i])
            stdFilter[i] = np.std(filteredY[(i-lag):i])

    return dict(signals = np.asarray(signals),
                avgFilter = np.asarray(avgFilter),
                stdFilter = np.asarray(stdFilter))

def find_nearidx(array,val):
	idx = (np.abs(array-val)).argmin()	
	return idx



def structure_para(filename,freql,freqh,snr,dm,onp):
    '''
    Copied from original to just measure rms across on-pulse
    '''

    #Using PyPulse
    '''
    data = pp.Archive(filename,prepare=False)
    data1 = data.data[0][0]
    #time = np.sum(data1,axis=0)	
    freq = data.freq[0]
    '''

    #Using PSRCHIVE
    #'''
    fpsr = psr.Archive_load(filename)
    fpsr.remove_baseline()	
    fpsr.set_dispersion_measure(0)
    fpsr.dedisperse()
    fpsr.set_dispersion_measure(dm)
    fpsr.dedisperse()
    ds = fpsr.get_data().squeeze()
    w = fpsr.get_weights().flatten()
    w = w/np.max(w) # Normalized it
    idx = np.where(w==0)[0]
    ds = np.multiply(ds, w[np.newaxis,:,np.newaxis]) # Apply it
    ds[:,idx,:] = np.nan
    data1 = ds[0,:,:]	
    #print np.shape(data1)
    freq=fpsr.get_frequencies()
    #''' 
    #print freq,freql,freqh
    #print find_nearidx(freq,freql),find_nearidx(freq,freqh)
    #plt.xlim(283,332)
    #plt.imshow(data1,aspect='auto')
    #plt.show()

    #time = np.sum(data1[find_nearidx(freq,freqh):find_nearidx(freq,freql),:],axis=0)
    #time = np.nansum(data1[find_nearidx(freq,freql):find_nearidx(freq,freqh),:],axis=0)

    time = np.nansum(data1,axis=0)	
    time = time - baseline(time,6)
    #time = (time - np.mean(time))/np.max(time)
 
    #tbin = list(data1.durations/data1[0].size)[0]
    #taxis = np.arange(0,data1.durations,tbin)

    # Settings 
    lag = 200
    #lag = 30
    threshold = snr
    influence = 0.4
    ext = 20 # Extra phase bins around the pulse
    #-----------
    
    y = time
    # Run algo with settings from above
    '''
    result = thresholding_algo(time, lag=lag, threshold=threshold, influence=influence)
  
    res = np.where(result['signals']==1)
    rescen = int(np.nanmean(res))

    low1 = rescen - ext
    high1 = rescen + ext
    
     
    #Orig	
    low = np.min(res)-1
    high = np.max(res)+1
    '''
    # Test
    #low = 301
    #high = 317

    #For TschFsch file
    #low1 = 282
    #high1 = 332  


    #For TscFsch2 file
    #low1 = 1
    #high1 = 128

    #For TschFsch3 file
    #low1 = 1
    #high1 = 256	
    
    low1 = onp[0]
    high1 = onp[1]
	
    #print time[low:high]
    #spec = np.sum(data1[:,low:high],axis=1)
    onpulse = time[low1:high1]  
    #onpulse = onpulse - np.mean(onpulse) 

    offpulse = time[low1+4*ext:high1+4*ext] 
    #print low1+4*ext,high1+4*ext

    #error = np.std(offpulse)	
    #strct_para = np.mean(abs(np.diff(abs(np.diff(onpulse)))))

    #Orig	
    #strct_para = np.mean(abs(np.diff(onpulse)))


    #DMgrad = abs(np.diff(time))
    DMgrad = abs(np.gradient(time,2))
    #DMgrad = time
    strct_para = np.mean(DMgrad)
  
    error = np.mean(abs(np.diff(offpulse)))
    #print dm,strct_para,error,low1,high1

    #plt.xlim(low1,high1)
    #plt.plot(time[low1:high1])
    #plt.plot(time)
    #plt.plot(onpulse)
    #plt.plot(np.diff(onpulse))	
    #plt.plot(offpulse)	
    #plt.show()	

    #oname = "".join(filename.split(".")[:-1]) + ".eps"
    #specname = "".join(filename.split(".")[:-1]) + "_spectra.txt"	
    '''	
    oname = filename + ".eps"
    specname = filename + "_spectra.txt"
	 
    np.savetxt(specname,spec,fmt="%.2f")	
    #np.savetxt(specname,np.c_[freq,spec],fmt="%.2f %.2f")	

    plt.rcParams['axes.linewidth'] = 2
    plt.subplots_adjust(hspace = .001)
    plt.subplots_adjust(wspace = .001)
    ax1 = plt.subplot2grid((6,5), (0,0), rowspan=2,colspan=4)

    #plt.subplot(211)
    plt.xlim(low-ext,high+ext)
    plt.setp(ax1.get_xticklabels(), visible=False)
    #plt.setp(ax1.get_yticklabels(), visible=False
    ax1.set_ylabel('Flux (mJy)',fontsize=24, fontweight='bold')
    plt.tick_params(axis='both', which='major', labelsize=22)
    ax1.yaxis.set_ticks(np.arange(0,max(y),max(y)/5))
    #plt.yticks(rotation=90)
    #plt.locator_params(axis='y', nticks=4)
    ax1.plot(np.arange(1, len(y)+1), y,linewidth=2,color='black')

    #plt.subplot(212)
    #plt.step(np.arange(1, len(y)+1), result["signals"], color="red", lw=2)
    ax2 = plt.subplot2grid((6,5), (2,0), rowspan=4,colspan=4)
    #plt.xlim(low-20,high+20)
    pdata = data1[:,low-ext:high+ext]
    ptime = taxis[low-ext:high+ext]
    lowedge = (taxis[low]-np.mean(ptime))*1000 # msec
    highedge = (taxis[high]-np.mean(ptime))*1000 # msec	
    ptime=(ptime-np.mean(ptime))*1000 # msec 
    plt.yticks(rotation=90)
    ax2.tick_params(length=4, width=2)
    plt.tick_params(axis='both', which='major', labelsize=22)
    ax2.set_ylabel('Frequency (MHz)',fontsize=24, fontweight='bold')
    ax2.set_xlabel('Time (msec)',fontsize=24, fontweight='bold')
    plt.axvline(lowedge, color='b', linestyle='dashed', linewidth=2)
    plt.axvline(highedge, color='b', linestyle='dashed', linewidth=2)
    plt.axhline(freqh,color='b', linestyle='dashed', linewidth=2)
    plt.axhline(freql, color='b', linestyle='dashed', linewidth=2)
    plt.imshow(pdata,aspect='auto',cmap='binary',extent=[min(ptime),max(ptime),min(freq),max(freq)],interpolation='none')

    ax3 = plt.subplot2grid((6,5), (2,4), rowspan=4,colspan=1)
    plt.setp(ax3.get_xticklabels(), visible=False)
    plt.setp(ax3.get_yticklabels(), visible=False)
    plt.ylim(min(freq),max(freq))
    plt.plot(spec,freq,linewidth=2,color='black')
  
    fig = plt.gcf()
    fig.set_size_inches(8,20) 
    plt.savefig(oname, bbox_inches='tight',dpi=300)	
    '''
    return strct_para,DMgrad,time

def gaus(x,a,x0,sigma,c):
    return a*exp(-(x-x0)**2/(2*sigma**2)) + c


if __name__ == "__main__":

	parser = optparse.OptionParser()

	parser.add_option("-f", action='store', dest='infile', type=str, help="Input archive file")
	
	parser.add_option("-o", action='store', dest='outdir', default="",type=str,help="Full Output directory (Default : .)")
	
	parser.add_option("--plot", action='store_true', dest='plot',help='Plot the final profile in eps file (Default: DO not plot)')

	parser.add_option("--fl", action='store', dest='freql', default=4000, type=float,help="Lower frequency for the pulse to add (this only affect the average pulse)")

	parser.add_option("--fh", action='store', dest='freqh', default=8000, type=float,help="Higher frequency for the pulse to add (this only affect the average pulse)")

	parser.add_option("-s", action='store', dest='snr', default=6, type=float,help="SNR threshold to detect pulse")

	options,args = parser.parse_args()

	infile = os.path.abspath(options.infile)
	plot = options.plot
	freql = options.freql
	freqh = options.freqh
	snr = options.snr

	if not options.outdir: outdir = os.getcwd()
        else:
                outdir = options.outdir
                if(os.path.isdir(outdir) is not True):
                        os.system("mkdir %s" % (outdir))	

        DM = 350
	DMarr = []
	DMgradArr = []
	DMgradArr_off = []
	TimeArr = []
	TimeArr_off = []

 	strct_para_array = []

	#Orig
	DMl = 330
	DMh = 370
	#best fit polynomial fit range
	#DMl = 340
	#DMh = 355
	#TcshFcsh2
	onp = [72,83]
	#TcshFcsh3
	#onp = [140,170]
	offp = [10,21]

	#plotfig2(infile,freql,freqh,snr)		
	#for DM in range(520,620,1):
	for DM in np.linspace(DMl,DMh,1000):
	#for DM in range(340,360,1):
		s,d,t = structure_para(infile,freql,freqh,snr,DM,onp)		
		strct_para_array.append(s)
		s1,d1,t1 = structure_para(infile,freql,freqh,snr,DM,offp)
		TimeArr.append(t)
		TimeArr_off.append(t)
		DMarr.append(DM)
		DMgradArr.append(d)
		DMgradArr_off.append(d1)
	
	fpsr2 = psr.Archive_load(infile)
	tbin = float(fpsr2.integration_length()/fpsr2.get_nbin())
	taxis = np.arange(0,fpsr2.integration_length(),tbin)*1000 - 300 # To get to msec

	DMgradArr = np.array(DMgradArr)
	DMgradArr_off = np.array(DMgradArr)
	TimeArr = np.array(TimeArr)
	TimeArr_off = np.array(TimeArr)
	
	#DMgradArr = nd.filters.gaussian_filter(DMgradArr,sigma=2,mode='constant')
	DMgradArr = nd.filters.uniform_filter(DMgradArr,size=(3,3),mode='constant')
	DMgradArr_off = nd.filters.uniform_filter(DMgradArr_off,size=(3,3),mode='constant')
	
	
	#TimeArr = nd.filters.gaussian_filter(TimeArr,sigma=1,mode='constant')
	TimeArr = nd.filters.uniform_filter(TimeArr,size=(3,3),mode='constant')
	TimeArr_off = nd.filters.uniform_filter(TimeArr_off,size=(3,3),mode='constant')

	DMgradArr/=DMgradArr.max()
	TimeArr/=TimeArr.max()

	DMgradArr_off/=DMgradArr_off.max()
	TimeArr_off/=TimeArr_off.max()


	#plt.plot(DMarr,strct_para_array)
	#plt.imshow(DMgradArr.T,aspect='auto',extent=(0,len(DMgradArr[0]),DMl,DMh),origin='lower')
	#plt.imshow(DMgradArr,aspect='auto',origin='lower',extent=(0,len(DMgradArr[0]),DMl,DMh))
	
	
	'''
	plt.xlim(40,110)
	#plt.xlim(onp[0],onp[1])

	plt.imshow(TimeArr,aspect='auto',origin='lower',extent=(0,len(DMgradArr[0]),DMl,DMh),cmap='gray_r')	

	Y = np.linspace(DMl,DMh,1000)
	X = range(0,len(DMgradArr[0]))
	plt.contour(X,Y,DMgradArr,10,aspect='auto',origin='lower',extent=(0,len(DMgradArr[0]),DMl,DMh),cmap='rainbow')

	plt.xlabel('Time bins')
	plt.ylabel('DM')
	#plt.savefig('DM_SNR_struct.pdf')
	plt.figure(2)
	
	#TschFsch2 
	#avgDMgrad = np.mean(DMgradArr[:,50:100],axis=1)
	#avgTime = np.mean(TimeArr[:,50:100],axis=1)

	#TschFsch3 
	avgDMgrad = np.mean(DMgradArr[:,onp[0]:onp[1]],axis=1)
	avgTime = np.mean(TimeArr[:,onp[0]:onp[1]],axis=1)
	'''

	plt.rcParams['axes.linewidth'] = 2
	plt.subplots_adjust(hspace = 0.001)
	plt.subplots_adjust(wspace = 0.001) 

	plt.tick_params(axis='both', which='major', labelsize=14)

	ax1 = plt.subplot2grid((1,4),(0,0), rowspan=1,colspan=3)
	#ax1.set_xlim(onp[0]-(onp[1]-onp[0])*0.2,onp[1]+(onp[1]-onp[0])*0.2)
	ax1.set_xlim(taxis[55],taxis[100])
	ax1.set_ylim(DMl,DMh)
	ax1.imshow(TimeArr,aspect='auto',origin='lower',extent=[min(taxis),max(taxis),DMl,DMh],cmap='gray_r',vmin=0.1)
	Y = np.linspace(DMl,DMh,1000)
        X = range(0,len(DMgradArr[0]))
        ax1.contour(taxis,Y,DMgradArr,[0.35,0.5,0.6,0.7,0.8],aspect='auto',origin='lower',extent=[min(taxis),max(taxis),DMl,DMh],cmap='rainbow')
	ax1.set_xlabel('Time (msec)',fontsize=12, fontweight='bold')
	ax1.set_ylabel('DM',fontsize=12, fontweight='bold')
	ax1.tick_params(length=4, width=2)
	ax1.tick_params(axis='x',labelsize=12)
	ax1.tick_params(axis='y',labelsize=12)

	ax2 = plt.subplot2grid((1,4),(0,3), rowspan=1,colspan=1)

	avgDMgrad = np.mean(DMgradArr[:,onp[0]:onp[1]],axis=1)
        avgTime = np.mean(TimeArr[:,onp[0]:onp[1]],axis=1)
	
        avgDMgrad_off = np.mean(DMgradArr_off[:,offp[0]:offp[1]],axis=1)
        avgTime_off = np.mean(TimeArr_off[:,offp[0]:offp[1]],axis=1)

	avgDMgrad_off/=avgDMgrad.max()
	avgTime_off/=avgTime.max()
	
	avgDMgrad/=avgDMgrad.max()
        avgTime/=avgTime.max()
	
	DMgrad_rms = 3*np.std(avgDMgrad_off)
	SNR_rms = 3*np.std(avgTime_off)

	p0 = [1,350,5,1]
		
	print "DM_structure (Gauss fit)"
	popt,cov = curve_fit(gaus,Y,avgDMgrad,p0,sigma=np.full(len(Y),DMgrad_rms),absolute_sigma=True)
	#print popt,cov
	perr = np.sqrt(np.diag(cov))
	#print (popt[0]+popt[3]),3*(perr[0]+perr[3])
        print popt[1] , 3*perr[1]
        #print popt[2] , 3*perr[2]
	
	polyfit_l = 10
	polyfit_h = 990	
	
	pDMgrad = np.polyfit(Y[polyfit_l:polyfit_h],avgDMgrad[polyfit_l:polyfit_h],120)
	pDMgradline = np.poly1d(pDMgrad)
	
	pDMTime =  np.polyfit(Y[polyfit_l:polyfit_h],avgTime[polyfit_l:polyfit_h],120)
	pDMTimeline = np.poly1d(pDMTime)
	
	print "DM_SNR (Gauss fit)"
	popt1,cov1 = curve_fit(gaus,Y,avgTime,p0,sigma=np.full(len(Y),SNR_rms),absolute_sigma=True)
	perr1 = np.sqrt(np.diag(cov1))	
	print popt1[1] , 3*perr1[1]
        #print popt1[2] , 3*perr1[2]
	
	ax2.get_xaxis().set_visible(False)
	ax2.get_yaxis().set_visible(True)
	plt.setp(ax2.get_xticklabels(), visible=True)
	ax2.tick_params(length=4, width=2)
	ax2.yaxis.set_label_position("right")
	ax2.yaxis.tick_right()
	ax2.tick_params(axis='y',labelsize=12)
	ax2.set_ylim(DMl,DMh)
	ax2.plot(avgDMgrad,Y,color='r',linewidth=2)
	ax2.plot(avgTime,Y,color='b',linewidth=2)
	ax2.tick_params(length=4, width=2)

	'''
	plt.figure(2)
	#plt.plot(Y,avgDMgrad,color='r',linewidth=2)
	#plt.plot(Y,pDMgradline(Y),color='g',linewidth=2)
	#plt.plot(Y,avgTime,color='r',linewidth=2)
	#plt.plot(Y,pDMTimeline(Y),color='g',linewidth=2)
	plt.plot(Y,avgDMgrad_off,color='r')
	plt.plot(Y,avgTime_off,color='b')	
	
	print "DM_strcture (polyfit):" 
	print Y[np.argmax(pDMgradline(Y))]
	print "DM_SNR (polyfit):"  
	print Y[np.argmax(pDMTimeline(Y))]		
	'''

	plt.savefig("DM_SNR_struct.pdf", bbox_inches='tight',dpi=300)	
	#plt.plot(Y,avgDMgrad)
	#plt.plot(Y,avgTime)
	plt.show()

				
