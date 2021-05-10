#!/usr/bin/env python
import matplotlib
matplotlib.use("pdf")
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
import matplotlib.pyplot as plt
from scipy.signal import detrend
import psrchive as psr
from scipy.optimize import curve_fit
from scipy import asarray as ar,exp
import pandas as pd

def gaus(x,a,x0,sigma,c):
    return a*exp(-(x-x0)**2/(2*sigma**2)) + c

if __name__ == "__main__":
	
	filename = sys.argv[1]
	pngfile = filename + ".myrmfit.png"
	rmfile = filename + ".myrm.csv"
	rmfitfile = filename + ".rmfit.csv"	

	#fileorig = psr.Archive_load(filename)
	RM = []
	LvsI = []
	
	for i in range(-500,500,1):
		file1 = psr.Archive_load(filename)	
		#print file1.get_rotation_measure()
		file1.set_rotation_measure(i)
		#print file1.get_rotation_measure()
		file1.defaraday()	
		file1.dedisperse()
		file1.fscrunch()
		file1.remove_baseline()
		data = file1.get_data()
		I = data[0,0,0,:]
		U = data[0,2,0,:] 
		Q = data[0,1,0,:]	
		L = np.sqrt(U**2+Q**2)
		#plt.plot(L)
		#print i,np.max(L)/np.max(I)	
		RM.append(i)
		LvsI.append(np.max(L)/np.max(I))
			
	p0 = [1,-150,100, 0.2]
	RMLvsI = pd.DataFrame({'RM':RM,'LvsI':LvsI,'fname':filename})
	#print RMLvsI
	#popt,cov = curve_fit(gaus,dmsnr['dm'],dmsnr['snr'],p0,bounds=((1,500,5,1),(100,650,100,100)))
	popt,cov = curve_fit(gaus,RM,LvsI,p0)
	perr = np.sqrt(np.diag(cov))
	#print perr

	#result = (popt[0]+popt[3]),3*(perr[0]+perr[3]),popt[1] , 3*perr[1], popt[2], 3*perr[2]
	result = pd.DataFrame({'Amp':(popt[0]+popt[3]),'AmpErr':3*(perr[0]+perr[3]),'Cent':popt[1],'CentErr':3*perr[1],'Wid':popt[2],'WidErr':3*perr[2]},index=[filename])
	#print result 
	#print result
	#print (popt[0]+popt[3]),3*(perr[0]+perr[3]
	#print popt[1] , 3*perr[1]
	#print popt[2] , 3*perr[2]
	plt.xlabel("RM")
	plt.ylabel("L/I")
	plt.scatter(RM,LvsI)
	plt.plot(RM,gaus(RM,*popt),label='fit')
	#plt.show()
	plt.savefig(pngfile)	
	RMLvsI.to_csv(rmfile)
	result.to_csv(rmfitfile)
	#plt.plot(I)
	#plt.show()
