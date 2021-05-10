#!/usr/bin/env python
import matplotlib
matplotlib.use('pdf')
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
from Coherent_dedisp import *

if __name__ == "__main__":
	
	parser = optparse.OptionParser()

        parser.add_option("-f", action='store', dest='infile', type=str, help="File with list of RAW files")

        parser.add_option("-o", action='store', dest='outdir', default="",type=str,help="Full Output directory (Default : .)")

	parser.add_option("-R", action='store', dest='RM', default="93500",type=float,help="RM for the burst (Default : 93500)")

        parser.add_option("-D", action='store', dest='DM', default=565, type=float,help="DM")

        parser.add_option("--fl", action='store', dest='freql', default=4000, type=float,help="Lower frequency for the pulse to add (this only affect the average pulse)")

        parser.add_option("--fh", action='store', dest='freqh', default=8000, type=float,help="Higher frequency for the pulse to add (this only affect the average pulse)")

	parser.add_option("--tl", action='store', dest='low', default=0, type=int,help="Lower time bin for the pulse to add ")

        parser.add_option("--th", action='store', dest='high', default=127, type=int,help="Higher time bin for the pulse to add")

	parser.add_option("-s", action='store', dest='snr', default=6, type=float,help="SNR threshold to detect pulse")

	parser.add_option("--ext", action='store', dest='ext', default=1, type=int,help="Extra phase bins")

        parser.add_option("--pol", action='store_true', dest='polplot',help="Polarization plot")

	parser.add_option("--PAF", action='store', default="",dest='PAfile',help="PA from PSRCHIVE as an ASCII file generated using 'pdv -FT -Z -t'")

	parser.add_option("--xlab", action='store_true', dest='minlabx',help="X lab")
	parser.add_option("--ylab", action='store_true', dest='minlaby',help="Y lab")
	
        options,args = parser.parse_args()

        infile = os.path.abspath(options.infile)
        freql = options.freql
        freqh = options.freqh
	low = options.low
	high = options.high
	snr = options.snr
	RM = options.RM
	DM = options.DM
	polplot = options.polplot
	ext = options.ext
	PAfile = options.PAfile
	minlabx = options.minlabx
	minlaby = options.minlaby

	if not options.outdir: outdir = os.getcwd()
        else:
                outdir = options.outdir
                if(os.path.isdir(outdir) is not True):
                        os.system("mkdir %s" % (outdir))

	#print "Plotting " + infile
	
	plotfig(infile,freql,freqh,snr,low,high,RM,polplot,ext,DM,PAfile,minlabx,minlaby)	

