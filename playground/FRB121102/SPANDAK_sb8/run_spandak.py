import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import pandas as pd

#Read in -tentative 'database'- .csv

def run_spandak(filpaths, hidm=2000, lodm=100):

	run_commands = []
	
	#Define csv name by fil path range
	#For running R3 blc73/74 individually
	#r = -17
	#l = 109
	#For running R3 spliced
	#r
	#l
	#For running 121102
	r = -25
	l = 40

	for fil in filpaths:
		spandak_run = 'SPANDAK ' + '--fil ' + fil + ' --hidm ' + str(hidm) + ' --lodm ' + str(lodm) + ' --dorfi ' + ' --subBand 8 ' +  '--logs=' + fil[l:r] + '.csv'	
		run_commands.append(spandak_run)

	return run_commands


if __name__ == '__main__':	
	
	database = sys.argv[1]
	filepaths = pd.read_csv(str(database))
	filpaths = filepaths.iloc[:, 0]
	run_commands = run_spandak(filpaths)

	for run in run_commands:
		#print('SPANDAK run ', run)
		os.system(run)
		os.system('rm *downsampled*')

