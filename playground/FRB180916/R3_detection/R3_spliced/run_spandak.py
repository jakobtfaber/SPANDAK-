import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import pandas as pd

#Read in -tentative 'database'- .csv
filepaths = pd.read_csv("database_r3_spliced.csv")

filpaths = filepaths.iloc[:, 0]

print(filpaths)

def run_spandak(filpaths, hidm=2000, lodm=100):

	run_commands = []
	
	#Define csv name by fil path range
	#For running R3 blc73/74 individually
	#r = -17
	#l = 50
	#For running R3 spliced
	r = -17 
	l = 59
	#For running 121102
	#r = -25
	#l = 40

	for fil in filpaths:
		spandak_run = 'SPANDAK ' + '--fil ' + fil + ' --hidm ' + str(hidm) + ' --lodm ' + str(lodm)  + ' ' +  '--logs=' + fil[l:r] + '.csv'	
		run_commands.append(spandak_run)

	return run_commands


def main():	

	run_commands = run_spandak(filpaths)

	for run in run_commands:
		#print('SPANDAK run ', run)
		os.system(run)

main()
