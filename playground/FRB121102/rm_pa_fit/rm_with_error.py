import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys

directory = sys.argv[1]

def plot_err(directory):
	rmLI = pd.read_csv(str(directory), sep=',', usecols=['LvsI','RM'])
	LI = rmLI['LvsI']
	rm = rmLI['RM']
	fig = plt.figure()
	plt.scatter(rm, LI, 
	
	
