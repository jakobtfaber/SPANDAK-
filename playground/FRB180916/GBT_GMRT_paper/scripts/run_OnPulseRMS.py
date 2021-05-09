import numpy as np
import matplotlib.pyplot
import os
import sys

def main():

        #os.system('python OnPulseRMS.py -f A_117/*.norm -o .')
	#os.system('python OnPulseRMS.py -f B_686/*.norm -o .')
	#os.system('python OnPulseRMS.py -f C_1164/*.norm -o .')
	os.system('python OnPulseRMS.py -f D_1549/*.norm -o .')
	os.system('python OnPulseRMS.py -f E_639/*.norm -o .')
	#os.system('python OnPulseRMS.py -f F_267/*.norm -o .')
        os.system('python OnPulseRMS.py -f G_579/*.norm -o .')

main()
