import os
import sys
import numpy as np

#calib = sys.argv[1]
#pre_calib = '.zoom'
calib = sys.argv[1]
#fhi = sys.argv[1]
fhi = 8199
#flo = sys.argv[2]
flo = 4401
dm = sys.argv[2]
rm = sys.argv[3]
 
def pafit(calib, fhi, flo, dm, rm, ts=128):
	#os.system('pam --setnbin ' + str(ts) + ' -e .tsch128 *' + str(pre_calib))
	#print('TScrunched to 128 bins')
	print('pam -S -e .4p ' + '*' + str(calib))
    	os.system('pam -S -e .4p ' + '*' + str(calib))
    	print('pam --DD -e .4p.DD ' + '*' + str(calib))
    	os.system('pam --DD -e .4p.DD ' + '*' + str(calib))
   	print('psredit -c "dm=' + str(dm) + '" -e DM ' + '*' + '.4p.DD')
    	os.system('psredit -c "dm=' + str(dm) + '" -e DM ' + '*' +  '.4p.DD')
    	print('pam -D -m ' + '*' +  '.4p.DM')
    	os.system('pam -D -m ' + '*' +  '.4p.DM')
    	print('pam -R ' + str(rm) + ' -e RM' + ' ' + '*' + '.4p.DM')
    	os.system('pam -R ' + str(rm) + ' -e RM' + ' ' + '*' + '.4p.DM')
    	print('paz -F "8200 ' + str(fhi) + '" -F "' + str(flo) + ' 4400" -e RM.zoom ' + '*' + '.4p.RM')
    	os.system('paz -F "8200 ' + str(fhi) + '" -F "' + str(flo) + ' 4400" -e RM.zoom ' + '*' + '.4p.RM')
    	print('pdv -FT -Z -t ' +  '*' +  '.4p.RM.zoom > ' + 'pulse.RM.zoom.paprof')
    	os.system('pdv -FT -Z -t ' + '*' + '.4p.RM.zoom > ' + 'pulse.RM.zoom.paprof')
    	print('python ../../rm_pa_fit/PlotFigure2.py -f ' + '*' +  '.4p --pol -R ' + str(rm) + ' --fl ' + str(flo) + ' --fh ' + str(fhi) + ' --ext 10 --PAF ' + '*' + '.RM.zoom.paprof')
    	os.system('python /datax/scratch/jfaber/FLITS/playground/FRB121102/rm_pa_fit_2/PlotFigure2.py -f ' + '*' + '.4p --pol -R ' + str(rm) + ' --fl ' + str(flo) + ' --fh ' + str(fhi) + ' --ext 10 --PAF ' + '*' + '.RM.zoom.paprof')

	return

pafit(calib, fhi, flo, dm, rm)
