import os
import sys
import numpy as np

calib = sys.argv[1]
fhi = sys.argv[2]
flo = sys.argv[3]
dm = sys.argv[4]
rm = sys.argv[5]
 
def pafit(calib, fhi, flo, dm, rm):
	os.system('source /home/vgajjar/spandakenv/bin/activate')
	print('pam -S -e .4p ' + str(calib))
    	os.system('pam -S -e .4p ' + str(calib))
    	print('pam --DD -e calib.4p.DD ' + str(calib))
    	os.system('pam --DD -e calib.4p.DD ' + str(calib))
   	print('psredit -c "dm=' + str(dm) + '" -e DM ' + str(calib) + '.4p.DD')
    	os.system('psredit -c "dm=' + str(dm) + '" -e DM ' + str(calib) + '.4p.DD')
    	print('pam -D -m ' + str(calib) + '.4p.DM')
    	os.system('pam -D -m ' + str(calib) + '.4p.DM')
    	print('pam -R ' + str(rm) + ' -e RM' + ' ' + str(calib) + '.4p.DM')
    	os.system('pam -R ' + str(rm) + ' -e RM' + ' ' + str(calib) + '.4p.DM')
    	print('paz -F "8200 ' + str(fhi) + '" -F "' + str(flo) + ' 4000" -e RM.zoom ' + str(calib) + '.4p.RM')
    	os.system('paz -F "8200 ' + str(fhi) + '" -F "' + str(flo) + ' 4000" -e RM.zoom ' + str(calib) + '.4p.RM')
    	print('pdv -FT -Z -t ' + str(calib) + '.4p.RM.zoom >' + str(calib) + '.RM.zoom.paprof')
    	os.system('pdv -FT -Z -t ' + str(calib) + '.4p.RM.zoom >' + str(calib) + '.RM.zoom.paprof')
    	print('python ../../../../rm_pa_fit/PlotFigure2.py -f ' + str(calib) + '.4p --pol -R ' + str(rm) + ' --fl ' + str(flo) + ' --fh ' + str(fhi) + ' --ext 10 --PAF ' + str(calib) + '.RM.zoom.paprof')
    	os.system('python ../../../../rm_pa_fit/PlotFigure2.py -f ' + str(calib) + '.4p --pol -R ' + str(rm) + ' --fl ' + str(flo) + ' --fh ' + str(fhi) + ' --ext 10 --PAF ' + str(calib) + '.RM.zoom.paprof')

	return

pafit(calib, fhi, flo, dm, rm)
