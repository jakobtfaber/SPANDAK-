import os
import sys
import numpy as np

#calib = sys.argv[1]
pre_calib = '.pazi'
calib = '.tsch256'
#fhi = sys.argv[1]
fhi = 999
#flo = sys.argv[2]
flo = 601
dm = sys.argv[1]
rm = sys.argv[2]
 
def pafit(pre_calib, calib, fhi, flo, dm, rm, ts=256):
	os.system('pam --setnbin ' + str(ts) + ' -e .tsch256 *' + str(pre_calib))
	print('TScrunched to 256 bins')
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
    	print('paz -F "1000 ' + str(fhi) + '" -F "' + str(flo) + ' 600" -e RM.zoom ' + '*' + '.4p.RM')
    	os.system('paz -F "1000 ' + str(fhi) + '" -F "' + str(flo) + ' 600" -e RM.zoom ' + '*' + '.4p.RM')
    	print('pdv -FT -Z -t ' +  '*' +  '.4p.RM.zoom > ' + 'pulse.RM.zoom.paprof')
    	os.system('pdv -FT -Z -t ' + '*' + '.4p.RM.zoom > ' + 'pulse.RM.zoom.paprof')
    	print('python ../../../../../rm_pa_fit/PlotFigure2.py -f ' + '*' +  '.4p --pol -R ' + str(rm) + ' --fl ' + str(flo) + ' --fh ' + str(fhi) + ' --ext 10 --PAF ' + '*' + '.RM.zoom.paprof')
    	os.system('python ../../../../../rm_pa_fit_2/PlotFigure2.py -f ' + '*' + '.4p --pol -R ' + str(rm) + ' --fl ' + str(flo) + ' --fh ' + str(fhi) + ' --ext 10 --PAF ' + '*' + '.RM.zoom.paprof')

	return

pafit(pre_calib, calib, fhi, flo, dm, rm)
