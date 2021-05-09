import os
import sys

direct = sys.argv[1]

def sch():
	for i in os.listdir(str(direct)):
		if i.endswith('.fits'):
			os.system('pam -f 128 -setnbin 128 -e .fschtsch ' + str(i))
	return

sch()

