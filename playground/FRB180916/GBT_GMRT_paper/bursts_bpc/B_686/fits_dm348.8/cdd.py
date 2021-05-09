import os
os.system('dspsr -U 1024 -F 7680:D -K -d 4 -b 4096 -E  ../par/*.par -s -a psrfits -e fits ../raw/*.raw')
