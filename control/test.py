#!/usr/bin/env python

import xaosim as xs
from ciao_wfs import WFS
from ciao_wfc import *
import threading
import os
import time

instru = xs.instrument('CIAO')
instru.start()

wfs = WFS()
wfc = TT_WFC()

#os.popen("shmview /tmp/SHcam.im.shm &")
#os.popen("shmview /tmp/xslp.im.shm &")
#os.popen("shmview /tmp/yslp.im.shm &")

time.sleep(5)

# ------------------
#   WFS thread
# ------------------
t = threading.Thread(target=wfs.loop, args=())
t.start() # start the WFS monitoring thread

# ------------------
#  TT calibration
# ------------------

wfs.calc_SH_data(ref=True) # establish current position as reference
wfc.calibrate(0.4)

time.sleep(0.5)


















# ------------------
# the bar is closing
# ------------------
'''

wfs.keepgoing = False
instru.close()

'''
