#!/usr/bin/env python

import xaosim as xs
from ciao_wfs import WFS
from ciao_wfc import *
import threading
import os
import time

#instru = xs.instrument('CIAO')
#instru.start(delay=0.05)

wfs = WFS(shmf="/tmp/ixon.im.shm")
wfc = ZER_WFC(4, 15)

#os.popen("shmview /tmp/ixon.im.shm &")
#os.popen("shmview /tmp/xslp.im.shm &")
#os.popen("shmview /tmp/yslp.im.shm &")
time.sleep(5)

def close():
    wfs.keepgoing = False
    wfc.keepgoing = False
    #instru.close()
    
# ------------------
#   WFS thread
# ------------------
t = threading.Thread(target=wfs.loop, args=())
t.start() # start the WFS monitoring thread

# ------------------
#  calibration
# ------------------

wfs.calc_SH_data(ref=True) # establish current position as reference
wfc.nav = 50
wfc.calibrate(a0=0.02, reform=False)

# ------------------
#   WFC thread
# ------------------
wfc.gain = 0.001
#t1 = threading.Thread(target=wfc.cloop, args=())
#t1.start() # start the WFS monitoring thread
#wfc.cloop()


'''

close()

'''
