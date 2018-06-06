#!/usr/bin/env python

''' ===========================================================================
Module to be imported by the CIAO_MASTER control GUI of CIAO for a more
flexible use of the DM.

WORK IN PROGRESS !!! CURRENTLY INCOMPLETE 

=========================================================================== '''

import numpy as np
import os, sys, time
import struct

from xaosim.shmlib import shm

if (8 * struct.calcsize("P")) == 32:
    print("Use x86 libraries.")
    sys.path.append(home+"/python/libs/alpao/Lib/")
else:
    print("Use x86_64 libraries.")
sys.path.append(home+"/python/libs/alpao/Lib64/")

from asdk import DM

# =============================================================================
# =============================================================================

class ALPAO_DM():
    # -------------------------------------------------------------------------
    def __init__(self,):
        self.seri_num = "BOL115"
        self.alpao   = DM(self.seri_num)
        self.nba     = int(alpao.Get('NBOfActuator'))
        self.alpao.Send([0.] * nbAct)  # init the DM: zeroes everywhere
        self.dms     = 11              # DM size
        dms = self.dms
        xdm, ydm     = np.meshgrid(np.arange(dms)-dms/2, np.arange(dms)-dms/2)
        self.xdm     = xdm.T.astype('float32')
        self.ydm     = ydm.T.astype('float32')
        self.dmmask  = np.ones((dms, dms), dtype=np.int) # corner mask
        self.dmmask[np.abs(self.xdm) + np.abs(self.ydm) > 7] = 0.0
        self.flat_mask = self.dmmask.flatten()

        self.nch   = 8 # number of DM channels
        self.dmd0  = np.zeros((dms,dms), dtype=np.float32)
        self.cntrs = np.zeros(self.nch+1) - 1 # shared mem counters for channels

        for ii in range(self.nch):
            exec "disp%d = shm(fname='/tmp/dmdisp%d.im.shm', data=self.dmd0, verbose=False)" % (ii,ii)

    # -------------------------------------------------------------------------
    def map2dm(self, cmap):
        ''' --------------------------------------------
        Sends a 2D command to the DM, using the mapping
        knowledge provided by the flat_mask array.
        -------------------------------------------- ''' 
        temp   = cmap.flatten()
        values = [0.0] * self.nba # data to be uploaded to the DM
        
        jj = 0
        for ii in xrange(cmap.size):
            if not np.isnan(self.flat_mask[ii]):
                values[jj] = temp[i]
                jj += 1
        self.alpao.Send(values)

    # -------------------------------------------------------------------------
    def reset(self, chn=0):
        ''' ----------------------------------------------
        Reset one or (default) all channels of the DM.
        ---------------------------------------------- '''
        if chn == 0: # send the flat-map!
            try:
                data = np.loadtxt('/home/ciaodev/.config/ciao/zygo_flat.txt')
                disp0.set_data(list2map2D(data))
            except:
                disp0.set_data(self.dmd0)
        else:
            try:
                exec 'self.dmdisp%d.set_data(self.dmd0)' % 
            except:
                print("error reset channel #%d" % (chn))

    # -------------------------------------------------------------------------
    def quit(self, chn=None):

    # -------------------------------------------------------------------------
    def quit(self, chn=None):
        self.alpao.Reset()
        sys.exit()

# =============================================================================
# =============================================================================

if __name__ == "__main__":
    mydm = ALPAO_DM()
