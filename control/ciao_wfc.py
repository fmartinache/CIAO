#!/usr/bin/env python

import numpy as np
from xaosim.shmlib import shm
import sys
import time
import pdb

# =============================================================================
# =============================================================================
def get_coadd(shm1, shm2=None, nav=20, reform=True):
    ''' -----------------------------------------------------------------------
    A useful tool that produces a time average of the last *nav* maps provided
    by one or two *shm* shared memory structure
    ----------------------------------------------------------------------- '''
    try:
        im_cnt1 = shm1.get_counter()
        if shm2 is not None:
            im_cnt2 = shm2.get_counter()
    except:
        print("get_coadd error: arguments is not a valid shm structure?")
        return
    
    myim1   = shm1.get_data(check=im_cnt1, reform=reform)
    im_cnt1 = shm1.get_counter()
    if shm2 is not None:
        myim2   = shm2.get_data(check=im_cnt2, reform=reform)
        im_cnt2 = shm2.get_counter()

    for ii in range(nav-1):
        myim1  += shm1.get_data(check=im_cnt1, reform=reform)
        im_cnt1 = shm1.get_counter()
        if shm2 is not None:
            myim2  += shm2.get_data(check=im_cnt2, reform=reform)
            im_cnt2 = shm2.get_counter()
            
    myim1 /= nav
    if shm2 is not None:
        myim2 /= nav
        return((myim1, myim2))
    else:
        return(myim1)
    
# =============================================================================
# =============================================================================

class DM(object):
    ''' -----------------------------------------------------------------------
    High level information regarding the DM
    ----------------------------------------------------------------------- '''
    def __init__(self,):
        self.dms = 11
        dms = self.dms
        xdm, ydm = np.meshgrid(np.arange(dms)-dms/2, np.arange(dms)-dms/2)
        self.xdm = xdm.T.astype('float32')
        self.ydm = ydm.T.astype('float32')
        self.dmmask   = np.ones((dms, dms), dtype=np.int) # corner mask
        self.dmmask[np.abs(self.xdm) + np.abs(self.ydm) > 7] = 0.0

    def list_2_map(self, data):
        ''' -------------------------------------------------------------------
        Convert a 1D list of 97 voltages back into a 2D map for an
        easier 2D representation.
        ------------------------------------------------------------------- '''
        dms = self.dms
        res = np.zeros((dms, dms), dtype=np.float32)

        i0 = 0
        for k in xrange(dms*dms):
            if (dmmask.flatten()[k] > 0):
                i = k / dms
                j = k % dms
                res[i,j] = data[i0]
                i0 += 1
        return(res)

# =============================================================================
# =============================================================================

class WFC(object):
    ''' -----------------------------------------------------------------------
    Generic wavefront control class
    ----------------------------------------------------------------------- '''
    def __init__(self, shm_cor, shm_cal, nav):
        self.DM = DM()
        self.alpao_cor = shm_cor # correction
        self.alpao_cal = shm_cal # calibration
        self.nav = nav           # number of calibration coadds
        
        self.shm_xslp = shm('/tmp/xslp.im.shm', verbose=False)
        self.shm_yslp = shm('/tmp/yslp.im.shm', verbose=False)
        self.shm_phot = shm('/tmp/phot_inst.im.shm', verbose=False)

        self.modes = np.array([self.DM.xdm, self.DM.ydm]) # default modes (ttilt)
        
    def calibrate(self, a0 = 0.1):
        dm0 = self.alpao_cal.get_data() # DM starting position
        phot = self.shm_phot.get_data()

        RESP = [] # empty response matrix holder
        
        self.nmodes = self.modes.shape[0]
        self.mode_mult = np.zeros(self.nmodes)
        
        # go over the modes to be used for this WFC loop
        for ii in range(self.nmodes):
            self.alpao_cal.set_data((dm0 + self.modes[ii] * a0))
            sys.stdout.write("\rmode %d" % (ii+1,))
            sys.stdout.flush()
            temp = get_coadd(self.shm_xslp, self.shm_yslp,
                             nav=self.nav, reform=False)
            pdb.set_trace()
            
        self.alpao_cal.set_data(dm0) # back to DM starting position
        

    def cloop(self,):
        pass

# =============================================================================
# =============================================================================

class TT_WFC(WFC):
    ''' -----------------------------------------------------------------------
    Tip-tilt wavefront control data structure
    ----------------------------------------------------------------------- '''

    def __init__(self,):
        cor = shm('/tmp/dmdisp2.im.shm', verbose=False) # correction
        cal = shm('/tmp/dmdisp6.im.shm', verbose=False) # calibration
        nav = 5
        super(TT_WFC, self).__init__(cor, cal, nav)

    def calibrate(self, a0 = 0.1):
        super(TT_WFC, self).calibrate(a0=a0)
        
        '''
        dm0 = self.alpao_cal.get_data()
        phot = self.shm_phot.get_data()
        
        self.alpao_cal.set_data(dm0 + self.DM.xdm * a0)
        print("moving X")
        self.tt_y = get_coadd(self.shm_yslp, nav=self.nav, reform=True)
        
        self.alpao_cal.set_data(dm0 + self.DM.ydm * a0)
        print("moving Y")
        self.tt_x = get_coadd(self.shm_xslp, nav=self.nav, reform=True)
        
        self.alpao_cal.set_data(dm0)
        print("Back to base")

        pdb.set_trace()
        self.ttx_mult = a0 / np.mean(self.tt_x[phot > 1.0])
        self.tty_mult = a0 / np.mean(self.tt_y[phot > 1.0])
        '''
        
    def cloop(self,):
        pass

# =============================================================================
# =============================================================================

if __name__ == "__main__":
    wfc = TT_WFC()
