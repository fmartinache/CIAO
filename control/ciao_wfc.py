#!/usr/bin/env python

import numpy as np
from xaosim.shmlib import shm
from xaosim.zernike import mkzer1
import astropy.io.fits as pf
from numpy.linalg import solve
from numpy.linalg import pinv

import sys
import time

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
        self.nba = 97  # number of actuators for the DM 97-15
        self.dms = 11  # equivalent grid size of the DM 97-15 (11x11 - corners)
        dms = self.dms
        xdm, ydm = np.meshgrid(np.arange(dms)-dms/2, np.arange(dms)-dms/2)
        self.xdm = xdm.T.astype('float32')
        self.ydm = ydm.T.astype('float32')
        self.dmmask   = np.ones((dms, dms), dtype=np.int) # corner mask
        self.dmmask[np.abs(self.xdm) + np.abs(self.ydm) > 7] = 0.0

    # -------------------------------------------------------------------------
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
    
    # -------------------------------------------------------------------------
    def zer_mode_bank_2D(self, i0, i1):
        ''' -------------------------------------------------------------------
        Produce a bank of Zernike modes for the ALPAO 97-15 DM.

        Parameters:
        ----------
        - i0: index of the first zernike mode to be added
        - i1: index of the last zernike mode to be included
        ------------------------------------------------------------------- '''
        dZ = 1 + i1 - i0
        res = np.zeros((dZ, self.dms, self.dms)).astype('float32')
        for i in range(i0, i1+1):
            res[i-i0] = mkzer1(i, self.dms, self.dms/2 +1) * self.dmmask
        return(res)
    
    # -------------------------------------------------------------------------
    def poke_mode_bank_2D(self):
        ''' -------------------------------------------------------------------
        Produce a bank of poke modes for the ALPAO 97-15 DM
        ------------------------------------------------------------------- '''
        res = np.zeros((self.nba, self.dms, self.dms)).astype('float32')

        kk = 0
        for ii in range(self.dms):
            for jj in range(self.dms):
                if self.dmmask[jj,ii] > 0.5:
                    res[kk, jj,ii] = 1.0
                    kk += 1
        return res
    
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

        self.shm_comb = shm('/tmp/comb.im.shm', verbose=False)
        
        self.shm_xslp = shm('/tmp/xslp.im.shm', verbose=False)
        self.shm_yslp = shm('/tmp/yslp.im.shm', verbose=False)
        self.shm_phot = shm('/tmp/phot_inst.im.shm', verbose=False)
        
        self.modes     = np.array([self.DM.xdm, self.DM.ydm]) # (ttilt)
        self.keepgoing = False
        self.gain      = 0.001 # default gain
        self.tsleep    = 1e-6 # for tests
        self.verbose   = False
        self.abort     = False
        self.calib_on  = False
        
    # -------------------------------------------------------------------------
    def get_slopes_comb(self, nav=20, reform=True):
        ''' -------------------------------------------------------------------
        test
        ------------------------------------------------------------------- '''
        cnt = self.shm_comb.get_counter()
        tmp = self.shm_comb.get_data(check=cnt, reform=True)
        cnt = self.shm_comb.get_counter()
        
        x_sig = tmp[0]
        y_sig = tmp[1]
        
        for ii in range(nav-1):
            tmp = self.shm_comb.get_data(check=cnt, reform=True)
            cnt = self.shm_comb.get_counter()
            x_sig += tmp[0]
            y_sig += tmp[1]
        x_sig /= nav
        y_sig /= nav

        return np.concatenate((x_sig.flatten(), y_sig.flatten()))
    
    # -------------------------------------------------------------------------
    def get_slopes(self, nav=20, reform=True):
        ''' -------------------------------------------------------------------
        test
        ------------------------------------------------------------------- '''
        x_cnt = self.shm_xslp.get_counter()
        x_sig = self.shm_xslp.get_data(check=x_cnt, reform=reform)
        x_cnt = self.shm_xslp.get_counter()
        
        y_cnt = self.shm_yslp.get_counter()
        y_sig = self.shm_yslp.get_data(check=y_cnt, reform=reform)
        y_cnt = self.shm_yslp.get_counter()
        
        for ii in range(nav-1):
            x_sig += self.shm_xslp.get_data(check=x_cnt, reform=reform)
            x_cnt  = self.shm_xslp.get_counter()
            
            y_sig += self.shm_yslp.get_data(check=y_cnt, reform=reform)
            y_cnt  = self.shm_yslp.get_counter()
            
        x_sig /= nav
        y_sig /= nav

        return np.concatenate((x_sig, y_sig))

    # -------------------------------------------------------------------------
    def reset(self,):
        self.alpao_cor.set_data(0.0 * self.alpao_cor.get_data())
        
    # -------------------------------------------------------------------------
    def calibrate(self, a0 = 0.1, reform=False):
        ''' -------------------------------------------------------------------
        Generic calibration procedure for set of modes attached to the object
        ------------------------------------------------------------------- '''
        self.calib_on = True
        dm0 = self.alpao_cal.get_data() # DM starting position
        #phot = self.shm_phot.get_data()
        phot = self.shm_comb.get_data()[2]
        
        RESP = [] # empty response matrix holder
        
        self.nmodes = self.modes.shape[0]
        self.mode_mult = np.zeros(self.nmodes)
        
        # go over the modes to be used for this WFC loop
        for ii in range(self.nmodes):
            if self.abort:
                self.abort = False
                return
            self.alpao_cal.set_data((dm0 + self.modes[ii] * a0))
            time.sleep(self.tsleep) # for simulation only!
            sys.stdout.write("\rmode %d" % (ii+1,))
            sys.stdout.flush()

            RESP.append(self.get_slopes_comb(self.nav, reform=reform))

        self.alpao_cal.set_data(dm0) # back to DM starting position
        self.RR = np.array(RESP) / a0
        self.RR[np.abs(self.RR) < 1e-2] = 0.0 # matrix clean up
        self.RTR = self.RR.dot(self.RR.T)
        self.calib_on = False
    # -------------------------------------------------------------------------
    def cloop(self,):
        ''' -------------------------------------------------------------------
        Generic closed-loop procedure for the modes attached to the object.
        ------------------------------------------------------------------- '''

        self.keepgoing = True
        while self.keepgoing:
            sig = self.get_slopes_comb(1, reform=False)
            dm0 = self.alpao_cor.get_data()
            #ww  = solve(self.RTR, np.dot(self.RR, sig)) + 1e-9
            ww = self.RRinv.dot(sig) + 1e-9
            tmp = np.average(self.modes, weights=ww, axis=0)
            cor = tmp * self.nmodes * ww.sum()
            dm1 = 0.999 * (dm0 - self.gain * cor.astype('float32'))
            self.alpao_cor.set_data(dm1)

            time.sleep(self.tsleep)

            if self.verbose:
                sys.stdout.write("\rcoeffs: "+ "%+.3f " * self.nmodes % (tuple(ww)))
                sys.stdout.flush()
        print("\nWFC loop opened.")
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
        self.modes = np.array([self.DM.xdm, self.DM.ydm]) # (ttilt)
        self.verbose = True

    def calibrate(self, a0 = 0.1, reform=False):
        super(TT_WFC, self).calibrate(a0=a0, reform=reform)
        self.RRinv = pinv(self.RR.T, rcond=0.1)
        
    def cloop(self,):
        super(TT_WFC, self).cloop()

    def reset(self,):
        super(TT_WFC, self).reset()
        
# =============================================================================
# =============================================================================

class ZER_WFC(WFC):
    ''' -----------------------------------------------------------------------
    Zernike based wavefront control data structure.
    ----------------------------------------------------------------------- '''

    def __init__(self, iz0=4, iz1=11):
        ''' -------------------------------------------------------------------
        Class constructor.

        Specifies calibration and correction channels as well as the basis
        of modes to use to control the wavefront.
        
        Parameters:
        ----------
        - iz0: first index of the Zernike modes to control (default = 4)
        - iz1: last  index of the Zernike modes to control (default = 11)
        ------------------------------------------------------------------- '''
        
        cor = shm('/tmp/dmdisp3.im.shm', verbose=False) # correction
        cal = shm('/tmp/dmdisp7.im.shm', verbose=False) # calibration
        nav = 5
        super(ZER_WFC, self).__init__(cor, cal, nav)
        self.modes = self.DM.zer_mode_bank_2D(iz0, iz1)
        self.verbose = True
        
    def calibrate(self, a0 = 0.1, reform=False):
        super(ZER_WFC, self).calibrate(a0=a0, reform=reform)
        self.RRinv = pinv(self.RR.T, rcond=0.1)
        
    def cloop(self,):
        super(ZER_WFC, self).cloop()

    def reset(self,):
        super(ZER_WFC, self).reset()

# =============================================================================
# =============================================================================

class ZON_WFC(WFC):
    ''' -----------------------------------------------------------------------
    Zonal based wavefront control data structure.
    ----------------------------------------------------------------------- '''

    def __init__(self,):
        ''' -------------------------------------------------------------------
        Class constructor.

        Specifies calibration and correction channels as well as the basis
        of modes to use to control the wavefront.
        
        ------------------------------------------------------------------- '''
        
        cor = shm('/tmp/dmdisp3.im.shm', verbose=False) # correction
        cal = shm('/tmp/dmdisp7.im.shm', verbose=False) # calibration
        nav = 5
        super(ZON_WFC, self).__init__(cor, cal, nav)
        self.modes = self.DM.poke_mode_bank_2D()
        self.verbose = False

    def calibrate(self, a0 = 0.1, reform=False):

        super(ZON_WFC, self).calibrate(a0=a0, reform=reform)
        # filtering for the pseudo-inverse matrix
        self.RRinv = pinv(self.RR.T, rcond=0.1)
        
    def cloop(self,):
        super(ZON_WFC, self).cloop()

    def reset(self,):
        super(ZON_WFC, self).reset()

# =============================================================================
# =============================================================================

if __name__ == "__main__":
    wfc = TT_WFC()
