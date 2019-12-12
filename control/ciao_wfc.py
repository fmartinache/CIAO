#!/usr/bin/env python

import numpy as np
from xaosim.shmlib import shm
from xaosim.zernike import mkzer1
from numpy.linalg import solve
from numpy.linalg import pinv
import astropy.io.fits as pf

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
            if (self.dmmask.flatten()[k] > 0):
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
            res[i-i0] = mkzer1(i, self.dms, self.dms/2 +2) * self.dmmask
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
        self.a0 = 0.1            # calibration amplitude
        
        self.shm_comb = shm('/tmp/comb.im.shm', verbose=False)
        
        self.modes     = np.array([self.DM.xdm, self.DM.ydm]) # (ttilt)
        self.keepgoing = False
        self.gain      = 0.001 # default gain
        self.tsleep    = 1e-6 # for tests
        self.verbose   = False
        self.abort     = False
        self.calib_on  = False
        self.loop_on   = False
        self.avgtiming = 0.
        self.timingT   = 100
        
    # -------------------------------------------------------------------------
    def get_slopes(self, nav=20, reform=True):
        ''' -------------------------------------------------------------------
        test
        ------------------------------------------------------------------- '''
        cnt = self.shm_comb.get_counter()
	#time_start = time.time()
        tmp = self.shm_comb.get_data(check=cnt, reform=True)
	#print("Timing call get_data: %.3e" % (time.time() - time_start))
        cnt = self.shm_comb.get_counter()
        
        x_sig = tmp[0]
        y_sig = tmp[1]
        phot  = tmp[2]

       	for ii in range(nav-1):
            tmp = self.shm_comb.get_data(check=cnt, reform=True)
            cnt = self.shm_comb.get_counter()
            x_sig += tmp[0]
            y_sig += tmp[1]
            phot  += tmp[2]
        x_sig /= nav
        y_sig /= nav
        phot /= nav

        x_sig[phot == 0] = 0.0
        y_sig[phot == 0] = 0.0
        return np.concatenate((x_sig.flatten(), y_sig.flatten()))
    
    # -------------------------------------------------------------------------
    def reset(self,):
        self.alpao_cor.set_data(0.0 * self.alpao_cor.get_data())
        
    # -------------------------------------------------------------------------
    def calibrate(self, a0 = 0.1, reform=False):
        ''' -------------------------------------------------------------------
        Generic calibration procedure for set of modes attached to the object
        ------------------------------------------------------------------- '''
        self.calib_on = True
        self.abort = False
        self.a0 = a0
        
        dm0 = self.alpao_cal.get_data() # DM starting position

        phot = self.shm_comb.get_data()[2]
        
        RESP = [] # empty response matrix holder
        
        self.nmodes = self.modes.shape[0]
        self.mode_mult = np.zeros(self.nmodes)
        
        # go over the modes to be used for this WFC loop
        for ii in range(self.nmodes):
            if self.abort:
                self.abort = False
                self.calib_on = False
                return
            self.alpao_cal.set_data((dm0 + self.modes[ii] * a0))
            #time.sleep(self.tsleep) # for simulation only!
            sys.stdout.write("\rmode %d" % (ii+1,))
            sys.stdout.flush()

            RESP.append(self.get_slopes(self.nav, reform=reform))

        self.alpao_cal.set_data(dm0) # back to DM starting position
        self.RR = np.array(RESP)
        #self.RR[np.abs(self.RR) < 1e-2] = 0.0 # matrix clean up
        self.RTR = self.RR.dot(self.RR.T)
        self.calib_on = False

    # -------------------------------------------------------------------------
    def reload_cal(self, fname=None):
        ''' -------------------------------------------------------------------
        Tries to reload a previously saved calibration?
        ------------------------------------------------------------------- '''

        try:
            self.RR = pf.getdata(fname)
        except:
            print("Calibration matrix %s not available" % (fname))
            return
        self.RRinv = pinv(self.RR.T, rcond=0.1)

    # -------------------------------------------------------------------------
    def cloop(self,):
        ''' -------------------------------------------------------------------
        Generic closed-loop procedure for the modes attached to the object.
        ------------------------------------------------------------------- '''

        self.keepgoing = True
        self.loop_on = True
        self.avgtiming = 0.
        counter = 0
        time_start = time.time()
        time_end = time_start
        time_A = 0.
        time_B = 0.
        time_C = 0.
        time_D = 0.
        time_E = 0.
        time_F = 0.
        time_G = 0.
	
	time_tmp = time_start
        diff_time = 0.
        while self.keepgoing:
	    time_tmp = time.time()
            sig = self.get_slopes(1, reform=False) # WFS slopes
	    time_A += time.time() - time_tmp
            
	    time_tmp = time.time()
	    dm0 = self.alpao_cor.get_data()        # DM shape B4 correction  
	    time_B += time.time() - time_tmp

	    time_tmp = time.time()
            ee  = self.a0 * self.RRinv.dot(sig)    # error signal
	    time_C += time.time() - time_tmp

	    time_tmp = time.time()
            cor = np.average(self.modes, weights=ee, axis=0) 
	    time_D += time.time() - time_tmp

	    time_tmp = time.time()
            cor *= self.nmodes * ee.sum()
	    time_E += time.time() - time_tmp

	    time_tmp = time.time()
            dm1 = 0.999 * (dm0 - self.gain * cor.astype('float32'))
	    time_F += time.time() - time_tmp
            
 	    time_tmp = time.time()
            self.alpao_cor.set_data(dm1)
	    time_G += time.time() - time_tmp

            '''time.sleep(self.tsleep)'''

            counter += 1
            if counter % self.timingT == 0:
                time_end = time.time()
                self.avgtiming = (time_end - time_start)/counter
                '''print("cloop average timing: %.3f %.3f %.3f %.3f %.3f %.3f %.3f"% (self.avgtiming, 1./self.avgtiming))'''
                print("cloop average timing: %.3f %.3f %.3f %.3f %.3f %.3f %.3f"% (time_A/counter,time_B/counter,time_C/counter,time_D/counter,time_E/counter,time_F/counter,time_G/counter))
                time_start = time_end
                counter = 0
        	time_A = 0.
        	time_B = 0.
        	time_C = 0.
        	time_D = 0.
        	time_E = 0.
        	time_F = 0.
        	time_G = 0.
		

            if self.verbose:
                sys.stdout.write(
                    "\rcoeffs: "+ "%+.3f " * self.nmodes % (tuple(ee)))
                sys.stdout.flush()
        self.loop_on = False
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

    def reload_cal(self, fname=None):
        super(TT_WFC, self).reload_cal(fname=fname)

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

    def reload_cal(self, fname=None):
        super(ZER_WFC, self).reload_cal(fname=fname)

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
        
    def reload_cal(self, fname=None):
        super(ZON_WFC, self).reload_cal(fname=fname)

# =============================================================================
# =============================================================================

if __name__ == "__main__":
    wfc = TT_WFC()
