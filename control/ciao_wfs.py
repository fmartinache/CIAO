#!/usr/bin/env python

import sys
import numpy as np
from xaosim.shmlib import shm
import datetime

def centroid_position_0(img):
    ''' --------------------------------------------------------------
    Returns the (x,y) position of the centroid of the array *img*
    -------------------------------------------------------------- '''
    xxc = np.arange(img.shape[1])
    yyc = np.arange(img.shape[0])
    mprofx = img.mean(0)
    mprofy = img.mean(1)
    mprofx -= mprofx.min()
    mprofy -= mprofy.min()
    denomx = np.sum(mprofx)
    denomy = np.sum(mprofy)

    if denomx > 1:
        xc = np.sum(mprofx * xxc) / denomx
    else:
        xc = 0.0
    if denomy > 1:
        yc = np.sum(mprofy * yyc) / denomy
    else:
        yc = 0.0
    return((yc, xc))

def centroid_position(img):
    ''' --------------------------------------------------------------
    Returns the (x,y) position of the centroid of the array *img*
    -------------------------------------------------------------- '''
    return np.unravel_index(np.argmax(img, axis=None), img.shape)

# =============================================================================
# =============================================================================

class WFS():
    # =========================================================
    def __init__(self, shmf="/tmp/SHcam.im.shm"):

        # -------------------------
        # to be provided externally
        # -------------------------
        self.SH_x0, self.SH_y0 = 6, 0   # properties of the grid
        self.SH_dx, self.SH_dy = 12.8, 12.8 # properties of the grid
        self.threshold = 110.0
        
        # -----------------------------
        # shared memory data structures
        # -----------------------------
        self.shm_im = shm(shmf)                  # shared mem. image
        self.isz = self.shm_im.mtdata['size'][0] # should be 128
        self.im_cnt = self.shm_im.get_counter()  # image counter
        self.update_cells()
        self.define_SH_data()

        self.ttx_mean = 0.0
        self.tty_mean = 0.0
        self.log_len = 200 # length of the tip-tilt log
        self.ttx_log = []
        self.tty_log = []
        # ------------------
        # control structures
        # ------------------
        self.keepgoing = False

    # =========================================================
    def update_cells(self,):
        ''' -------------------------------------------------------
        Updates the SH cells positions for new SH grid parameters
        ------------------------------------------------------- '''
        nx = int((self.isz - self.SH_x0) / self.SH_dx + 2)
        xx = np.round(self.SH_x0+np.arange(nx-1)*self.SH_dx).astype(int)
        if (self.isz - xx[-1] > self.SH_dx/2+2):
            xx = np.append(xx, self.isz)
        self.SH_xx = xx
        
        ny = int((self.isz - self.SH_y0) / self.SH_dy + 2)
        yy = np.round(self.SH_y0+np.arange(ny-1)*self.SH_dy).astype(int)
        if (self.isz - yy[-1] > self.SH_dy/2+2):
            yy = np.append(yy, self.isz)
        self.SH_yy = yy

    # =========================================================
    def define_SH_data(self):
        xx, yy   = self.SH_xx , self.SH_yy        
        ncx, ncy = xx.size - 1, yy.size - 1

        self.SH_phot = np.zeros((ncy, ncx)) # photometry
        self.SH_xslp = np.zeros((ncy, ncx)) # x-slope
        self.SH_yslp = np.zeros((ncy, ncx)) # y-slope
        self.SH_xref = np.zeros((ncy, ncx)) # x-reference
        self.SH_yref = np.zeros((ncy, ncx)) # y-reference
        self.SH_comb = np.zeros((3, ncy, ncx)) # xy slopes + phot combined
        # additional arrays for convenience
        self.SH_xtmp = np.zeros((ncy, ncx))
        self.SH_ytmp = np.zeros((ncy, ncx))

        self.shm_SNR = shm('/tmp/SNR.im.shm', data=self.SH_phot, 
                           verbose=False) # experiment!

        self.shm_phot_inst = shm('/tmp/phot_inst.im.shm', data=self.SH_phot, 
                                 verbose=False) # experiment!

        #self.shm_xslp = shm('/tmp/xslp.im.shm',data=self.SH_xslp,verbose=False)
        #self.shm_yslp = shm('/tmp/yslp.im.shm',data=self.SH_yslp,verbose=False)
        self.shm_comb = shm('/tmp/comb.im.shm',data=self.SH_comb,verbose=False)
        
        for j in xrange(ncy):
            y0, y1 = int(np.round(yy[j])), int(np.round(yy[j+1]))

            for i in xrange(ncx):
                x0, x1 = int(np.round(xx[i])), int(np.round(xx[i+1]))

                self.SH_xref[j,i] = 0.5 * (x1 - x0)
                self.SH_yref[j,i] = 0.5 * (y1 - y0)

    # =========================================================
    def update_grid(self, **kwargs):
        ''' -------------------------------------------------------
        update the parameters that define the SH cells
        ------------------------------------------------------- '''
        valid_keys = {"x0": "SH_x0", "y0": "SH_y0",
                      "dx": "SH_dx", "dy": "SH_dy",
                      "i0": "SH_threshold"}
        if kwargs is not None:
            for key, value in kwargs.items():
                if key in valid_keys:
                    setattr(self, valid_keys[key], value)
                else:
                    print("%s: invalid key" % (key,))
                    print("valid keys are: " +
                          len(valid_keys.keys()) * "%s " % \
                          tuple(valid_keys.keys()))
        else:
            print("No parameter was updated")

    # =========================================================
    def calc_SH_data(self, data=None, ref=True):
        xx, yy   = self.SH_xx , self.SH_yy        
        ncx, ncy = xx.size - 1, yy.size - 1

        if data is None:
            self.live_img = self.shm_im.get_data(check=self.im_cnt,reform=True)
        else:
            self.live_img = data

        self.im_cnt = self.shm_im.get_counter()  # image counter
        bckgd = self.live_img.mean()

        #self.live_img[self.live_img <= self.vmin] = self.vmin

        for j in xrange(ncy):
            y0, y1 = int(np.round(yy[j])), int(np.round(yy[j+1]))

            for i in xrange(ncx):
                x0, x1 = int(np.round(xx[i])), int(np.round(xx[i+1]))

                sub_arr           = self.live_img[y0:y1,x0:x1]
                self.SH_phot[j,i] = sub_arr.max()# - bckgd# max()

                sub_arr[sub_arr < self.threshold] = self.threshold
                #if self.SH_phot[j,i] > self.threshold: # original line
                if self.SH_phot[j,i] > 1.3 * bckgd: # 30 % above background?
                    (yc, xc) = centroid_position_0(sub_arr)
                else:
                    (yc, xc) = (self.SH_xref[j,i], self.SH_yref[j,i])

                self.SH_xtmp[j,i] = xc
                self.SH_ytmp[j,i] = yc

        if ref is True:
            self.SH_xref = self.SH_xtmp.copy()
            self.SH_yref = self.SH_ytmp.copy()

        self.SH_xslp = self.SH_xtmp - self.SH_xref
        self.SH_yslp = self.SH_ytmp - self.SH_yref

        self.SH_phot -= self.threshold
        self.SH_xslp[self.SH_phot < 0] = 0.0
        self.SH_yslp[self.SH_phot < 0] = 0.0

        self.SH_comb[0] = self.SH_xslp
        self.SH_comb[1] = self.SH_yslp
        self.SH_comb[2] = self.SH_phot / self.SH_phot.max()

        self.shm_comb.set_data(self.SH_comb)
        #self.shm_xslp.set_data(self.SH_xslp)
        #self.shm_yslp.set_data(self.SH_yslp)
        self.shm_phot_inst.set_data(self.SH_phot)
        
        # here is information about the tip-tilt in pixels!
        # weighted mean version!
        self.ttx_mean = np.average(self.SH_xslp, weights=self.SH_phot)
        self.tty_mean = np.average(self.SH_yslp, weights=self.SH_phot)
        #self.ttx_mean = np.sum(self.SH_xslp * self.SH_phot) / np.sum(self.SH_phot)
        #self.tty_mean = np.sum(self.SH_yslp * self.SH_phot) / np.sum(self.SH_phot)

        # original version below
        #self.ttx_mean = np.median(self.SH_xslp[self.SH_phot > self.threshold])
        #self.tty_mean = np.median(self.SH_yslp[self.SH_phot > self.threshold])

        self.update_log()
    # =========================================================
    def update_log(self,):
        self.ttx_log.append(self.ttx_mean)
        self.tty_log.append(self.tty_mean)
        if len(self.ttx_log) > self.log_len:
            self.ttx_log.pop(0)
            self.tty_log.pop(0)

    # =========================================================
    def loop(self,):
        self.keepgoing = True
        while self.keepgoing:
            self.calc_SH_data(ref=False)
        print("WFS monitoring is now off")

    # =========================================================
    def start(self,):
        print("start %s" % ("YEAH!", ))
        print(self.SH_xx)
        print(self.SH_yy)
        self.calc_SH_data(ref=True)
        
# ==========================================================
# ==========================================================
if __name__ == "__main__":
    mon = WFS()
    mon.start()
    
    
