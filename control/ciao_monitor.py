#!/usr/bin/env python

import sys
import numpy as np
from xaosim.shmlib import shm

def centroid_position(img):
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
        xc = np.sum(mprofx * xxc) / np.sum(mprofx)
    else:
        xc = 0.0
    if denomy > 1:
        yc = np.sum(mprofy * yyc) / np.sum(mprofy)
    else:
        yc = 0.0
    return((yc, xc))

class monitor():
    # =========================================================
    def __init__(self, shmf="/tmp/SHcam.im.shm"):

        # -------------------------
        # to be provided externally
        # -------------------------
        self.SH_x0, self.SH_y0 = 0, 0   # properties of the grid
        self.SH_dx, self.SH_dy = 13, 13 # properties of the grid
        self.threshold = 10.0
        
        # -----------------------------
        # shared memory data structures
        # -----------------------------
        self.im_shm = shm(shmf)                  # shared mem. image
        self.isz = self.im_shm.mtdata['size'][0] # should be 128

        self.update_grid()
        self.define_SH_data()
        
    # =========================================================
    def update_grid(self,):
        ''' -------------------------------------------------------
        Updates the SH cells positions for new SH grid parameters
        ------------------------------------------------------- '''
        nx = int((self.isz - self.SH_x0) / self.SH_dx + 2)
        xx = self.SH_x0 + np.arange(nx-1) * self.SH_dx
        if (self.isz - xx[-1] > self.SH_dx/2+2):
            xx = np.append(xx, self.isz)
        self.SH_xx = xx
        
        ny = int((self.isz - self.SH_y0) / self.SH_dy + 2)
        yy = self.SH_y0 + np.arange(ny-1) * self.SH_dy
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
        # additional arrays for convenience
        self.SH_xtmp = np.zeros((ncy, ncx))
        self.SH_ytmp = np.zeros((ncy, ncx))
        self.calc_SH_data(ref=True)

        self.shm_phot_gain = shm('/tmp/phot_gain.im.shm', data=self.SH_phot, 
                                 verbose=False) # experiment!

        self.shm_phot_inst = shm('/tmp/phot_inst.im.shm', data=self.SH_phot, 
                                 verbose=False) # experiment!

        self.shm_xslp = shm('/tmp/xslp.im.shm', data=self.SH_xslp, verbose=False)
        self.shm_yslp = shm('/tmp/yslp.im.shm', data=self.SH_yslp, verbose=False)
        
    # =========================================================
    def calc_SH_data(self, data=None, ref=True):
        xx, yy   = self.SH_xx , self.SH_yy        
        ncx, ncy = xx.size - 1, yy.size - 1

        if data is None:
            self.live_img = self.im_shm.get_data(True, True)
        else:
            self.live_img = data
            
        #self.live_img[self.live_img <= self.vmin] = self.vmin

        for j in xrange(ncy):
            y0, y1 = int(np.round(yy[j])), int(np.round(yy[j+1]))

            for i in xrange(ncx):
                x0, x1 = int(np.round(xx[i])), int(np.round(xx[i+1]))

                self.SH_xref[j,i] = 0.5 * (x1 - x0) #self.SH_dx
                self.SH_yref[j,i] = 0.5 * (y1 - y0) #self.SH_dy

                sub_arr           = self.live_img[y0:y1,x0:x1]
                self.SH_phot[j,i] = sub_arr.max() #- self.threshold

                sub_arr[sub_arr < self.threshold] = self.threshold
                if self.SH_phot[j,i] > self.threshold:
                    (yc, xc) = centroid_position(sub_arr)
                else:
                    (yc, xc) = (self.SH_xref[j,i], self.SH_yref[j,i])

                self.SH_xtmp[j,i] = xc
                self.SH_ytmp[j,i] = yc

        if ref is True:
            self.SH_xref = self.SH_xtmp.copy()
            self.SH_yref = self.SH_ytmp.copy()

        self.SH_xslp = self.SH_xtmp - self.SH_xref
        self.SH_yslp = self.SH_ytmp - self.SH_yref

        self.SH_xslp[self.SH_phot <= self.threshold] = 0.0
        self.SH_yslp[self.SH_phot <= self.threshold] = 0.0

        # here is information about the tip-tilt in pixels!
        # weighted mean version!
        #self.ttx_mean = np.sum(self.SH_xslp * self.SH_phot) / np.sum(self.SH_phot)
        #self.tty_mean = np.sum(self.SH_yslp * self.SH_phot) / np.sum(self.SH_phot)

        # original version below
        self.ttx_mean = np.median(self.SH_xslp[self.SH_phot > self.threshold])
        self.tty_mean = np.median(self.SH_yslp[self.SH_phot > self.threshold])

    # =========================================================
    def start(self,):
        print("start %s" % ("YEAH!", ))
        print(self.SH_xx)
        print(self.SH_yy)

        
# ==========================================================
# ==========================================================
if __name__ == "__main__":
    mon = monitor()
    mon.start()
    
