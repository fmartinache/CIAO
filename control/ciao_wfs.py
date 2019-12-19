#!/usr/bin/env python3

import sys
import numpy as np
from xaosim.shmlib import shm
import datetime
import time
import fitsio
from fitsio import FITS,FITSHDR
import astropy.io.fits as pf

from scipy.ndimage import center_of_mass as comass
# would be worth giving this function a shot!

import multiprocessing as mp

# !!! multi-CPU not working yet !!!

def ext_centroid_position_cell(args, **kwargs):
    return WFS.centroid_position_cell(*args, **kwargs)

def ext_max_cell(args, **kwargs):
    return WFS.max_cell(*args, **kwargs)

def centroid_position(img):
    ''' --------------------------------------------------------------
    Wrapper function to make it easy to change different algorithms
    -------------------------------------------------------------- '''
    #return comass(img) # the scipy center_of_mass?
    return maximum_position(img) # another possibility
    #return centroid_position_0(img)

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
    return (yc, xc)

def maximum_position(img):
    ''' --------------------------------------------------------------
    Returns the (x,y) position of the centroid of the array *img*
    -------------------------------------------------------------- '''
    return np.unravel_index(np.argmax(img, axis=None), img.shape)

# =============================================================================
# =============================================================================

class WFS():
    # =========================================================
    def __init__(self, shmf="/tmp/SHcam.im.shm", ncpu=False):
        ''' -------------------------------------------------------
        Default Wavefront Sensor (WFS) class constructor.

        Parameters:
        ----------
        - shmf: shared memory data structure where to read images
        - ncpu: (integer) number of cpus to use for computation
        ------------------------------------------------------- '''
        # -------------------------
        # to be provided externally
        # -------------------------
        self.SH_x0, self.SH_y0 = 0, 0   # properties of the grid
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
        self.keepgoing = False # close-loop flag

        self.mproc = False # multi-processor option
        if ncpu is not False:
            ncpumax = mp.cpu_count()
            print("Multi processor option selected")
            print("%d CPUs available: %d requested" % (ncpumax, ncpu))
            self.mproc = True
            self.ncpu = ncpu
            self.pool = mp.Pool(ncpu) # split computation on CPU pool

    
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

        # -------------------------------------------------------------
        # with these coordinates computed, we can now define the slices
        # of the image that will correspond to the different cells
        # Hopefully, these slices will accelerate the computation !!
        # -------------------------------------------------------------
        self.SH_cells = []
        for jj in range(ny-2):
            for ii in range(nx-2):
                self.SH_cells.append([slice(yy[jj], yy[jj+1]), slice(xx[ii], xx[ii+1])])

    # =========================================================
    def define_SH_data(self):
        ''' -------------------------------------------------------
        Creates empty arrays that will be used to store the information
        extracted from the analysis of the SH images: photometry, slope
        in the x and y directions, position of the reference

        Creates the shared memory data structures that will be used
        to communicate this information with other processes:

        - shm_phot_inst: current illumination of the apertures
        - shm_comb     : xy slopes & photometry combo

        Currently not used: shm_SNR
        ------------------------------------------------------- '''
        xx, yy   = self.SH_xx , self.SH_yy  # position of the SH cell edges
        ncx, ncy = xx.size - 1, yy.size - 1 # number of cells
        # -------------
        # empty arrays
        # -------------
        self.SH_phot = np.zeros((ncy, ncx)) # photometry
        self.SH_xslp = np.zeros((ncy, ncx)) # x-slope
        self.SH_yslp = np.zeros((ncy, ncx)) # y-slope
        self.SH_xref = np.zeros((ncy, ncx)) # x-reference
        self.SH_yref = np.zeros((ncy, ncx)) # y-reference
        self.SH_comb = np.zeros((3, ncy, ncx)) # xy slopes + phot combined
        # ----------------------------
        # additional utility arrays
        # ----------------------------        
        self.SH_xtmp = np.zeros((ncy, ncx))
        self.SH_ytmp = np.zeros((ncy, ncx))

        self.shm_SNR = shm(
            '/tmp/SNR.im.shm', data=self.SH_phot, verbose=False) # experiment!
        self.shm_phot_inst = shm(
            '/tmp/phot_inst.im.shm', data=self.SH_phot, verbose=False)

        self.shm_comb = shm('/tmp/comb.im.shm',data=self.SH_comb,verbose=False)
        
        for jj in xrange(ncy):
            #y0, y1 = int(np.round(yy[jj])), int(np.round(yy[jj+1]))
            y0, y1 = yy[jj], yy[jj+1]
            for ii in xrange(ncx):
                #x0, x1 = int(np.round(xx[ii])), int(np.round(xx[ii+1]))
                x0, x1 = xx[ii], xx[ii+1]
                self.SH_xref[jj,ii] = np.round(0.5 * (x1 - x0)) + 1
                self.SH_yref[jj,ii] = np.round(0.5 * (y1 - y0)) + 1

        pf.writeto("SHrefx.fits", self.SH_xref, clobber=True)
        pf.writeto("SHrefy.fits", self.SH_yref, clobber=True)
        

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
    def centroid_position_cell(self, aslice):
        ''' -------------------------------------------------------
        Class specific wrapper function that returns the centroid
        for a given slice of the live image
        ------------------------------------------------------- '''
        return centroid_position(self.live_img[aslice])

    # =========================================================
    def max_cell(self, aslice):
        ''' -------------------------------------------------------
        Another such wrapper that should eventually be merged with
        the centroid since both info can easily be aggregated in 
        one go.
        ------------------------------------------------------- '''
        return self.live_img[aslice].max()
    
    # =========================================================
    def calc_SH_data(self, data=None, ref=True):
        ''' -------------------------------------------------------
        Calculates the Shack-Hartman data for the current image.

        Parameters:
        - ref: if True -> current image taken as new reference (bool)
        - data: if not None, info computed from the provided image
        ------------------------------------------------------- '''
        ncx, ncy = self.SH_xx.size - 1, self.SH_yy.size - 1

        if data is None:
            self.live_img = self.shm_im.get_data(check=self.im_cnt,reform=True)
        else:
            self.live_img = data
                    
        self.im_cnt = self.shm_im.get_counter()  # update image counter
        self.live_img[self.live_img < self.threshold] = 0.0#self.threshold
        #self.live_img -= self.threshold
        #self.live_img[self.live_img <=0] = 0.0

        # =======================================================
        # trying to write this the right way now...
        # optional self.pool.map(...)
        # despite my efforts, the pool trick doesn't seem to work with python2.7

        if self.mproc is True:
            # -------------------------------
            # multi-processor scenario (TBC)
            # -------------------------------
            #pool = mp.Pool(self.ncpu) # split computation on CPU pool
            cen_list = np.array(list(pool.map(ext_centroid_position_cell, zip([self]*len(self.SH_cells), self.SH_cells))))
            max_list = np.array(list(pool.map(ext_max_cell, zip([self]*len(self.SH_cells), self.SH_cells))))
            #pool.close()
        else:
            # -------------------------------
            # this is the single CPU scenario
            # -------------------------------

            # get the list of centroid positions
            cen_list = np.array(list(map(self.centroid_position_cell, self.SH_cells)))
            # and get the list of maximum values
            max_list = np.array(list(map(self.max_cell, self.SH_cells)))

        self.SH_phot = max_list.reshape((ncy,ncx)).astype('float')
        self.SH_ytmp = cen_list[:,0].reshape((ncy, ncx))
        self.SH_xtmp = cen_list[:,1].reshape((ncy, ncx))# * 0.0
        
        # populate the different data structures with the new info
        # =========================================================
        if ref is True:
            print("NEW REFERENCE!")
            self.SH_xref = self.SH_xtmp.copy()
            self.SH_yref = self.SH_ytmp.copy()
            pf.writeto("SHrefx.fits", self.SH_xref, clobber=True)
            pf.writeto("SHrefy.fits", self.SH_yref, clobber=True)

        self.SH_xslp = self.SH_xtmp - self.SH_xref
        self.SH_yslp = (self.SH_ytmp - self.SH_yref)

        #self.SH_phot -= self.threshold
        discard_them = self.SH_phot <= 1500#self.threshold * 2
        self.SH_phot[discard_them] = 0.0
        self.SH_xslp[discard_them] = 0.0
        self.SH_yslp[discard_them] = 0.0

        self.SH_comb[0] = self.SH_xslp
        self.SH_comb[1] = self.SH_yslp
        self.SH_comb[2] = self.SH_phot / self.SH_phot.max()

        self.shm_comb.set_data(self.SH_comb)
        self.shm_phot_inst.set_data(self.SH_phot)
        
        # here is information about the tip-tilt in pixels!
        # weighted mean version!
        self.ttx_mean = np.average(self.SH_xslp, weights=self.SH_phot)
        self.tty_mean = np.average(self.SH_yslp, weights=self.SH_phot)

    # =========================================================
    def calc_SH_data_old(self, data=None, ref=True):
        ''' -------------------------------------------------------
        Calculates the Shack-Hartman data for the current image.

        Parameters:
        - ref: if True -> current image taken as new reference (bool)
        - data: if not None, info computed from the provided image

        This version of the function was idenfied as limiting the
        frequency of the loop. It is for now kept as a reference.
        ------------------------------------------------------- '''
        xx, yy   = self.SH_xx , self.SH_yy        
        ncx, ncy = xx.size - 1, yy.size - 1

        if data is None:
            self.live_img = self.shm_im.get_data(check=self.im_cnt,reform=True)
        else:
            self.live_img = data
            
        #~ fitsio.write("./xx.fits",xx)
        #~ fitsio.write("./yy.fits",yy)
        #~ fitsio.write("./live_img.fits",self.live_img)
        
        self.im_cnt = self.shm_im.get_counter()  # image counter
        bckgd = self.live_img.mean()

        #self.live_img[self.live_img <= self.vmin] = self.vmin
        #~ mycount = 0
        #~ time_start = time.time()
        #~ time_diff = 0.
        for j in xrange(ncy):
            y0, y1 = int(np.round(yy[j])), int(np.round(yy[j+1]))

            for i in xrange(ncx):
                x0, x1 = int(np.round(xx[i])), int(np.round(xx[i+1]))

                sub_arr           = self.live_img[y0:y1,x0:x1]
                self.SH_phot[j,i] = sub_arr.max()# - bckgd# max()

                sub_arr[sub_arr < self.threshold] = self.threshold
                #if self.SH_phot[j,i] > self.threshold: # original line
                if self.SH_phot[j,i] > 1.3 * bckgd: # 30 % above background?
                    #~ time_tmp = time.time()
                    (yc, xc) = centroid_position(sub_arr)
                    #~ time_diff += time.time() - time_tmp
                    #~ mycount += 1
                else:
                    (yc, xc) = (self.SH_xref[j,i], self.SH_yref[j,i])

                self.SH_xtmp[j,i] = xc
                self.SH_ytmp[j,i] = yc
		
        #~ print("calc_SH_data timing: %.3e %.3e %d" % (time.time() - time_start, time_diff/mycount,mycount))

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
        self.shm_phot_inst.set_data(self.SH_phot)
        
        # here is information about the tip-tilt in pixels!
        # weighted mean version!
        self.ttx_mean = np.average(self.SH_xslp, weights=self.SH_phot)
        self.tty_mean = np.average(self.SH_yslp, weights=self.SH_phot)
	
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
    
    # =========================================================
    def stop(self,):
        # release the CPUs used by multi-processing
        if self.mproc is True:
            print("OK")
            #print("relieving the %d CPUs from duty" % (self.ncpu))
            #self.pool.close()
# ==========================================================
# ==========================================================
if __name__ == "__main__":
    mon = WFS(shmf="/tmp/ixon.im.shm")
    mon.start()
    
    
