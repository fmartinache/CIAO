#!/usr/bin/env python3

import sys
import numpy as np
from xaosim.shmlib import shm
import datetime
import time
#import fitsio
#from fitsio import FITS,FITSHDR
import astropy.io.fits as pf
import pdb
import functools

from scipy.ndimage import center_of_mass as comass

import multiprocessing as mp

def centroid_position_slice(aslice, img):
    ''' --------------------------------------------------------------
    Wrapper used for multiprocessing use case
    -------------------------------------------------------------- '''
    return centroid_position(img[aslice])

def maxof_slice(aslice, img):
    ''' --------------------------------------------------------------
    Wrapper used for multiprocessing use case
    -------------------------------------------------------------- '''
    return img[aslice].max()

def centroid_position(img):
    ''' --------------------------------------------------------------
    Wrapper function to make it easy to change different algorithms
    -------------------------------------------------------------- '''
    return comass(img) # the scipy center_of_mass?
    #return maximum_position(img) # another possibility

def maximum_position(img):
    ''' --------------------------------------------------------------
    Returns the (x,y) position of the centroid of the array *img*
    -------------------------------------------------------------- '''
    return np.unravel_index(np.argmax(img, axis=None), img.shape)

# =============================================================================
# =============================================================================

class WFS():
    # =========================================================
    def __init__(self, shmf="/dev/shm/SHcam.im.shm", ncpu=False):
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
        self.shm_im = shm(shmf)                 # shared mem. image
        self.isz = self.shm_im.mtdata['x']      # should be 128
        self.im_cnt = self.shm_im.get_counter() # image counter
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
        nx = int(np.round((self.isz - self.SH_x0) / self.SH_dx) + 2)
        xx = np.round(self.SH_x0+np.arange(nx-1)*self.SH_dx).astype(int)
        if (self.isz - xx[-1] > self.SH_dx/2+2):
            xx = np.append(xx, self.isz)
        self.SH_xx = xx
        
        ny = int(np.round((self.isz - self.SH_y0) / self.SH_dy) + 2)
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
            '/dev/shm/SNR.im.shm', data=self.SH_phot, verbose=False) # experiment!
        self.shm_phot_inst = shm(
            '/dev/shm/phot_inst.im.shm', data=self.SH_phot, verbose=False)

        self.shm_comb = shm('/dev/shm/comb.im.shm',data=self.SH_comb,verbose=True)
        
        for jj in range(ncy):
            y0, y1 = yy[jj], yy[jj+1]
            for ii in range(ncx):
                x0, x1 = xx[ii], xx[ii+1]
                self.SH_xref[jj,ii] = np.round(0.5 * (x1 - x0))
                self.SH_yref[jj,ii] = np.round(0.5 * (y1 - y0))

        pf.writeto("SHrefx.fits", self.SH_xref, overwrite=True)
        pf.writeto("SHrefy.fits", self.SH_yref, overwrite=True)
        

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
    def centroid_position_cell0(self, aslice):
        ''' -------------------------------------------------------
        Class specific wrapper function that returns the centroid
        for a given slice of the live image

        Original form - using the self.live_img data
        ------------------------------------------------------- '''
        return centroid_position(self.live_img[aslice])

    # =========================================================
    def centroid_position_cell(self, minimg):
        ''' -------------------------------------------------------
        Class specific wrapper function that returns the centroid
        for a given slice of the live image
        ------------------------------------------------------- '''
        return centroid_position(img[aslice])

    def cell_img(self, aslice, img):
        ''' Wrapper for a test '''
        return img[aslice]
    
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
        self.live_img[self.live_img < self.threshold] = 0.0

        # =======================================================
        if self.mproc is True:
            # -------------------------------
            # multi-processor scenario (TBC)
            # -------------------------------
            #pool = mp.Pool(self.ncpu) # split computation on CPU pool
            cen_pos = functools.partial(centroid_position_slice, img=self.live_img)
            cen_list = np.array(list(self.pool.map(cen_pos, self.SH_cells)))

            cell_max = functools.partial(maxof_slice, img=self.live_img)
            max_list = np.array(list(self.pool.map(cell_max, self.SH_cells)))
            #max_list = np.array(list(map(self.max_cell, self.SH_cells)))
            
            #pool.close()
        else:
            # -------------------------------
            # this is the single CPU scenario
            # -------------------------------
            # TEST
            # subimgs = list(map(lambda aslice: self.live_img[aslice], self.SH_cells))
            # cen_list = np.array(list(map(lambda x: centroid_position(x), subimgs)))

            cen_list = np.array(list(map(lambda x: centroid_position(self.live_img[x]), self.SH_cells)))

            # # get the list of centroid positions
            # cen_list = np.array(list(map(self.centroid_position_cell, self.SH_cells)))
            # # and get the list of maximum values
            max_list = np.array(list(map(self.max_cell, self.SH_cells)))

        # --------------------------------------------------------
        # IMPORTANT HERE: x and y slope measurements can be swapped
        # here to match actual DM orientation !!
        # --------------------------------------------------------            
        self.SH_phot = max_list.reshape((ncy,ncx)).astype('float')
        self.SH_ytmp = cen_list[:,0].reshape((ncy, ncx))
        self.SH_xtmp = cen_list[:,1].reshape((ncy, ncx))# * 0.0 # DEBUG TOOL
        
        # populate the different data structures with the new info
        # =========================================================
        if ref is True:
            print("NEW REFERENCE!")
            self.SH_xref = self.SH_xtmp.copy()
            self.SH_yref = self.SH_ytmp.copy()
            pf.writeto("SHrefx.fits", self.SH_xref, overwrite=True)
            pf.writeto("SHrefy.fits", self.SH_yref, overwrite=True)

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
    
    
