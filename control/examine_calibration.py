#!/usr/bin/env python

import numpy as np
from numpy.linalg import svd
import astropy.io.fits as pf
import os
from ciao_wfc import DM

home = os.getenv('HOME')
conf_dir = home+'/.config/ciao/'

MM = pf.getdata(conf_dir+"ZER_CAL.fits")

U,S,Vh = svd(MM)

plt.figure()
plt.plot(S)

mydm = DM()

plt.figure()
plt.imshow(mydm.list_2_map(U[0,:]), interpolation="nearest")
