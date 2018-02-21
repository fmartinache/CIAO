#!/usr/bin/env python

import pygame, sys, os
from pygame.locals import *
import numpy as np
import struct
import pdb

simulation = False# True # set to False when running the ALPAO DM

home = os.getenv('HOME')

sys.path.append(home+"/src/python/libs/")
from scexao_shm import shm

import PYG_GUI
from PYG_GUI import *

adelay = 50 # action delay (keys and other actions)
pygame.init()
pygame.key.set_repeat(adelay, 20)

FPS = 20
fpsClock = pygame.time.Clock()
XW, YW = 600, 800

screen = pygame.display.set_mode((XW, YW), 0, 32)
pygame.display.set_caption("POKE' EM")

# =============================================================
dms = 11
xx, yy = np.meshgrid(np.arange(dms)-dms/2, np.arange(dms)-dms/2)
xb = XW/2 + 50 * (np.arange(dms) - dms/2)
yb = YW/2 + 50 * (dms/2 - np.arange(dms))

actuators = []
vals = []

mask = np.ones((dms, dms),    dtype=np.float32) # to remove non-existing actuators
mask[:3,0] = np.nan
mask[0,:3] = np.nan
mask[1, 1] = np.nan
mask[-3:,0] = np.nan
mask[-1,:3] = np.nan
mask[-2,1]  = np.nan
mask[0,-3:] = np.nan
mask[1,-2:] = np.nan
mask[2,-1]  = np.nan
mask[-1,-3:] = np.nan
mask[-2,-2:] = np.nan
mask[-3,-1]  = np.nan

flat_mask = mask.flatten()

chn = 5                         # DM channel for this program
disp = shm('/tmp/dmdisp%d.im.shm' % (chn,), verbose=False)

# =============================================================

index = 0
for i in xrange(dms):
    for j in xrange(dms):
        if mask[i,j] > 0.0:
            actuators.append(pyg_button(screen, (xb[i], yb[j], 45, 45), "%d" % (index,)))
            vals.append(pyg_label(screen, (xb[i], yb[j]+15), "%+.2f" % (0.0), 20))
            index += 1

nbAct = actuators.__len__()
values = [0.] * nbAct # values sent to the DM
select = 48 # central actuator selected at startup
step = 0.001
amax, amin = 0.5, -0.5 # max commands for the DM

values = [0.] * nbAct
    
for i in xrange(nbAct):
    vals[i].text = "%+.3f" % (values[i],)

# ===========================================================
def list2map2D(data):
    ''' ------------------------------------------------
    Convert a list of 97 voltages back into a 2D map for
    easy 2D representation.
    ------------------------------------------------ '''
    global flat_mask, nbAct
    res = np.zeros((dms,  dms), dtype=np.float32)

    i0 = 0
    for k in xrange(dms*dms):
        if not np.isnan(flat_mask[k]):
            i = k / dms
            j = k % dms
            res[i,j] = data[i0]
            i0 += 1
    return(res)

# ===========================================================
def dm_send(values):
    for i in xrange(nbAct):
        vals[i].text = "%+.3f" % (values[i],)

    disp.set_data(list2map2D(values))


# =============================================================
dm_send(values) # send last saved map to DM at startup
# =============================================================

keepgoing = True
while keepgoing:

    for event in pygame.event.get():

        if event.type == QUIT:
            keepgoing = False
            break

        # ---------------------------------------------
        #               KEYBOARD actions
        # ---------------------------------------------
        elif event.type == KEYDOWN:
            mmods = pygame.key.get_mods()

            # ----- up to move one actuator up -----
            if event.key == K_UP:
                select += 1
                if (select >= nbAct):
                    select = 0

            # ----- down to do the opposite -----
            if event.key == K_DOWN:
                select -= 1
                if (select < 0):
                    select = nbAct-1

            # ----- left/right to move actuators ------
            if event.key == K_RIGHT:
                mult = 1
                if (mmods & KMOD_LSHIFT):
                    mult = 10.0
                values[select] += mult * step
                if values[select] > amax:
                    values[select] = amax
                dm_send(values)

            if event.key == K_LEFT:
                mult = 1
                if (mmods & KMOD_LSHIFT):
                    mult = 10.0
                values[select] -= mult * step
                if values[select] < amin:
                    values[select] = amin
                dm_send(values)

            if event.key == K_r:
                if (mmods & KMOD_LSHIFT):
                    values = [0.] * nbAct
                    dm_send(values)

            # ----- ESCAPE used to quit the GUI -----
            if event.key == K_ESCAPE:
                keepgoing = False
                break
        # ---------------------------------------------
        #               MOUSE actions
        # ---------------------------------------------
        elif event.type == MOUSEMOTION:
            mx, my = event.pos

        elif event.type == MOUSEBUTTONDOWN:
            mx, my = event.pos
            clicked = True
            mmods = pygame.key.get_mods()

            if event.button == 4:
                mult = 1
                if (mmods & KMOD_LSHIFT):
                    mult = 10.0
                values[select] += mult * step
                if values[select] > amax:
                    values[select] = amax
                dm_send(values)
                
            if event.button == 5:
                mult = 1
                if (mmods & KMOD_LSHIFT):
                    mult = 10.0
                values[select] -= mult * step
                if values[select] < amin:
                    values[select] = amin
                dm_send(values)

            # --------- CLICK used to switch modes ------
            for i in xrange(nbAct):
                mytest = actuators[i].thisRect.collidepoint(mx, my)
                if (mytest):
                    select = i


        else:
            break

    for i in xrange(nbAct):
        active = False
        if i == select:
            active = True

        actuators[i].updt(active)
        vals[i].updt(active)
    fpsClock.tick(FPS)
    pygame.display.flip()


#np.savetxt('flat_map.txt', values)

pygame.quit()
sys.exit()
# =============================================================
