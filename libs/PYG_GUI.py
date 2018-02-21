#!/usr/bin/env python

import pygame, sys
from pygame.locals import *
import time

import numpy as np
import matplotlib.cm as cm
import serial
import datetime
import pdb


WHITE = (255, 255, 255)
GREEN = (  0, 255,   0) 
DGREN = (  0, 102,  51)
TEAL  = (  0, 102, 102)
BLUE  = (  0,   0, 255)
RED   = (255,   0,   0)
BLK   = (  0,   0,   0)
DGREY = ( 50,  50,  50)
CYAN  = (  0, 255, 255)
GOLD  = (255, 215,   0)
ORANG = (255, 128,   0)
BLACK = (  0,   0,   0)

FGCOL = WHITE  # foreground color (text)
BGCOL = BLK    # background color
BTCOL = BLUE   # *button* color

# ==================================================
#  creates a convenient HELP message, in the window
# ==================================================
def help_menu(container, title, message, logo_file=False):
    FPS = 60
    end_it = False

    hlp_msg = message
    font1 = pygame.font.SysFont("default",   40)
    hfont = pygame.font.SysFont("monospace", 20)
    hrect = container.get_rect()
    XW = hrect[2]
    YW = hrect[3]
    HZ = 20

    if logo_file:
        speed    = [8, 1]
        try:
            logo     = pygame.image.load(logo_file)
            logo_rct = logo.get_rect()
        except:
            logo_file = False


    while (not end_it):
        container.fill(BGCOL)

        msg1 = font1.render(title, True, FGCOL)

        lines = hlp_msg.split('\n')
        nlines = lines.__len__()

        for event in pygame.event.get():
            if event.type == KEYDOWN:
                end_it = True

        if logo_file:
            logo_rct = logo_rct.move(speed)
            if logo_rct.left < 0 or logo_rct.right > XW:
                speed[0] = -speed[0]
            if logo_rct.top < 0 or logo_rct.bottom > YW:
                speed[1] = -speed[1]
            container.blit(logo, logo_rct)

        container.blit(msg1, msg1.get_rect(center=(XW/2, HZ)))

        for i in xrange(nlines):
            y0 = (i+4) * HZ
            msgx = hfont.render(lines[i], True, FGCOL)
            textpos = msgx.get_rect()
            textpos.left = 40
            textpos.centery = y0
            container.blit(msgx, textpos)
        pygame.display.flip()
        pygame.time.Clock().tick(FPS)

# ==================================================
class pyg_button:
    ''' ------------------------------------------------
        class to draw a "button" in a PYG GUI
    ------------------------------------------------ '''
    def __init__(self, container, (xc, yc, ww, hh), label):
        # storing information in the object structure
        self.xc, self.yc = xc, yc
        self.ww, self.hh = ww, hh
        self.x0, self.y0 = xc - ww/2, yc - hh/2
        self.label = label
        self.thisRect = pygame.Rect(self.x0, self.y0, self.ww, self.hh)

        self.myfont = pygame.font.SysFont("default",   28)
        self.BTCOL  = BLUE       # cursor color (can be updated)
        self.MYBGCOL  = self.BTCOL # internal cursor color (cannot be changed)
        self.MYFGCOL  = RED
        self.container = container

    def updt(self, active=False, pressed=False):
        if active:
            self.MYBGCOL = GOLD
            self.MYFGCOL = RED
        else:
            self.MYBGCOL = self.BTCOL
            self.MYFGCOL = FGCOL
        if pressed:
            self.MYBGCOL = RED
            self.MYFGCOL = WHITE

        pygame.draw.rect(self.container, self.MYBGCOL, self.thisRect, 0)
        pygame.draw.rect(self.container, self.MYFGCOL, self.thisRect, 1)

        lbl = self.myfont.render(self.label,  True, self.MYFGCOL)
        self.container.blit(lbl, lbl.get_rect(center=(self.xc, self.yc)))

# ==================================================
class pyg_label:
    def __init__(self, container, (xc, yc), text, fontsize=False):
        self.xc, self.yc = xc, yc
        if not fontsize:
            fontsize = 28
        self.myfont = pygame.font.SysFont("default", fontsize)
        self.FGCOLOR = WHITE
        self.BGCOLOR = DGREY
        self.container = container
        self.text = text

    def updt(self, active=False):
        if active:
            self.FGCOLOR = RED
        else:
            self.FGCOLOR = WHITE

        lbl = self.myfont.render(self.text, True, self.FGCOLOR)
        self.container.blit(lbl, lbl.get_rect(center=(self.xc, self.yc)))

# ==================================================
class pyg_cursor:
    ''' ------------------------------------------------
    class to draw a horizontal cursor in a PYG GUI
    ------------------------------------------------ '''
    def __init__(self, container, (xc, yc, hh, ww), active=False):
        # =====================
        self.vmin   = 0.0 # default behaviour for the cursor
        self.vmax   = 1.0 # need to be manually updated
        self.vcur   = 0.0 
        self.step   = 0.01 # normal increment

        # =====================
        self.xc, self.yc = xc, yc # center of the cursor on screen
        self.hh, self.ww = hh, ww # dimensions

        self.x0, self.y0 = xc - ww/2, yc - hh/2
        self.thisRect = pygame.Rect(self.x0, self.y0, self.ww, self.hh)

        self.x1 = self.x0 + ww * \
                  (self.vcur-self.vmin) / (self.vmax-self.vmin) - ww/20
        self.y1 = yc-self.hh

        self.myfont = pygame.font.SysFont("default",   28)
        self.mylbl = "X"

        self.BTCOL  = BLUE       # cursor color (can be updated)
        self.MYBGCOL  = self.BTCOL # internal cursor color (cannot be changed)
        self.MYFGCOL  = RED
        self.container = container

    # ===========================================
    def increment(self, step=False):
        temp = self.vcur
        if step:
            temp += step
        else:
            temp += self.step

        if (temp <= self.vmax):
            self.vcur = temp
        self.updt(active=True)

    # ===========================================
    def decrement(self, step=False):
        temp = self.vcur
        if step:
            temp -= step
        else:
            temp -= self.step

        if (temp >= self.vmin):
            self.vcur = temp
        self.updt(active=True)

    # ===========================================
    def updt(self, vcur=False, active=False, 
             vmin=False, vmax=False, label=False, step=False):
        if vcur:
            self.vcur = vcur
        if (vmin):
            self.vmin = vmin
        if (vmax):
            self.vmax = vmax
        if (step):
            self.step = step

        self.x1 = self.x0 + self.ww * \
                  (self.vcur-self.vmin) / (self.vmax-self.vmin) - self.ww/20

        if active:
            self.MYBGCOL = GOLD
            self.MYFGCOL = RED
        else:
            self.MYBGCOL = self.BTCOL
            self.MYFGCOL = FGCOL

        #pdb.set_trace()

        pygame.draw.rect(self.container, self.MYBGCOL, 
                         (self.x0, self.y0,self.ww, self.hh), 0)
        pygame.draw.rect(self.container, self.MYFGCOL, 
                         (self.x0, self.y0,self.ww, self.hh), 1)

        pygame.draw.rect(self.container, self.MYBGCOL, 
                         (self.x1, self.y1,self.hh, 2*self.hh), 0)
        pygame.draw.rect(self.container, self.MYFGCOL, 
                         (self.x1, self.y1,self.hh, 2*self.hh), 1)

        if label:
            self.mylbl = label
        else:
            if self.step > 0.9:
                self.mylbl = "%d" % (self.vcur)
            else:
                self.mylbl = "%+.2f" % (self.vcur)
        hdl = self.myfont.render(self.mylbl,  True, self.MYFGCOL)

        zn = hdl.get_rect()
        zn.center = (self.xc, self.yc)
        self.container.blit(hdl, zn)

# ==================================================
class pygimg:
    ''' ------------------------------------------------
    2D image plot in a pygame GUI container (typically screen)
    ------------------------------------------------ '''
    # ================
    def __init__(self, container, XW, YW, XC, YC):
        self.XW      = XW
        self.YW      = YW
        self.LW      = 1     # default line width
        self.BGCOLOR = BLUE  # default color theme
        self.FGCOLOR = WHITE # default color theme
        self.PLT_CLR = BLUE  # default color theme
        self.cmap    = cm.gray
        # ---
        self.container = container
        self.zone = pygame.Surface((XW, YW))
        self.rect = self.zone.get_rect(center=(XC, YC))
        self.zone.fill(BLUE)
        self.container.blit(self.zone, self.rect)

    def updt_img(self, arr, vmin=False, vmax=False, pwr=1.0):
        arr2 = arr.astype('float')**pwr

        if not vmin:
            mmin,mmax = arr2.min(), arr2.max()
        else:
            mmin,mmax = vmin,vmax

        arr2 -= mmin
        if np.abs(mmax-mmin) > 1e-4:
            arr2 /= (mmax - mmin)

        test = self.cmap(arr2)
        self.img = (255*test[:,:,:3]).astype('int')

        pygame.surfarray.blit_array(self.zone, self.img)
        self.container.blit(self.zone, self.rect)
        pygame.draw.rect(self.container, self.FGCOLOR, self.rect, 2)

# ==================================================
class pygplot:
    ''' ------------------------------------------------
    1D plot in a pygame GUI container (typically screen)
    ------------------------------------------------ '''
    # ================
    def __init__(self, container, XW, YW, XC, YC, COLOR):
        self.XW      = XW
        self.YW      = YW
        self.LW      = 1     # default line width
        self.BGCOLOR = COLOR # default color theme
        self.FGCOLOR = WHITE # default color theme
        self.PLT_CLR = BLUE  # default color theme
        # ---
        self.container = container
        self.zone = pygame.Surface((XW, YW))
        self.rect = self.zone.get_rect(center=(XC, YC))
        self.zone.fill(COLOR)
        self.container.blit(self.zone, self.rect)

        self.xmin =  0.0
        self.xmax =   XW
        self.ymin =  0.0
        self.ymax = 10.0
        self.xrng = self.xmax - self.xmin # for faster computation?
        self.yrng = self.ymax - self.ymin

        self.fsize = np.round(self.XW/20)
        self.font = pygame.font.SysFont("monospace", self.fsize)

    def xstretch(self, v):
        return((v - self.xmin) * self.XW / self.xrng)

    def ystretch(self, v): # includes the flip to make the plot up-right!
        return(self.YW - (v - self.ymin) * self.YW / self.yrng)

    # --------------------------------------------------------------
    def updt_plot(self, (xx, yy), PLT_CLR, LW=False, DECORATE=False, RANGE=False):
        XW = self.XW
        YW = self.YW
        self.PLT_CLR = PLT_CLR
        FGCOLOR = self.FGCOLOR
        BGCOLOR = self.BGCOLOR
        if LW:
            self.LW = LW
        if DECORATE:
            self.zone.fill(BGCOLOR)

        if not RANGE:
            self.xmin = min(xx)
            self.xmax = max(xx)
            self.ymin = min(yy)
            self.ymax = max(yy)
            self.xrng = self.xmax - self.xmin
            self.yrng = self.ymax - self.ymin

        xx1 = map(self.xstretch, xx)
        yy1 = map(self.ystretch, yy)

        pygame.draw.lines(self.zone, PLT_CLR, False, zip(*[xx1, yy1]), LW)

        if DECORATE:
            pygame.draw.line(self.zone, FGCOLOR, (0, YW/2), (XW, YW/2), 1)
            pygame.draw.line(self.zone, FGCOLOR, (0, YW/4), (XW, YW/4), 1)
            pygame.draw.line(self.zone, FGCOLOR, (0,3*YW/4), (XW,3*YW/4), 1)
            # --
            lbl = self.font.render("%+.1f" % (self.ymax), True, FGCOLOR)#, BGCOLOR)
            pos = lbl.get_rect(center=(lbl.get_width()/2,self.fsize/2))
            self.zone.blit(lbl, pos)
            # --
            lbl = self.font.render("%+.1f" % (self.ymin), True, FGCOLOR)#, BGCOLOR)
            pos = lbl.get_rect(center=(lbl.get_width()/2,self.YW - self.fsize/2))
            self.zone.blit(lbl, pos)

        self.container.blit(self.zone, self.rect)
        pygame.draw.rect(self.container, FGCOLOR, self.rect, 2)
