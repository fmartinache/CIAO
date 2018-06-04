#!/usr/bin/env python

import os

class Cam():
    ''' -----------------------------------------------------------------------
    Simple camera connection interface (via the fifo)

    ----------------------------------------------------------------------- '''

    # =========================================================================
    def __init__(self, fifo_dir="/home/ciaodev/bin/"):
        self.fifo_in  = fifo_dir + "ixon_fifo_in"
        self.fifo_out = fifo_dir + "ixon_fifo_out"
        self.connected = False
        self.streaming = False
        
        if os.path.exists(self.fifo_in):
            try:
                self.cmd_fifo = open(self.fifo_in, 'w')
                self.connected = True
            except:
                print("could not open the fifo in write mode")
        else:
            print("expected fifo does not exist")

        # exposure time control: only a finite number of possibilities
        self.cam_tints = [0.00001, 0.00002, 0.00005,
                          0.00010, 0.00020, 0.00050,
                          0.00100, 0.00200, 0.00500, 
                          0.01000, 0.02000, 0.05000,
                          0.10000, 0.20000, 0.50000]

        self.cam_itint = 6
        self.cam_tint = 0.001

        if self.connected:
            self.set_tint(self.cam_tints[self.cam_itint])

    # =========================================================================
    def stream(self,):
        if self.connected:
            self.cmd_fifo.write("stream")
            self.cmd_fifo.flush()
            self.streaming = True
        
    # =========================================================================
    def pause(self,):
        self.streaming = False
        if self.connected:
            self.cmd_fifo.write("abort")
            self.cmd_fifo.flush()
        
    # =========================================================================
    def quit(self,):
        self.streaming = False
        self.connected = False

        if self.connected:
            self.cmd_fifo.write("quit")
            self.cmd_fifo.flush()
        
    # =========================================================================
    def set_tint(self, tint):
        if self.streaming is True:
            was_streaming = True
            self.pause()
        self.cmd_fifo.write("tint %.4f" % (tint,))
        self.cmd_fifo.flush()
        if was_streaming:
            self.stream()
            self.streaming = True

    # =========================================================================
    def tint_dec(self,):
        self.cam_itint = max(self.cam_itint-1, 0)
        self.cam_tint = self.cam_tints[self.cam_itint]
        if self.connected:
            self.set_tint(self.cam_tint)
        
    # =========================================================================
    def tint_inc(self,):
        self.cam_itint = min(self.cam_itint+1, len(self.cam_tints)-1)
        self.cam_tint = self.cam_tints[self.cam_itint]
        if self.connected:
            self.set_tint(self.cam_tint)

# ==========================================================
# ==========================================================
if __name__ == "__main__":
    test = Cam(fifo_dir="./")
    if test.connected:
        test.stream()
    else:
        print("Not connected to camera fifo")
