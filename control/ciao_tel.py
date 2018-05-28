#!/usr/bin/env python

import socket

class Tel():
    ''' -----------------------------------------------------------------------
    Simple telescope connection interface

    Primarily used to get pointing information from the telescope
    and to send tip-tilt offload commands
    ----------------------------------------------------------------------- '''

    # =========================================================================
    def __init__(self,):
        self.TCS_IP   = "10.150.10.18"
        self.TCS_PORT = 3737
        self.sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.sock.settimeout(1)
        self.connected = False
        
    # =========================================================================
    def connect(self,):
        try:
            self.sock.connect((self.TCS_IP, self.TCS_PORT))
            self.connected = True
        except:
            print("Socket connection refused")
        return self.connected
    
    # =========================================================================
    def close(self,):
        self.sock.close()
        self.connected = False
        
    # =========================================================================
    def get_hour_angle(self):        
        self.sock.send('#TCS_GET_POINTING_INFO:?\r')
        reply = sock.recv(1024)[1:-1].split(' ')
        self.tcs_hour_angle = float(reply[5])
        return self.tcs_hour_angle
    
    # =========================================================================
    def set_guiding_offsets(self, ra=0.0, dec=0.0):
        cmd = '#TCS_SET_GUIDING_OFFSETS: %.2f %.2f:?\r' % (ra, dec)
        self.sock.send(cmd)
        print(cmd)
        reply = self.sock.recv(1024)[1:-1].split(' ')
        print(reply)
        
    # =========================================================================
    def offload_tiptilt(self, ttx, tty, h0=10):
        if not self.connected:
            self.connect()

        ha = self.get_hour_angle()
        
        h0 = 10*np.pi/180.0    # mystery angle
        houra = ha + h0        # hour angle in rad

        myrot = np.array([[-np.sin(houra), -np.cos(houra)],
                          [-np.cos(houra), np.sin(houra)]])

        mycmd = np.dot(myrot, [ttx, tty]) # in arcseconds
        self.set_guiding_offsets(ra=mycmd[0], dec=mycmd[1])

# ==========================================================
# ==========================================================
if __name__ == "__main__":
    test = Tel()
    test.connect()
    if test.connected:
        print("Connected to TCS!")
    else:
        print("Not connected to TCS")
