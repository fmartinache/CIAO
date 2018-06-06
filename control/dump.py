''' ===========================================================================
collection of a few pieces taken out of the other modules that can be of use
later.
=========================================================================== '''

self.shm_xslp = shm('/tmp/xslp.im.shm', verbose=False)
self.shm_yslp = shm('/tmp/yslp.im.shm', verbose=False)
self.shm_phot = shm('/tmp/phot_inst.im.shm', verbose=False)
        
    # -------------------------------------------------------------------------
    def get_slopes(self, nav=20, reform=True):
        ''' -------------------------------------------------------------------
        test
        ------------------------------------------------------------------- '''
        x_cnt = self.shm_xslp.get_counter()
        x_sig = self.shm_xslp.get_data(check=x_cnt, reform=reform)
        x_cnt = self.shm_xslp.get_counter()
        
        y_cnt = self.shm_yslp.get_counter()
        y_sig = self.shm_yslp.get_data(check=y_cnt, reform=reform)
        y_cnt = self.shm_yslp.get_counter()
        
        for ii in range(nav-1):
            x_sig += self.shm_xslp.get_data(check=x_cnt, reform=reform)
            x_cnt  = self.shm_xslp.get_counter()
            
            y_sig += self.shm_yslp.get_data(check=y_cnt, reform=reform)
            y_cnt  = self.shm_yslp.get_counter()
            
        x_sig /= nav
        y_sig /= nav

        return np.concatenate((x_sig, y_sig))

