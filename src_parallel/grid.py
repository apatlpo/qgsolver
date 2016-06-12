#!/usr/bin/python
# -*- encoding: utf8 -*-


class grid_uniform():
    """ Uniform grid object
    """
        
    def __init__(self, **kwargs):
        # fill in user provided parameters
        for key, value in kwargs.items():
            setattr(self, key, value)
        # compute metric terms
        self.dx=self.Lx/self.Nx
        self.dy=self.Ly/self.Ny
        self.dz=self.H/self.Nz
        # print out grid parameters ... moved to qg.__init__
        #print 'Lx=%e ,Ly= %e, H=%e \n' % (self.Lx, self.Ly, self.H)
        #print 'Nx=%i ,Ny= %i, Nz=%i \n' % (self.Nx, self.Ny, self.Nz)
        #print 'dx=%e ,dy= %e, dz=%e \n' % (self.dx, self.dy, self.dz)
        

class grid_lonlat():
    """ lon/lat grid object, parallel run
    """
    
    def __init__(self):
        pass
  
    