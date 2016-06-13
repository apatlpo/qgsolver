#!/usr/bin/python
# -*- encoding: utf8 -*-

import numpy as np

class grid_uniform():
    """ Uniform grid object
    """
        
    def __init__(self, **kwargs):
        # fill in user provided parameters
        for key, value in kwargs.items():
            setattr(self, key, value)
        # compute metric terms
        self.dx=self.Lx/(self.Nx-1.)
        self.dy=self.Ly/(self.Ny-1.)
        self.dz=self.H/(self.Nz-1.)
        # print out grid parameters ... moved to qg.__init__
        #print 'Lx=%e ,Ly= %e, H=%e \n' % (self.Lx, self.Ly, self.H)
        #print 'Nx=%i ,Ny= %i, Nz=%i \n' % (self.Nx, self.Ny, self.Nz)
        #print 'dx=%e ,dy= %e, dz=%e \n' % (self.dx, self.dy, self.dz)
        
    def get_xyz(self):
        x=np.linspace(0,self.Lx,self.Nx)
        y=np.linspace(0,self.Ly,self.Ny)
        z=np.linspace(0,self.H,self.Nz)
        return x,y,z
    
class grid_lonlat():
    """ lon/lat grid object, parallel run
    """
    
    def __init__(self):
        pass
  
    