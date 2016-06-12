#!/usr/bin/python
# -*- encoding: utf8 -*-


class grid_uniform():
    """ Uniform grid object
    """
    
#     # Default parameters
#     Lx = 1.e7
#     Ly = 1.5e7
#     H = 4.e3
#     Nx = 100
#     Ny = 150
#     Nz = 10
    
    def __init__(self, **kwargs):
        # fill in user provided parameters
        for key, value in kwargs.items():
            setattr(self,key, value)
        # compute metric terms
        self.dx=self.Lx/self.Nx
        self.dy=self.Ly/self.Ny
        self.dz=self.H/self.Nz
        # print out grid parameters
        print 'Lx=%e ,Ly= %e, H=%e \n' % (self.Lx, self.Ly, self.H)
        print 'Nx=%i ,Ny= %i, Nz=%i \n' % (self.Nx, self.Ny, self.Nz)
        print 'dx=%e ,dy= %e, dz=%e \n' % (self.dx, self.dy, self.dz)
        

class grid_lonlat():
    """ lon/lat grid object, parallel run
    """
    
    def __init__(self):
        pass
  
    