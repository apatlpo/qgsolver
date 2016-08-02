#!/usr/bin/python
# -*- encoding: utf8 -*-

import numpy as np
    
class grid(object):
    """ Grid object
    """
    
#
#==================== BUilders ============================================
# 
    
    def __init__(self, hgrid = None, vgrid = None):
        
        ### horizontal global grids
        hgrid_uniform_default = {'Lx':3.e2*1.e3, 'Ly':2e2*1.e3, 'H':4.e3, 
                                      'Nx':150, 'Ny':100, 'Nz':10}
        self._flag_hgrid_uniform = False
        if hgrid is None or isinstance(hgrid,dict):
            #
            self._flag_hgrid_uniform = True            
            #
            hgrid_input = hgrid_uniform_default
            for key, value in hgrid.items():
                hgrid_input[key]=value
            #
            self._build_hgrid_uniform(**hgrid_input)
        else:
            self._build_hgrid(hgrid)
            
        ### vertical grid
        vgrid_uniform_default = {'H':4.e3, 'Nz':10}
        self._flag_vgrid_uniform = False
        if vgrid is None or isinstance(vgrid,dict):
            self._flag_vgrid_uniform = True
            #
            vgrid_input = vgrid_uniform_default
            for key, value in vgrid.items():
                hgrid_input[key]=value
            #
            self._build_vgrid_uniform(**vgrid_input)
        else:
            self._build_vgrid_stretched(vgrid)

    def _build_hgrid_uniform(self, **kwargs):
        for key, value in kwargs.items():
            setattr(self, key, value)
        # compute metric terms
        self.dx=self.Lx/(self.Nx-1.)
        self.dy=self.Ly/(self.Ny-1.)
        
    def _build_vgrid_uniform(self, **kwargs):
        for key, value in kwargs.items():
            setattr(self, key, value)
        # compute metric terms
        self.dz=self.H/(self.Nz-1.)
        
#
#==================== Grid information ============================================
#             
              
    def __str__(self):
        if self._flag_hgrid_uniform:
            out = 'The horizontal grid is uniform with:\n' \
                + '  Nx = %i , Ny = %i \n' % (self.Nx, self.Ny) \
                + '  Lx = %e km , Ly = %e km \n' % (self.Lx/1e3, self.Ly/1e3) \
                + '  dx = %e , dy = %e \n' % (self.dx, self.dy)
        if self._flag_hgrid_uniform:
            out += 'The vertical grid is uniform with:\n' \
                + '  Nz = %i' % (self.Nz) \
                + ' , H = %e m' % (self.H) \
                + ' , dz = %e \n' % (self.dz)
        return out
      
                  
#
#==================== extract grid data ============================================
#

    def get_xyz(self):
        x,y = self.get_xy()
        z = self.get_z()
        return x,y,z
                  
    def get_xy(self):
        if self._flag_hgrid_uniform:
            x=np.linspace(0,self.Lx,self.Nx)
            y=np.linspace(0,self.Ly,self.Ny)
        return x,y

    def get_z(self):
        if self._flag_vgrid_uniform:
            z=np.linspace(0,self.H,self.Nz)
        return z
   
