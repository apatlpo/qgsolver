#!/usr/bin/python
# -*- encoding: utf8 -*-

import numpy as np
from .io import read_nc
 
class grid(object):
    """ Grid object
    """
    
#
#==================== Builders ============================================
# 
    
    def __init__(self, hgrid = None, vgrid = None):
        
        ### horizontal global grids
        hgrid_uniform_default = {'Lx':3.e2*1.e3, 'Ly':2e2*1.e3, 'H':4.e3, 
                                      'Nx':150, 'Ny':100, 'Nz':10}
        self._flag_hgrid_uniform = False
        if hgrid is None or isinstance(hgrid,dict):
            # uniform grid
            self._flag_hgrid_uniform = True            
            #
            hgrid_input = hgrid_uniform_default
            for key, value in hgrid.items():
                hgrid_input[key]=value
            #
            self._build_hgrid_uniform(**hgrid_input)
        else:
            # curvilinear grid
            self._build_hgrid_curvilinear(hgrid)
            
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

    #
    # Uniform grids
    #
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
    # Curvilinear horizontal grid
    #
    def _build_hgrid_curvilinear(self, hgrid_file):
        # store metric file but metric terms are loaded later
        self.hgrid_file = hgrid_file
        
    def load_metric_terms(self,):
        # load curvilinear metric terms
        pass
    
    #
    # Vertically stretched grid
    #
    def _build_vgrid_stretched(self,vgrid_filename):
        V = read_nc(['zc','zf'], vgrid_filename)
        self.zc = V[0]
        self.zf = V[0]
        self.dzc = np.diff(self.zc)
        self.dzf = np.diff(self.zf)
        self.Nz = len(self.zc)


#
#==================== Grid information ============================================
#             
              
    def __str__(self):
        
        if self._flag_hgrid_uniform:
            out = 'The horizontal grid is uniform with:\n' \
                + '  Nx = %i , Ny = %i \n' % (self.Nx, self.Ny) \
                + '  Lx = %e km , Ly = %e km \n' % (self.Lx/1e3, self.Ly/1e3) \
                + '  dx = %e , dy = %e \n' % (self.dx, self.dy)
        else:
            out = 'The horizontal grid is curvlinear with:\n' \
                + '  Nx = %i , Ny = %i \n' % (self.Nx, self.Ny) \
                + '  min(dx) = %e , mean(dx) = %e, max(dx) = %e \n' % (np.min(self.dx), np.mean(self.dx), np.max(self.dx)) \
                + '  min(dy) = %e , mean(dy) = %e, max(dy) = %e \n' % (np.min(self.dy), np.mean(self.dy), np.max(self.dy))
                
        if self._flag_hgrid_uniform:
            out += 'The vertical grid is uniform with:\n' \
                + '  Nz = %i' % (self.Nz) \
                + ' , H = %e m' % (self.H) \
                + ' , dz = %e \n' % (self.dz)
        else:
            out += 'The vertical grid is stretched with:\n' \
                + '  Nz = %i' % (self.Nz) \
                + '  min(dz) = %e , mean(dz) = %e, max(dz) = %e \n' % (np.min(self.dz), np.mean(self.dz), np.max(self.dz))
        
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
   
