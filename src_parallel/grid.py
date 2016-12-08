#!/usr/bin/python
# -*- encoding: utf8 -*-

import sys
import numpy as np
from .io import read_nc, read_hgrid_dimensions
# for curvilinear grids
from netCDF4 import Dataset

class grid(object):
    """ Grid object
    """



#
#==================== Builders ============================================
# 
    
    #
    # object init
    #
    def __init__(self, hgrid = None, vgrid = None, vdom={} , hdom={}, verbose=1):

        self.verbose = verbose
        #
        # horizontal global grids
        #
        hgrid_uniform_default = {'Lx':3.e2*1.e3, 'Ly':2e2*1.e3, 'Nx0':150, 'Ny0':100}
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

        #   
        # vertical grid
        #
        vgrid_uniform_default = {'H':4.e3, 'Nz0':10}
        self._flag_vgrid_uniform = False
        if vgrid is None or isinstance(vgrid,dict):
            self._flag_vgrid_uniform = True
            #
            vgrid_input = vgrid_uniform_default
            for key, value in vgrid.items():
                vgrid_input[key]=value
            #
            self._build_vgrid_uniform(**vgrid_input)
        else:
            self._build_vgrid_stretched(vgrid)

        self.kdown = 0
        self.kup = self.Nz0 - 1
        self.k0 = 0
        for key, value in vdom.items():
            exec ('self.' + key + '=' + str(value))
        self.kmargin = self.kdown - self.k0

        #
        self.istart = 0
        self.iend = self.Nx0 - 1
        self.i0 = 0
        self.jstart = 0
        self.jend = self.Ny0 - 1
        self.j0 = 0
        for key, value in hdom.items():
            exec ('self.' + key + '=' + str(value))
        self.imargin = self.istart - self.i0
        self.jmargin = self.jstart - self.j0

        self.Nx = min(self.Nx0, self.iend - self.istart + 1 + 2 * self.imargin)
        self.Ny = min(self.Ny0, self.jend - self.jstart + 1 + 2 * self.jmargin)
        self.Nz = min(self.Nz0, self.kup - self.kdown + 1 + 2 * self.kmargin)
        print

    #
    # Uniform grids
    #
    def _build_hgrid_uniform(self, **kwargs):
        for key, value in kwargs.items():
            setattr(self, key, value)
        self.hgrid_file = None
        # compute metric terms
        self.dx=self.Lx/(self.Nx0-1.)
        self.dy=self.Ly/(self.Ny0-1.)
        
    def _build_vgrid_uniform(self, **kwargs):
        for key, value in kwargs.items():
            setattr(self, key, value)
        # compute metric terms
        self.dz=self.H/(self.Nz0-1.)
    
    #
    # Curvilinear horizontal grid
    #
    def _build_hgrid_curvilinear(self, hgrid_file):
        # store metric file but metric terms are loaded later
        self.hgrid_file = hgrid_file
        # loads dimensions for dmda creation
        self.Nx0, self.Ny0 = read_hgrid_dimensions(self.hgrid_file)
        
    def load_metric_terms(self, da, comm):
        # create a 3D vector containing metric terms
        self.D = da.createGlobalVec()
        # load curvilinear metric terms
        v = da.getVecArray(self.D)
        (xs, xe), (ys, ye), (zs, ze) = da.getRanges()
        # indexes along the third dimension of 
        self._k_dx =zs
        self._k_dy =zs+1
        self._k_lon=zs+2
        self._k_lat=zs+3

        if self.hgrid_file is None:
            # roms input
            for j in range(ys, ye):
                for i in range(xs, xe):
                    v[i, j, self._k_lon] = i*self.dx
                    v[i, j, self._k_lat] = j*self.dy
            v[xs:xe, ys:ye, self._k_dx] = self.dx
            v[xs:xe, ys:ye, self._k_dy] = self.dy
        else:
            # open and read netcdf file
            rootgrp = Dataset(self.hgrid_file, 'r')

            # curvilinear metric
            for j in range(ys, ye):
                for i in range(xs, xe):
                    v[i, j, self._k_dx] = rootgrp.variables['e1'][j+self.j0,i+self.i0]
                    v[i, j, self._k_dy] = rootgrp.variables['e2'][j+self.j0,i+self.i0]
                    v[i, j, self._k_lon] = rootgrp.variables['lon'][j+self.j0,i+self.i0]
                    v[i, j, self._k_lat] = rootgrp.variables['lat'][j+self.j0,i+self.i0]
            rootgrp.close()
        #

        if self._flag_vgrid_uniform:
            for k in range(zs,ze):
                self.zc[k]=(k-0.5)*self.dz
                self.zf[k]=k*self.dz
                self.dzc = np.diff(self.zc)
                self.dzf = np.diff(self.zf)
        else:
            # open netdc file
            rootgrp = Dataset(self.vgrid_file, 'r')
            self.zc = rootgrp.variables['zc'][zs+self.k0:ze+self.k0]
            self.zf = rootgrp.variables['zf'][zs+self.k0:ze+self.k0]
            rootgrp.close()

            self.dzc = np.diff(self.zc)
            self.dzf = np.diff(self.zf)

        comm.barrier()
        pass
    
    def load_coriolis_parameter(self, coriolis_file, da, comm):
        v = da.getVecArray(self.D)
        (xs, xe), (ys, ye), (zs, ze) = da.getRanges()
        # indexes along the third dimension 
        self._k_f=zs+4       
        # open and read netcdf file
        rootgrp = Dataset(coriolis_file, 'r')
        for j in range(ys, ye):
            for i in range(xs, xe):
                v[i, j, self._k_f] = rootgrp.variables['f'][j+self.j0,i+self.i0]
        rootgrp.close()
        #
        comm.barrier()
        pass
    
    #
    # Vertically stretched grid
    #
    def _build_vgrid_stretched(self,vgrid_file):

        # store metric file but metric terms are loaded later
        self.vgrid_file = vgrid_file
        # open netcdf file
        rootgrp = Dataset(vgrid_file, 'r')
        self.Nz0 = len(rootgrp.dimensions['zc'])

        # V = read_nc(['zc','zf'], vgrid_filename)
        # self.zc = V[0]
        # self.zf = V[1]
        # #
        # self.dzc = np.diff(self.zc)
        # self.dzf = np.diff(self.zf)
        # #
        # self.Nz0 = len(self.zc)
        # if self.verbose>0 : print self.dzc



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
            # get stats about metric terms
            # not trivial to implement as min/max needs to be taken across tiles ...
            out = 'The horizontal grid is curvlinear with:\n' \
                + '  Nx = %i , Ny = %i \n' % (self.Nx, self.Ny) 
                #+ '  min(dx) = %e , mean(dx) = %e, max(dx) = %e \n' % (np.min(self.dx), np.mean(self.dx), np.max(self.dx)) \
                #+ '  min(dy) = %e , mean(dy) = %e, max(dy) = %e \n' % (np.min(self.dy), np.mean(self.dy), np.max(self.dy))
                
        if self._flag_vgrid_uniform:
            out += 'The vertical grid is uniform with:\n' \
                + '  Nz = %i' % (self.Nz) \
                + ' , H = %e m' % (self.H) \
                + ' , dz = %e \n' % (self.dz)
        else:
            out += 'The vertical grid is stretched with:\n' \
                + '  Nz = %i' % (self.Nz) \
                + '  min(dzf) = %e , mean(dzf) = %e, max(dzf) = %e \n' \
                    % (np.min(self.dzf), np.mean(self.dzf), np.max(self.dzf)) \
                + '  min(dzc) = %e , mean(dzc) = %e, max(dzc) = %e \n' \
                    % (np.min(self.dzc), np.mean(self.dzc), np.max(self.dzc))

            # print if a subdomain is considered
            if self.kdown>0 or self.kup<self.Nz0-1:
                print 'Vertical subdomain: kdown=%d, kup=%d' %(self.kdown, self.kup)
            if self.istart>0 or self.iend<self.Nx0-1 or self.jstart>0 or self.jend<self.Ny0-1:
                print 'Horizontal subdomain: (istart, iend) = (%d, %d), (jstart, jend) = (%d, %d)' \
                         %(self.istart, self.iend, self.jstart, self.jend)
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
        else:
            x=np.arange(0., float(self.Nx))
            y=np.arange(0., float(self.Ny))
        return x,y

    def get_z(self):
        if self._flag_vgrid_uniform:
            z=np.linspace(0,self.H,self.Nz)
        else:
            z=self.zc
        return z
   


