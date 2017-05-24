#!/usr/bin/python
# -*- encoding: utf8 -*-

import sys
import numpy as np
from .inout import read_nc, read_hgrid_dimensions
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
    def __init__(self, hgrid_in = None, vgrid_in = None, vdom_in={} , hdom_in={}, verbose=1):

        self._verbose = verbose
        
        #
        # horizontal global grids
        #
        self._flag_hgrid_uniform = False
        if hgrid_in is None or isinstance(hgrid_in,dict):
            # uniform grid
            self._flag_hgrid_uniform = True            
            #
            hgrid = {'Lx':3.e2*1.e3, 'Ly':2e2*1.e3, 'Nx':150, 'Ny':100}
            hgrid.update(hgrid_in)
            self._build_hgrid_uniform(**hgrid)

        else:
            # curvilinear grid
            #print '!!! need to determine Nx and Ny from files'
            self._build_hgrid_curvilinear(hgrid_in)

        #   
        # vertical grid
        #
        self._flag_vgrid_uniform = False
        if vgrid_in is None or isinstance(vgrid_in,dict):
            # uniform grid
            self._flag_vgrid_uniform = True
            #
            vgrid = {'H':4.e3, 'Nz':10}
            vgrid.update(vgrid_in)
            self._build_vgrid_uniform(**vgrid)
        else:
            # curvilinear grid
            #print '!!! need to determine Nz from files'
            self._build_vgrid_stretched(vgrid_in)


        # check that Nx, Ny, Nz can be derived or that they are provided
        if 'Nx' in hdom_in.keys():
            self.Nx=hdom_in['Nx']
        else:
            if not hasattr(self,'Nx'):
                try:
                    self.Nx = hdom_in['iend']-hdom_in['istart']+1
                except:
                    print '!!! you need to prescribe one of the two variables: Nx, iend'
                    sys.exit()
        if 'Ny' in hdom_in.keys():
            self.Ny=hdom_in['Ny']      
        else:      
            if not hasattr(self,'Ny'):
                try:
                    self.Ny = hdom_in['jend']-hdom_in['jstart']+1
                except:
                    print '!!! you need to prescribe one of the two variables: Ny, jend'
                    sys.exit()
        if 'Nz' in vdom_in.keys():
            self.Nz=vdom_in['Nz']
        else:
            if not hasattr(self,'Nz'):
                try:
                    self.Nz = vdom_in['kup']-vdom_in['kdown']+1
                except:
                    print '!!! you need to prescribe one of the two variables: Nz, kup'
                    sys.exit()


        #
        # deals with subdomains
        #
        self._flag_vdom=False
        if vdom_in:
            self._flag_vdom=True
        vdom = {'kdown': 0, 'kup': self.Nz-1, 'k0': 0}
        vdom.update(vdom_in)
        for key, value in vdom.items():
            exec ('self.' + key + '=' + str(value))

        self._flag_hdom=False
        if hdom_in:
            self._flag_hdom=True
        hdom = {'istart': 0, 'iend': self.Nx-1, 'i0': 0,
                'jstart': 0, 'jend': self.Ny-1, 'j0': 0}
        hdom.update(hdom_in)
        for key, value in hdom.items():
            exec ('self.' + key + '=' + str(value))

        # fills in jend, iend, or kup if necessary
        if 'iend' not in hdom_in:
            try:
                self.iend=self.istart+self.Nx-1
            except:
                print '!!! iend cannot be determined'
        if 'jend' not in hdom_in:
            try:
                self.jend=self.jstart+self.Ny-1
            except:
                print '!!! jend cannot be determined'
        if 'kup' not in vdom_in:
            try:
                self.kup=self.kdown+self.Nz-1
            except:
                print '!!! kup cannot be determined'                
                 

        # check consistency between subdomain indices and Nx, Ny and Nz
        if self.iend-self.istart+1!=self.Nx:
            print '!!! iend-istart+1 not equal to Nx'
            sys.exit()
        elif self.jend-self.jstart+1!=self.Ny:
            print '!!! jend-jstart+1 not equal to Ny'
        #     sys.exit()
        elif self.kup-self.kdown+1!=self.Nz:
            print '!!! kup-kdown+1 not equal to Nz'
            sys.exit()



    #
    # Uniform grids
    #
    def _build_hgrid_uniform(self, **kwargs):
        for key, value in kwargs.items():
            setattr(self, key, value)
        self.hgrid_file = None
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
        # loads dimensions for dmda creation
        #self.Nx0, self.Ny0 = read_hgrid_dimensions(self.hgrid_file)
    
    
    def load_metric_terms(self, da, comm):
        
        # create a 3D vector containing metric terms
        self.D = da.createGlobalVec()
        # load curvilinear metric terms
        v = da.getVecArray(self.D)
        v[:] = 0.
        (xs, xe), (ys, ye), (zs, ze) = da.getRanges()
        # indexes along the third dimension of 
        self._k_dxt =zs
        self._k_dyt =zs+1
        self._k_dxu =zs+2
        self._k_dyu =zs+3
        self._k_dxv =zs+4
        self._k_dyv =zs+5
        self._k_lon=zs+6
        self._k_lat=zs+7


        # Initialize xt,yt,dxt,dyt
        if self.hgrid_file is None:
            # roms input
            v[:, :, self._k_dxt] = self.dx
            v[:, :, self._k_dyt] = self.dy                   
            v[:, :, self._k_lon] = i*self.dx
            v[:, :, self._k_lat] = j*self.dy        
                    
        else:
            # open and read netcdf file
            rootgrp = Dataset(self.hgrid_file, 'r')

            # curvilinear metric
            v[:,:, self._k_dxt] = np.transpose(rootgrp.variables['dxt'][ys+self.j0:ye+self.j0,xs+self.i0:xe+self.i0],(1,0))
            v[:,:, self._k_dyt] = np.transpose(rootgrp.variables['dyt'][ys+self.j0:ye+self.j0,xs+self.i0:xe+self.i0],(1,0))
            v[:,:, self._k_lon] = np.transpose(rootgrp.variables['lon'][ys+self.j0:ye+self.j0,xs+self.i0:xe+self.i0],(1,0))
            v[:,:, self._k_lat] = np.transpose(rootgrp.variables['lat'][ys+self.j0:ye+self.j0,xs+self.i0:xe+self.i0],(1,0))
            try:
                v[:, :, self._k_dxu] = np.transpose(rootgrp.variables['dxu'][ys+self.j0:ye+self.j0,xs+self.i0:xe+self.i0],(1,0))
            except:
                print '!!! must init dxu'
                sys.exit()
            try:
                v[:, :, self._k_dyu] = np.transpose(rootgrp.variables['dyu'][ys+self.j0:ye+self.j0,xs+self.i0:xe+self.i0],(1,0))
            except:                        
                print '!!! must init dyu' 
                sys.exit()    
            try:
                v[:, :, self._k_dxv] = np.transpose(rootgrp.variables['dxv'][ys+self.j0:ye+self.j0,xs+self.i0:xe+self.i0],(1,0))
            except:
                print '!!! must init dxv ' 
                sys.exit()    
            try:
                v[:, :, self._k_dyv] = np.transpose(rootgrp.variables['dyv'][ys+self.j0:ye+self.j0,xs+self.i0:xe+self.i0],(1,0))
            except:
                print '!!!  must init dyv' 
                sys.exit()
 
        rootgrp.close()


        # Initialize dxu,dyu,dyv,dyv
        # v[xs:xe-1,:, self._k_dxu] = 0.5*(v[xs:xe-1,:, self._k_dxt]+v[xs+1:xe,:, self._k_dxt])
        # v[xe-1,:, self._k_dxu] = v[xe-2,:, self._k_dxu]
        # v[:,:, self._k_dyu] = v[:,:, self._k_dyt]

        # v[:,:, self._k_dxv] = v[:,:, self._k_dxt]
        # v[:,ys:ye-1, self._k_dyv] = 0.5*(v[:,ys:ye-1, self._k_dyt]+v[:,ys+1:ye, self._k_dyt])
        # v[:,ye-1, self._k_dyv] = v[:,ye-2, self._k_dyv]
        

        if self._flag_vgrid_uniform:
            self.zt = np.ones(self.Nz)
            self.zw = np.ones(self.Nz)
            for k in range(zs,ze):
                self.zt[k]=(k-0.5)*self.dz
                self.zw[k]=k*self.dz
        else:
            # open netdc file
            rootgrp = Dataset(self.vgrid_file, 'r')
            self.zt = rootgrp.variables['zt'][zs+self.k0:ze+self.k0]
            self.zw = rootgrp.variables['zw'][zs+self.k0:ze+self.k0]
            try:
                self.dzt = rootgrp.variables['dzt'][zs+self.k0:ze+self.k0] 
            except:
                print '!!! must init dzt ' 
                sys.exit()   
            try:
                self.dzw = rootgrp.variables['dzw'][zs+self.k0:ze+self.k0] 
            except:
                print '!!! must init dzw ' 
                sys.exit()   

            rootgrp.close()

        # self.dzt = np.diff(self.zw)
        # self.dzw = np.diff(self.zt)

        comm.barrier()
        pass
    
    def load_coriolis_parameter(self, coriolis_file, da, comm):
        v = da.getVecArray(self.D)
        (xs, xe), (ys, ye), (zs, ze) = da.getRanges()
        # indexes along the third dimension 
        self._k_f=zs+9       
        # open and read netcdf file
        rootgrp = Dataset(coriolis_file, 'r')
        v[:, :, self._k_f] = np.transpose(rootgrp.variables['f'][ys+self.j0:ye+self.j0,xs+self.i0:xe+self.i0],(1,0))
        rootgrp.close()
        #
        comm.barrier()
        pass
     
    def load_mask(self, mask_file, da, comm, mask3D=False):
        """
        load reference mask from metrics file
        input:
        - mask_file : netcdf file containning the mask
        - da : instance of data management object
        - comm : The communicator for the DMDA object da
        - mask3D: flag for 3D masks (default is False)
        output:
        - grid.D[grid._k_mask,:,:] : contains the mask
        """
        self.mask3D = mask3D
        if not mask3D:
            v = da.getVecArray(self.D)
        else:
            self.mask3D = da.createGlobalVec()
            v = da.getVecArray(self.mask3D)

        (xs, xe), (ys, ye), (zs, ze) = da.getRanges()
        # index of the mask along the third dimension 
        self._k_mask=zs+8
        if not mask3D:
            try:
                # open the netcdf file and read the mask
                rootgrp = Dataset(mask_file, 'r')
                v[:, :, self._k_mask] = np.transpose(rootgrp.variables['mask'][ys+self.j0:ye+self.j0,xs+self.i0:xe+self.i0],(1,0))
                rootgrp.close()
                if self._verbose:
                    print 'The mask is 2D and loaded'
            except:
                # no mask found, only sea
                v[:, :, self._k_mask] = 1.   
                if self._verbose:
                    print 'The mask is 2D but no data was found'
        else:
            try:
                # open the netcdf file and read the mask
                rootgrp = Dataset(mask_file, 'r')
                for k in range(zs, ze):
                    for j in range(ys, ye):
                        for i in range(xs, xe):
                            v[i, j, k] = rootgrp.variables['mask'][k+self.k0,j+self.j0,i+self.i0]               
                rootgrp.close()
                if self._verbose:
                    print 'The mask is 3D and loaded'
            except:
                # no mask found, only sea
                v[:, :, :] = 1.
                if self._verbose:
                    print 'The mask is 3D but no data was found'
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
                + '  min(dzw) = %e , mean(dzw) = %e, max(dzw) = %e \n' \
                    % (np.min(self.dzw), np.mean(self.dzw), np.max(self.dzw)) \
                + '  min(dzt) = %e , mean(dzt) = %e, max(dzt) = %e \n' \
                    % (np.min(self.dzt), np.mean(self.dzt), np.max(self.dzt))

        if self._flag_hdom:
            print 'Horizontal subdomain: (istart, iend) = (%d, %d), (jstart, jend) = (%d, %d)' \
                         %(self.istart, self.iend, self.jstart, self.jend)
                         
        if self._flag_vdom:
            print 'Vertical subdomain: kdown=%d, kup=%d' %(self.kdown, self.kup)

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
            z=self.zt
        return z
   


