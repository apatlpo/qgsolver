#!/usr/bin/python
# -*- encoding: utf8 -*-

import sys
import numpy as np
from .inout import read_nc, read_hgrid_dimensions
# for curvilinear grids

from netCDF4 import Dataset
import xarray as xr

class grid(object):
    """
    Grid object
    """

#
#==================== Builders ============================================
#

    #
    # object init
    #
    def __init__(self,
                 hgrid,
                 vgrid,
                 hdom,
                 vdom,
                 mask=False,
                 verbose=1,
                 ):
        """ Builds a grid object

        Parameters
        ----------
        hgrid : str, dict or None
            horizontal grid file name or analytical grid if dict or None
            Example:
            hgrid = {'Lx':300.*1.e3, 'Ly':200.*1.e3, 'Nx':150, 'Ny':100}
        vgrid_in : str, dict or None
            vertical grid file name or analytical grid if dict or None
            Example:
            vgrid = {'H':4.e3, 'Nz':10}
        hdom : dict
            horizontal grid dimension description
            Example:
            hdom = {'Nx': 100, 'Ny': 200}
            hdom = {'Nx': 100, 'Ny': 200, 'i0': 10, 'j0': 20}
            i0 and j0 are start indices in grid input netcdf file
            missing parameters are deduced but one should have: Nx=iend-istart+1, Ny=jend-jstart+1
        vdom : dict
            vertical grid dimension description
            Example:
            vdom = {'Nz': 10, 'k0': 10}
            k0 is the start index in grid input netcdf file
            missing parameters are deduced but one should have: kup-kdown+1
        mask : boolean, optional
            activates the use of a mask, default is false
        verbose : int, optional
            degree of verbosity, 0 means no outputs, default is 1

        """
        self._verbose = verbose

        #
        # horizontal global grids
        #
        self._flag_hgrid_uniform = False
        if hgrid is None or isinstance(hgrid,dict):
            # uniform grid
            self._flag_hgrid_uniform = True
            #
            _hgrid = {'Lx':300.*1.e3, 'Ly':200.*1.e3, 'Nx':200, 'Ny':100}
            _hgrid.update(hgrid)
            self._build_hgrid_uniform(**_hgrid)
        else:
            # curvilinear grid
            #print('!!! need to determine Nx and Ny from files')
            self._build_hgrid_curvilinear(hgrid)
        # mask
        self.mask=mask

        #
        # vertical grid
        #
        self._flag_vgrid_uniform = False
        if vgrid is None or isinstance(vgrid,dict):
            # uniform grid
            self._flag_vgrid_uniform = True
            #
            _vgrid = {'H':4.e3, 'Nz':10}
            _vgrid.update(vgrid)
            self._build_vgrid_uniform(**_vgrid)
        else:
            # curvilinear grid
            #print('!!! need to determine Nz from files')
            self._build_vgrid_stretched(vgrid)


        # check that Nx, Ny, Nz can be derived or that they are provided
        if 'Nx' in hdom:
            self.Nx=hdom['Nx']
        else:
            if not hasattr(self,'Nx'):
                try:
                    self.Nx = hdom['iend']-hdom['istart']+1
                except:
                    print('!Error: you need to prescribe one of the two variables: Nx, iend')
                    sys.exit()
        if 'Ny' in hdom:
            self.Ny=hdom['Ny']
        else:
            if not hasattr(self,'Ny'):
                try:
                    self.Ny = hdom['jend']-hdom['jstart']+1
                except:
                    print('!Error: you need to prescribe one of the two variables: Ny, jend')
                    sys.exit()
        if 'Nz' in vdom:
            self.Nz=vdom['Nz']
        else:
            if not hasattr(self,'Nz'):
                try:
                    self.Nz = vdom['kup']-vdom['kdown']+1
                except:
                    print('!Error: you need to prescribe one of the two variables: Nz, kup')
                    sys.exit()

        #
        # deals with subdomains
        #
        self._flag_vdom=True if vdom else False
        _vdom = {'kdown': 0, 'kup': self.Nz-1, 'k0': 0}
        _vdom.update(vdom)

        self._flag_hdom=True if hdom else False
        _hdom = {'istart': 0, 'iend': self.Nx-1, 'i0': 0,
                'jstart': 0, 'jend': self.Ny-1, 'j0': 0}
        _hdom.update(hdom)

        # set as attributes
        for key, value in {**_hdom, **_vdom}.items():
            setattr(self, key, value)

        # fills in jend, iend, or kup if necessary
        if 'iend' not in hdom:
            try:
                self.iend=self.istart+self.Nx-1
            except:
                print('!Error: iend cannot be determined')
        if 'jend' not in hdom:
            try:
                self.jend=self.jstart+self.Ny-1
            except:
                print('!Error: jend cannot be determined')
        if 'kup' not in vdom:
            try:
                self.kup=self.kdown+self.Nz-1
            except:
                print('!Error: kup cannot be determined')

        # check consistency between subdomain indices and Nx, Ny and Nz
        if self.iend-self.istart+1!=self.Nx:
            print('!Error: iend-istart+1 not equal to Nx')
            sys.exit()
        elif self.jend-self.jstart+1!=self.Ny:
            print('!Error: jend-jstart+1 not equal to Ny')
            sys.exit()
        elif self.kup-self.kdown+1!=self.Nz:
            print('!Error: kup-kdown+1 not equal to Nz')
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
        #self.zt = np.ones(self.Nz)
        #self.zw = np.ones(self.Nz)
        self.dz = self.H/(self.Nz-1.)
        self.dzt = np.ones(self.Nz) * self.dz
        self.dzw = np.ones(self.Nz) * self.dz
        self.zt = np.cumsum(self.dzt) - .5*self.dz
        self.zw = np.cumsum(self.dzw) - self.dz

    #
    # Curvilinear horizontal grid
    #
    def _build_hgrid_curvilinear(self, hgrid_file):
        # store metric file but metric terms are loaded later
        self.hgrid_file = hgrid_file
        # loads dimensions for dmda creation
        #self.Nx0, self.Ny0 = read_hgrid_dimensions(self.hgrid_file)


    def load_metric_terms(self, da):
        """ Load metric terms from self.hgrid_file

        Parameters
        ----------
        da : petsc DMDA
            holds the petsc grid

        """

        #if not hasattr(self, D):
        # create a 3D vector containing metric terms
        #self.D = da.createGlobalVec()

        # load curvilinear metric terms
        v = da.getVecArray(self.D)
        v[:] = 0.
        (xs, xe), (ys, ye), (zs, ze) = da.getRanges()
        # indexes along the third dimension of
        self._k_dxt = zs
        self._k_dyt = zs+1
        self._k_dxu = zs+2
        self._k_dyu = zs+3
        self._k_dxv = zs+4
        self._k_dyv = zs+5
        self._k_lon = zs+6
        self._k_lat = zs+7


        # Initialize xt,yt,dxt,dyt
        if self.hgrid_file is None:
            # roms input
            v[:, :, self._k_dxt] = self.dx
            v[:, :, self._k_dyt] = self.dy
            lx, ly = np.meshgrid(self.dx*np.arange(xs,xe,1.),self.dy*np.arange(ys,ye,1.),indexing='ij')
            v[:, :, self._k_lon] = lx
            v[:, :, self._k_lat] = ly

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
                print('!Error: must init dxu')
                sys.exit()
            try:
                v[:, :, self._k_dyu] = np.transpose(rootgrp.variables['dyu'][ys+self.j0:ye+self.j0,xs+self.i0:xe+self.i0],(1,0))
            except:
                print('!Error: must init dyu')
                sys.exit()
            try:
                v[:, :, self._k_dxv] = np.transpose(rootgrp.variables['dxv'][ys+self.j0:ye+self.j0,xs+self.i0:xe+self.i0],(1,0))
            except:
                print('!Error: must init dxv ')
                sys.exit()
            try:
                v[:, :, self._k_dyv] = np.transpose(rootgrp.variables['dyv'][ys+self.j0:ye+self.j0,xs+self.i0:xe+self.i0],(1,0))
            except:
                print('!Error:  must init dyv')
                sys.exit()

        rootgrp.close()

        #if self._flag_vgrid_uniform:
        #    self.zt = np.ones(self.Nz)
        #    self.zw = np.ones(self.Nz)
        #    self.dzt = np.ones(self.Nz)*self.dz
        #    self.dzw = np.ones(self.Nz)*self.dz
        #    for k in range(zs,ze):
        #        self.zt[k]=(k-0.5)*self.dz
        #        self.zw[k]=k*self.dz
        #else:
        # open netdc file
        rootgrp = Dataset(self.vgrid_file, 'r')
        self.zt = rootgrp.variables['zt'][zs+self.k0:ze+self.k0]
        self.zw = rootgrp.variables['zw'][zs+self.k0:ze+self.k0]
        try:
            self.dzt = rootgrp.variables['dzt'][zs+self.k0:ze+self.k0]
        except:
            print('!Error: must init dzt ')
            sys.exit()
        try:
            self.dzw = rootgrp.variables['dzw'][zs+self.k0:ze+self.k0]
        except:
            print('!Error: must init dzw ')
            sys.exit()

            rootgrp.close()
        #
        comm = da.getComm()
        comm.barrier()

    def load_coriolis_parameter(self, coriolis_file, da):
        """ Load the Coriolis parameter

        Parameters
        ----------
        coriolis_file : str
            netcdf file containing the Coriolis parameter
        da : petsc DMDA
            holds the petsc grid

        """
        if not hasattr(self, D):
            # create a 3D vector containing metric terms
            self.D = da.createGlobalVec()
        v = da.getVecArray(self.D)
        (xs, xe), (ys, ye), (zs, ze) = da.getRanges()
        # indexes along the third dimension
        self._k_f=zs+9
        # open and read netcdf file
        rootgrp = Dataset(coriolis_file, 'r')
        v[:, :, self._k_f] = np.transpose(rootgrp.variables['f'][ys+self.j0:ye+self.j0,xs+self.i0:xe+self.i0],(1,0))
        rootgrp.close()
        #
        comm=da.getComm()
        comm.barrier()

    def load_mask(self,
                  mask_file,
                  da,
                  mask3D=False,
                  dims=("x", "y"),
                  ):
        """Load reference mask from metrics file
        grid.D[grid._k_mask,:,:] will contain the mask

        Parameters
        ----------
        mask_file : str
            netcdf file containing the mask
        da : petsc DMDA
            holds the petsc grid
        mask3D: boolean
            flag for 3D masks, default is False

        """
        if not hasattr(self, "D"):
            # create a 3D vector containing metric terms
            self.D = da.createGlobalVec()
        if not mask3D:
            self.mask3D = False
            v = da.getVecArray(self.D)
        else:
            self.mask3D = da.createGlobalVec()
            v = da.getVecArray(self.mask3D)
        (xs, xe), (ys, ye), (zs, ze) = da.getRanges()
        # index of the mask along the third dimension
        self._k_mask=zs+8
        ds = xr.open_dataset(mask_file)
        assert "mask" in ds, \
            "No mask variable found in {}".format(mask_file)
        if not mask3D:
            _ma = (ds["mask"]
                  .isel(x=slice(xs+self.i0, xe+self.i0),
                        y=slice(ys+self.j0, ye+self.j0))
                  .transpose("x", "y")
                  )
            v[:, :, self._k_mask] = _ma
            if self._verbose>0:
                print('    The mask is 2D and loaded')
            # no mask found, only sea
            #v[:, :, self._k_mask] = 1.
        else:
            try:
                # open the netcdf file and read the mask
                rootgrp = Dataset(mask_file, 'r')
                for k in range(zs, ze):
                    for j in range(ys, ye):
                        for i in range(xs, xe):
                            v[i, j, k] = rootgrp.variables['mask'][k+self.k0,j+self.j0,i+self.i0]
                rootgrp.close()
                if self._verbose>0:
                    print('    The mask is 3D and loaded')
            except:
                # no mask found, only sea
                v[:, :, :] = 1.
                if self._verbose>0:
                    print('    The mask is 3D but no data was found')
        #
        comm=da.getComm()
        comm.barrier()
        pass

    #
    # Vertically stretched grid
    #
    def _build_vgrid_stretched(self,vgrid_file):

        # store metric file but metric terms are loaded in load_metric_terms
        self.vgrid_file = vgrid_file
        # open netcdf file
        #rootgrp = Dataset(vgrid_file, 'r')




#
#==================== Grid information ============================================
#

    def __str__(self):

        if self._flag_hgrid_uniform:
            out = '  The horizontal grid is uniform with:\n' \
                + '    Nx = %i , Ny = %i \n' % (self.Nx, self.Ny) \
                + '    Lx = %.0f km , Ly = %.0f km \n' % (self.Lx/1e3, self.Ly/1e3) \
                + '    dx = %.0f m, dy = %.0f m\n' % (self.dx, self.dy)
        else:
            # get stats about metric terms
            # not trivial to implement as min/max needs to be taken across tiles ...
            out = '  The horizontal grid is curvlinear with:\n' \
                + '    Nx = %i , Ny = %i\n' % (self.Nx, self.Ny)
            #    + '  min(dx) = %e , mean(dx) = %e, max(dx) = %e \n' % (np.min(self.dx), np.mean(self.dx), np.max(self.dx)) \
            #    + '  min(dy) = %e , mean(dy) = %e, max(dy) = %e \n' % (np.min(self.dy), np.mean(self.dy), np.max(self.dy))

        if self._flag_vgrid_uniform:
            out += '  The vertical grid is uniform with:\n' \
                +  '    Nz = %i' % (self.Nz) \
                +  ' , H = %.0f m' % (self.H) \
                +  ' , dz = %.0f m' % (self.dz)
        else:
            out += '  The vertical grid is stretched with:\n' \
                +  '    Nz = %i \n' % (self.Nz) \
                +  '    min(dzw) = %.1f m, mean(dzw) = %.1f m, max(dzw) = %.1f m\n' \
                    % (np.min(self.dzw), np.mean(self.dzw), np.max(self.dzw)) \
                +  '    min(dzt) = %.1f m, mean(dzt) = %.1f m, max(dzt) = %.1f m' \
                    % (np.min(self.dzt), np.mean(self.dzt), np.max(self.dzt))

        if self._flag_hdom:
            print('\n  Horizontal subdomain: (istart, iend) = (%d, %d), (jstart, jend) = (%d, %d)' \
                         %(self.istart, self.iend, self.jstart, self.jend))

        if self._flag_vdom:
            print('\n  Vertical subdomain: kdown=%d, kup=%d' %(self.kdown, self.kup))

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
