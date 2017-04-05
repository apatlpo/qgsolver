#!/usr/bin/python
# -*- encoding: utf8 -*-

import sys

from .grid import *
from .solver import *

import petsc4py
#from Cython.Compiler.Main import verbose
petsc4py.init(sys.argv)
from petsc4py import PETSc

import numpy as np
from .io import read_nc_petsc, read_nc_petsc_2D
from .io import write_nc




class qg_model():
    """ QG object
    """
    
    def __init__(self,
                 hgrid = None, vgrid=None,
                 N2 = 1e-3, f0 = 7e-5, K = 1.e2,
                 f0N2_file = None,
                 dt = 86400.e-1,
                 vdom={}, hdom={},
                 ncores_x=None, ncores_y=None,
                 bdy_type_in={},
                 verbose = 1,
                 substract_fprime=False
                 ):
        """ QG object creation
        Parameters:
        """
        self.g = 9.81
        self.rho0 = 1000.

        #
        # Build grid object
        #
        self.grid = grid(hgrid, vgrid, vdom, hdom, verbose=verbose)

        #
        # init petsc
        #
        
        # test whether tiling is consistent with dimensions
        if self.grid.Nx%ncores_x!=0 or self.grid.Ny%ncores_y!=0:
            print 'Tiling does not match dimensionts: Nx/ncores_x=%f, Ny/ncores_y=%f' \
                    %(float(self.grid.Nx)/ncores_x, float(self.grid.Ny)/ncores_y) 
            sys.exit()
            
        # setup tiling
        self.da = PETSc.DMDA().create(sizes = [self.grid.Nx, self.grid.Ny, self.grid.Nz],
                                      proc_sizes = [ncores_x,ncores_y,1],
                                      stencil_width = 2)
        # http://lists.mcs.anl.gov/pipermail/petsc-dev/2016-April/018889.html
        self.comm = self.da.getComm()
        self.rank = self.comm.getRank()
        # print tiling information
        if self.rank is 0 and verbose>0:
            print 'PETSc DMDA created'
            print 'The 3D grid is tiled according to (nproc_x, nproc_y, nproc_z) : '\
                    +str(self.da.proc_sizes) 
            #print 'rank='+str(self.rank)+' ranges='+str(self.da.ranges)


        # print out grid information
        if self.rank is 0 and verbose>0:
            self._verbose=verbose
        else:
            self._verbose=0

        #
        # finalize grid/metric loading
        #

        # for lon/lat grids should load metric terms over tiles
        if not self.grid._flag_hgrid_uniform or not self.grid._flag_vgrid_uniform:
            self.grid.load_metric_terms(self.da, self.comm)
        
        # initialize mask
        self.grid.load_mask(self.grid.hgrid_file, self.da, self.comm)

        #
        if self._verbose>0:
            # general information
            print 'A QG model object is being created'
            # print out grid parameters
            print self.grid
            # # print if a subdomain is considered
            # if self.kdown==0 or self.kup<self.grid.Nz-1:
            #     print 'Vertical subdomain: kdown=%d, kup=%d' %(self.kdown, self.kup)
            # if self.istart==0 or self.iend<self.grid.Nx-1 or self.jstart==0 or self.jend<self.grid.Ny-1:
            #     print 'Horizontal subdomain: (istart, iend) = (%d, %d), (jstart, jend) = (%d, %d)' \
            #              %(self.istart, self.iend, self.jstart, self.jend)


        #
        # vertical stratification and Coriolis
        #
        # N2 is at w points (cell faces), N2[k] is between q[k] and q[k+1]
        if f0N2_file is not None:
            if self._verbose:
                print 'Reads N2 from '+f0N2_file
            #
            self.N2 = read_nc('N2', f0N2_file, self)
        else:
            if self._verbose:
                print 'Set N2 from user prescribed value = '+str(N2)+' 1/s^2'
            #
            self.N2 = N2*np.ones(self.grid.Nz)

        #
        if f0N2_file is not None:
            if self._verbose:
                print 'Reads f0 from '+f0N2_file
            #
            self.f0 = read_nc('f0', f0N2_file, self)
            #
            if self._verbose:
                print 'Reads Coriolis parameter f from '+f0N2_file
            self.grid.load_coriolis_parameter(f0N2_file, self.da, self.comm)
        else:
            self.f0 = f0
            if self._verbose:
                print 'Sets f0 to %.3e' %f0
                
        #
        self._sparam = self.f0**2/self.N2
        self.K = K
        #
        # declare petsc vectors
        #
        # PV
        self.Q = self.da.createGlobalVec()
        # streamfunction
        self.PSI = self.da.createGlobalVec()
        self.PSI_BG = self.da.createGlobalVec()
        self.PSI_OT = self.da.createGlobalVec()
        # density
        self.RHO = self.da.createGlobalVec()

        # default top and bottom boudary condition = 'N' pour Neumann. 
        # Other possibility 'D' for Direchlet
        self.bdy_type = {'top':'N','bottom':'N','lateral':'D'}
        self.bdy_type.update(bdy_type_in)

        # initiate pv inversion solver
        self.pvinv = pvinversion(self, substract_fprime=substract_fprime)

        # initiate time stepper
        #
        self.tstepper = time_stepper(self, dt)
        #print 'debug: time stepper object created'


    def tstep(self, nt=1):
        """ Time step wrapper
        """
        self.tstepper.go(self, nt)

    
    def get_uv(self):
        """ Compute horizontal velocities
        """

    def set_identity(self):
    	ONE = self.da.createGlobalVec()
        one = self.da.getVecArray(ONE)
        (xs, xe), (ys, ye), (zs, ze) = self.da.getRanges()
        one[:] = 1.
        return ONE

