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
from .io import read_nc_petsc
#from netCDF4 import Dataset


class qg_model():
    """ QG object
    """
    
    def __init__(self,
                 hgrid = None, vgrid=None,
                 N2 = 1e-3, f0 = 7e-5, K = 1.e2,
                 f0N2_file = None,
                 dt = 86400.e-1,
                 kdown=0,
                 kup=5000,
                 verbose = 1,
                 ):
        """ QG object creation
        Parameters:
        """
        self.g = 9.81
        self.rho0 = 1000.

        #
        # Build grid object
        #
        self.grid = grid(hgrid, vgrid, verbose=verbose)
        self.kdown=max(kdown,1)
        self.kup=min(kup,self.grid.Nz)

        #
        # init petsc
        #
        
        #OptDB = PETSc.Options()

        # setup tiling
        #self.da = PETSc.DMDA().create([self.grid.Nx, self.grid.Ny, self.grid.Nz],
        #                              stencil_width=2)
        self.da = PETSc.DMDA().create(sizes = [self.grid.Nx, self.grid.Ny, self.grid.Nz],
                                      proc_sizes = [2,4,1],
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

        # for lon/lat grids should load metric terms over tiles
        # if not self.grid._flag_hgrid_uniform:
        self.grid.load_metric_terms(self.da, self.comm)

        # print out grid information
        if self.rank is 0 and verbose>0:
            self._verbose=verbose
        else:
            self._verbose=0
        #
        if self._verbose>0:
            # general information
            print 'A QG model object is being created'
            # print out grid parameters
            print self.grid
        #
        # vertical stratification and Coriolis
        #
        # N2 is at w points (cell faces), N2[i] is between q[i-1] and q[i]
        if f0N2_file is not None:
            if self._verbose:
                print 'Reads N2 from '+f0N2_file
            #
            self.N2 = read_nc('N2', f0N2_file)
        else:
            if self._verbose:
                print 'Set N2 from user prescribed value = '+str(N2)+' 1/s^2'
            #
            self.N2 = N2*np.ones(self.grid.Nz+1)
        #
        if f0N2_file is not None:
            if self._verbose:
                print 'Reads f0 from '+f0N2_file
            #
            self.f0 = read_nc('f0', f0N2_file)
            #
            if self._verbose:
                print 'Reads Coriolis parameter f from '+f0N2_file
            self.grid.load_coriolis_parameter(f0N2_file, self.da, self.comm)
        else:
            self.f0 = f0
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
        # density
        self.RHO = self.da.createGlobalVec()
        #
        # initiate pv inversion solver
        #
        self.pvinv = pvinversion(self)

        #
        # initiate time stepper
        #
        self.tstepper = time_stepper(self, dt)
        
    
        
    def set_psi(self, analytical_psi=True, file_psi=None):
        """
        Set psi to a given value
        """
        if file_psi is not None:
            if self._verbose:
                print 'Set psi from file '+file_psi+' ...'
            read_nc_petsc(self.PSI, 'psi', file_psi, self)
        elif analytical_psi:
            if self._verbose:
                print 'Set psi analytically '
            self.set_psi_analytically()


    def set_psi_analytically(self):
        """ Set psi analytically
        """
        psi = self.da.getVecArray(self.PSI)
        mx, my, mz = self.da.getSizes()
        (xs, xe), (ys, ye), (zs, ze) = self.da.getRanges()
        #
        if self._verbose:
            pass
            #print 'Set psi analytically \n'
        for k in range(zs, ze):
            for j in range(ys, ye):
                for i in range(xs, xe):
                    psi[i, j, k] = 0.

    def set_q(self, analytical_q=True, file_q=None):
        """ Set q to a given value
        """
        #
        if file_q is not None:
            if self._verbose:
                print 'Set q from file '+file_q+' ...'
            read_nc_petsc(self.Q, 'q', file_q, self)
        elif analytical_q:
            if self._verbose:
                print 'Set q analytically '
            self.set_q_analytically()


    def set_q_analytically(self):
        """ Set q analytically
        """
        q = self.da.getVecArray(self.Q)
        mx, my, mz = self.da.getSizes()
        (xs, xe), (ys, ye), (zs, ze) = self.da.getRanges()
        #
        if self._verbose:
            pass
            #print 'Set q analytically \n'
        for k in range(zs, ze):
            for j in range(ys, ye):
                for i in range(xs, xe):
                    q[i, j, k] = 1.e-5*np.exp(-((i/float(mx-1)-0.5)**2 
                                              + (j/float(my-1)-0.5)**2)/0.1**2)
                    q[i, j, k] *= np.sin(i/float(mx-1)*np.pi) 
                    q[i, j, k] *= np.sin(2*j/float(my-1)*np.pi)


    def set_rho(self, analytical_rho=True, file_rho=None):
        """ Set rho to a given value
        """
        #
        if file_rho is not None:
            if self._verbose:
                print 'Set rho from file '+file_rho+' ...'
            read_nc_petsc(self.RHO, 'rho', file_rho, self)
        elif analytical_rho:
            if self._verbose:
                print 'Set rho analytically '
            self.set_rho_analytically()

    def set_rho_analytically(self):
        """ Set rho analytically
        """
        rho = self.da.getVecArray(self.RHO)
        mx, my, mz = self.da.getSizes()
        (xs, xe), (ys, ye), (zs, ze) = self.da.getRanges()
        #
        if self._verbose:
            pass
            #print 'Set rho analytically \n'
        for k in range(zs, ze):
            for j in range(ys, ye):
                for i in range(xs, xe):
                    rho[i, j, k] = 0.

    def invert_pv(self):
        """ wrapper around solver solve method
        """
        self.pvinv.solve(self)


    def tstep(self, nt=1):
        """ Time step wrapper
        """
        self.tstepper.go(self, nt)

    
    def get_uv(self):
        """ Compute horizontal velocities
        """



