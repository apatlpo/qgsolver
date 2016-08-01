#!/usr/bin/python
# -*- encoding: utf8 -*-

from .grid import *
from .solver import *

import petsc4py
#from Cython.Compiler.Main import verbose
petsc4py.init(sys.argv)
from petsc4py import PETSc

import numpy as np
from .io import read_nc_petsc
#from netCDF4 import Dataset

class qg():
    """ QG object
    """
    
    def __init__(self,
                 grid_uniform_input = None,
                 grid_lonlat_file = None,
                 N2 = 1e-3, f0=7e-5, K=1.e2,
                 dt = 86400.e-1,
                 verbose = 1,
                 ):
        """ QG object creation
        Parameters:
        """
                
        ### horizontal and vertical global grids
        grid_uniform_input_default = {'Lx':3.e2*1.e3, 'Ly':2e2*1.e3, 'H':4.e3, 
                                      'Nx':150, 'Ny':100, 'Nz':10}
        grid_input = grid_uniform_input_default
        #if grid_uniform_input is None:
        #    grid_uniform_input = grid_uniform_input_default
        #else:
        #for key in grid_uniform_input:
        #    grid_input[key]=grid_uniform_input[key]
        for key, value in grid_uniform_input.items():
            grid_input[key]=value
        if grid_lonlat_file is None:
            #print 'The grid is uniform'
            self.flag_grid_uniform = True
            self.grid = grid_uniform(**grid_input)
        else:
            print 'The grid is in lon/lat space, parallel run'
            print 'The grid file is: '+grid_lonlat_file
            self.grid = grid_lonlat(grid_lonlat_file)


        ### init petsc

        #OptDB = PETSc.Options()

        # setup tiling
        self.da = PETSc.DMDA().create([self.grid.Nx, self.grid.Ny, self.grid.Nz],
                                      stencil_width=2)
        self.comm = self.da.getComm()
        self.rank = self.comm.getRank()
        
        # for lon/lat grids should load metric terms over tiles
        # self.grid.load_metric()

        # print out grid information
        if self.rank is 0 and verbose>0:
            self._verbose=verbose
        else:
            self._verbose=0
        #
        if self._verbose>0 and self.flag_grid_uniform:
            # general information
            print '\nA QG model object is being created \n'
            # print out grid parameters
            print 'The grid is uniform with:'            
            print 'Nx=%i ,Ny= %i, Nz=%i' % (self.grid.Nx, self.grid.Ny, self.grid.Nz)
            print 'Lx=%e km ,Ly= %e km , H=%e m' % (self.grid.Lx/1e3, self.grid.Ly/1e3, self.grid.H)
            print 'dx=%e , dy= %e , dz=%e \n' % (self.grid.dx, self.grid.dy, self.grid.dz)


        ### vertical stratification and Coriolis
        self.N2 = N2
        self.f0 = f0
        self._sparam = self.f0**2/self.N2
        self.K = K

        ### declare petsc vectors
        # PV
        self.Q = self.da.createGlobalVec()
        # streamfunction
        self.PSI = self.da.createGlobalVec()

        ### initiate pv inversion solver
        self.pvinv=pvinversion(self)

        ### initiate time stepper
        self.tstepper = time_stepper(self, dt)
        
    
        
    def set_psi(self, analytical_psi=True, file_psi=None):
        """ Set psi to a given value
        """
        if self._verbose:
            print 'Set psi to ... \n'
        #if analytical_psi:
            

    def set_q(self, analytical_q=True, file_q=None):
        """ Set q to a given value
        """
        #
        if file_q is not None:
            if self._verbose:
                print 'Set q from file '+file_q+' ...\n'
            read_nc_petsc(self.Q, 'q', file_q, self)
        elif analytical_q:
            if self._verbose:
                print 'Set q analytically \n'
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
            

    def set_q_bdy(self):
        """ Reset q at boundaries such that dq/dn=0 """
        #
        q = self.da.getVecArray(self.Q)
        #
        mx, my, mz = self.da.getSizes()
        (xs, xe), (ys, ye), (zs, ze) = self.da.getRanges()
        #
        # south bdy
        if (ys==0):
            j=0
            for k in range(zs, ze):
                for i in range(xs, xe):
                    q[i, j, k] = q[i,j+1,k]
        # north bdy
        if (ye==my):
            j=my-1
            for k in range(zs, ze):
                for i in range(xs, xe):
                    q[i, j, k] = q[i,j-1,k]
        # west bdy
        if (xs==0):
            i=0
            for k in range(zs, ze):
                for j in range(ys, ye):
                    q[i, j, k] = q[i+1,j,k]
        # east bdy
        if (xe==mx):
            i=mx-1
            for k in range(zs, ze):
                for j in range(ys, ye):
                    q[i, j, k] = q[i-1,j,k]                


    def invert_pv(self):
        """ wrapper around solver solve method
        """
        self.pvinv.solve(self.Q,self.PSI,self.da)


    def tstep(self, nt=1):
        """ Time step wrapper
        """
        self.tstepper.go(self, nt)

    
    def get_uv(self):
        """ Compute horizontal velocities
        """



