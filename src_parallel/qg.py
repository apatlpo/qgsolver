#!/usr/bin/python
# -*- encoding: utf8 -*-

from .grid import *
from .solver import *

import petsc4py
#from Cython.Compiler.Main import verbose
petsc4py.init(sys.argv)
from petsc4py import PETSc

import numpy as np


class qg():
    """ QG object
    """
    
    def __init__(self,
                 grid_uniform_input = None,
                 grid_lonlat_file = None,
                 N2 = 1e-4,
                 verbose = 1,
                 ):
        """ QG object creation
        Parameters:
        """
        
        # determine if run is parallel
        #print 'The libraries assumes the code is run in parallel'
        #self.parallel=True
        
        
        ### horizontal and vertical global grids
        grid_uniform_input_default = {'Lx':1.e7, 'Ly':1.5e7, 'H':4.e3, 
                                      'Nx':100, 'Ny':150, 'Nz':10}
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
        
        # get petsc options from command line
        OptDB = PETSc.Options()

        # determine the tile decomposition        
        #n  = OptDB.getInt('n', 16)
        #nx = OptDB.getInt('nx', n)
        #ny = OptDB.getInt('ny', n)
        #nz = OptDB.getInt('nz', n)        
        #kplt = OptDB.getInt('kplt', nz//2)
        
        # setup tiling
        #da = PETSc.DMDA().create([nx, ny, nz], stencil_width=2)
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

        ### vertical stratification
        self.N2 = N2

        ### declare petsc vectors
        # PV
        self.Q = self.da.createGlobalVec()
        # streamfunction
        self.PSI = self.da.createGlobalVec()
        # local vectors
        #self._localQ  = self.da.createLocalVec()
        #self._localPSI  = self.da.createLocalVec()
#        # for plotting purposes
#        PSIn = da.createNaturalVec()
#        Qn = da.createNaturalVec()

        ### initiate pv inversion solver
        self.pvinv=pvinversion(self)

        
    def set_psi(self, analytical_psi=True, file_psi=None):
        """ Set psi to a given value
        """
        if self._verbose:
            print 'Set psi to ... \n'
        #if analytical_psi:
            

    def set_q(self, analytical_q=True, file_q=None):
        """ Set q to a given value
        """
        #""" Initialized potential vorticity along with boundary conditions """
        #
        #
        q = self.da.getVecArray(self.Q)
        mx, my, mz = self.da.getSizes()
        (xs, xe), (ys, ye), (zs, ze) = self.da.getRanges()
        #print 'xs=', xs, 'xe=', xe, 'mx=', mx, 'ys=', ys, 'ye=', ye, 'my=', my, 'zs=', zs, 'ze=', ze, 'mz=', mz
        #
        if analytical_q:
            if self._verbose:
                print 'Set q analytically \n'
            for k in range(zs, ze):
                for j in range(ys, ye):
                    for i in range(xs, xe):
                        q[i, j, k] = 10.*np.exp(-((i/float(mx-1)-0.5)**2 
                                                  + (j/float(my-1)-0.5)**2)/0.1**2)
                        q[i, j, k] *= np.sin(i/float(mx-1)*np.pi) 
                        q[i, j, k] *= np.sin(2*j/float(my-1)*np.pi) 
    
    def setup_time_stepping(self):
        """ Setup the time stepping
        Create additional vectors
        """
        pass
        



