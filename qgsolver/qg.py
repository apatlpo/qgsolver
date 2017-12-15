#!/usr/bin/python
# -*- encoding: utf8 -*-

import sys

from .grid import *
from .state import *
from .pvinv import *
from .omegainv import *
from .timestepper import *

import petsc4py
#from Cython.Compiler.Main import verbose
petsc4py.init(sys.argv)
from petsc4py import PETSc

import numpy as np
from .inout import read_nc_petsc, read_nc_petsc_2D
from .inout import write_nc




class qg_model():
    """
    QG model
    """


#
#==================== Object init ============================================
#

    def __init__(self,
                 ncores_x=None, ncores_y=None,
                 hgrid = None, vgrid=None,
                 vdom={}, hdom={}, mask=False,
                 bdy_type_in={},
                 N2 = 1e-3, f0 = 7e-5,
                 f0N2_file = None,
                 dt = None, K = 1.e2,
                 verbose = 1,
                 substract_fprime=False,
                 flag_pvinv=True,
                 flag_omega=False
                 ):
        """
        QG model initializer

        Parameters
        ----------
        ncores_x : int
            number of MPI tilings in x direction
        ncores_y : int
            number of MPI tilings in y direction
        hgrid : dict or str
            defines horizontal grid choice
        vgrid : dict or str
            defines vertical grid choice
        bdy_type_in : dict
            may be used to turn on periodic boundary conditions {'periodic'}
        N2 : float
            Brunt Vaisala frequency
        f0 : float
            Coriolis frequency
        f0N2_file : str
            netcdf file containing N2 and f0

        """

        # useful parameters
        self.g = 9.81
        self.rho0 = 1000.

        #
        # Build grid object
        #
        self.grid = grid(hgrid, vgrid, vdom, hdom, mask, verbose=verbose)

        # set boundary conditions
        if ('periodic' in bdy_type_in.keys()) and (bdy_type_in['periodic']):
            self.BoundaryType = 'periodic'
        else:
            self.BoundaryType = None
        # default top and bottom boudary condition = 'N' pour Neumann.
        # Other possibility 'D' for Direchlet
        self.bdy_type = {'top':'N','bottom':'N'}
        self.bdy_type.update(bdy_type_in)


        #
        # init petsc
        #
        self._init_petsc(ncores_x,ncores_y)

        # print tiling information
        if self.rank is 0 and verbose>0:
            print('A QG model object is being created')
            print('  PETSc DMDA created')
            print('  The 3D grid is tiled according to (nproc_x, nproc_y, nproc_z) : '\
                    +str(self.da.proc_sizes))
            if verbose>1:
                print('rank='+str(self.rank)+' ranges='+str(self.da.ranges))

        # set verbose variable for log and debug
        if self.rank is 0 and verbose>0:
            self._verbose=verbose
        else:
            self._verbose=0
        self.grid._verbose=self._verbose

        #
        # finalize grid/metric loading
        #

        # for lon/lat grids should load metric terms over tiles
        if not self.grid._flag_hgrid_uniform or not self.grid._flag_vgrid_uniform:
            self.grid.load_metric_terms(self.da, self.comm)

        if self.grid.mask:
            # initialize mask
            self.grid.load_mask(self.grid.hgrid_file, self.da, self.comm)

        #
        if self._verbose>0:
            # print out grid parameters
            print(self.grid)
            # periodicity
            if self._verbose and self.BoundaryType is 'periodic':
                print('Boundaries are periodic')

        #
        # create an ocean state
        #
        self.state = state(self.da, N2=N2, f0=f0, f0N2_file=f0N2_file)


        # initiate pv inversion solver
        if flag_pvinv:
            self.pvinv = pvinversion(self, substract_fprime=substract_fprime)

        # initiate omega inversion
        if flag_omega:
            self.W = self.da.createGlobalVec()
            self.omegainv = omegainv(self)

        # initiate time stepper
        #self.K = K
        if dt is not None:
            self.tstepper = time_stepper(self, dt, K)



    def _init_petsc(self,ncores_x,ncores_y):
        ''' Initate Petsc environement

        Parameters
        ----------
        ncores_x: int
            Number of MPI tiles in x direction
        ncores_y: int
            Number of MPI tiles in y direction

        '''
        # test whether tiling is consistent with dimensions
        if self.grid.Nx % ncores_x != 0 or self.grid.Ny % ncores_y != 0:
            print('!Error: MPI tiling does not match dimensions: Nx/ncores_x=%f, Ny/ncores_y=%f' \
                  % (float(self.grid.Nx) / ncores_x, float(self.grid.Ny) / ncores_y))
            sys.exit()

        # setup tiling
        self.da = PETSc.DMDA().create(sizes=[self.grid.Nx, self.grid.Ny, self.grid.Nz],
                                      proc_sizes=[ncores_x, ncores_y, 1],
                                      stencil_width=2, boundary_type=self.BoundaryType)
        # http://lists.mcs.anl.gov/pipermail/petsc-dev/2016-April/018889.html

        self.comm = self.da.getComm()
        self.rank = self.comm.getRank()


#
#==================== Wrappers to set values of critical variables ===================================
#

    def set_psi(self, **kwargs):
        """ Set psi
        """
        self.state.set_psi(self, **kwargs)


    def set_q(self, **kwargs):
        """ Set q
        """
        self.state.set_q(self, **kwargs)

    def set_rho(self, **kwargs):
        """ Set rho
        """
        self.state.set_rho(self, **kwargs)

    def set_w(self, **kwargs):
        """ Set w
        """
        self.state.set_w(self, **kwargs)

    #
#==================== useful wrappers ============================================
#
                 
    def invert_pv(self):
        """ wrapper around solver solve method
        """
        self.pvinv.solve(self)

    def invert_omega(self):
        """ wrapper around solver solve method
        """
        self.omegainv.solve(self)

    def tstep(self, nt=1):
        """ Time step wrapper
        """
        self.tstepper.go(self, nt)

#
#==================== utils ============================================
#
            
    def update_rho(self, PSI=None, RHO=None):
        """ update rho from psi
        """
        
        if PSI is None:
            PSI=self.PSI
        if RHO is None:
            RHO=self.RHO
        psi = self.da.getVecArray(PSI)
        rho = self.da.getVecArray(RHO)
        
        #
        idzt = 1./self.grid.dzt
        idzw = 1./self.grid.dzw
        #
        (xs, xe), (ys, ye), (zs, ze) = self.da.getRanges()
        istart = self.grid.istart
        iend = self.grid.iend
        jstart = self.grid.jstart
        jend = self.grid.jend
        kdown = self.grid.kdown
        kup = self.grid.kup
        
        for k in range(kdown+1, kup):
            for j in range(ys, ye):
                for i in range(xs, xe):
                    rho[i,j,k] = -self.rho0*self.f0/self.g * \
                            0.5* ( (psi[i,j,k+1]-psi[i,j,k])*idzw[k] \
                                  +(psi[i,j,k]-psi[i,j,k-1])*idzw[k-1] )
        # extrapolate top and bottom
        k=kdown
        for j in range(ys, ye):
            for i in range(xs, xe):
                rho[i,j,k] = -self.rho0*self.f0/self.g * (psi[i,j,k+1]-psi[i,j,k])*idzw[k]     
        k=kup
        for j in range(ys, ye):
            for i in range(xs, xe):
                rho[i,j,k] = -self.rho0*self.f0/self.g * (psi[i,j,k]-psi[i,j,k-1])*idzw[k-1]
        return

    def get_uv(self, PSI=None):
        """ Compute horizontal velocities
        Compute U & V from Psi
        U = -dPSIdy
        V =  dPSIdx
        """

        ### create global vectors
        self._U = self.da.createGlobalVec()
        self._V = self.da.createGlobalVec()

        ### create local vectors
        local_PSI  = self.da.createLocalVec()
        local_D = self.da.createLocalVec()

        #### load vector PSI used to compute U and V
        if PSI is None:
            self.da.globalToLocal(self.PSI, local_PSI)
        else:
            self.da.globalToLocal(PSI, local_PSI)

        self.da.globalToLocal(self.grid.D, local_D)

        #
        u = self.da.getVecArray(self._U)
        v = self.da.getVecArray(self._V)
        psi = self.da.getVecArray(local_PSI)
        D = self.da.getVecArray(local_D)


        mx, my, mz = self.da.getSizes()
        (xs, xe), (ys, ye), (zs, ze) = self.da.getRanges()

        kmask = self.grid._k_mask
        kdxu = self.grid._k_dxu
        kdyu = self.grid._k_dyu
        kdxv = self.grid._k_dxv
        kdyv = self.grid._k_dyv
        kdxt = self.grid._k_dxt
        kdyt = self.grid._k_dyt

        # Initialize u=-dpsidy and v=dpsidx

        for k in range(zs,ze):
            for j in range(ys, ye):
                for i in range(xs, xe): 
                    if (i==0    or j==0 or
                        i==mx-1 or j==my-1):
                        # lateral boundaries
                        u[i, j, k] = 0.
                        v[i, j, k] = 0.
                    else:
                        u[i,j,k] = - 1. /D[i,j,kdyu] * \
                             ( 0.25*(psi[i+1,j,k]+psi[i+1,j+1,k]+psi[i,j+1,k]+psi[i,j,k]) - \
                               0.25*(psi[i+1,j-1,k]+psi[i+1,j,k]+psi[i,j,k]+psi[i,j-1,k]) )
                        v[i,j,k] =   1. /D[i,j,kdxv] * \
                             ( 0.25*(psi[i+1,j,k]+psi[i+1,j+1,k]+psi[i,j+1,k]+psi[i,j,k]) - \
                               0.25*(psi[i,j,k]+psi[i,j+1,k]+psi[i-1,j+1,k]+psi[i-1,j,k]) )

    def compute_CFL(self, PSI=None):
        """ 
        Compute CFL = max (u*dt/dx)
        """

        # compute U from psi
        self.get_uv(PSI=PSI)

        # compute abs(u*dt/dx)
        self.compute_dudx(PSI=PSI)

        CFL=self._U.max()[1]
        self._U.destroy()
        return CFL

    def compute_dudx(self, PSI=None):
        """
        Compute abs(u*dt/dx)

        """
        # get u
        u = self.da.getVecArray(self._U)
        # get dx
        D = self.da.getVecArray(self.grid.D)

        dt = self.tstepper.dt
        kdxu = self.grid._k_dxu
        (xs, xe), (ys, ye), (zs, ze) = self.da.getRanges()
    
        for k in range(zs,ze):
            u[:,:,k] = u[:,:,k]*dt/D[:,:,kdxu]

    def compute_KE(self, PSI=None):
        """ 
        Compute kinetic energy = 0.5 * sum (u**2+v**2)
        """
        
        # compute local kinetic energy
        self.compute_local_KE(PSI=PSI)
        
        # average spatially
        KE=self._lKE.sum()
        Vol=self._Vol.sum()
        self._lKE.destroy()
        self._Vol.destroy()
        
        return KE/Vol
        
    def compute_local_KE(self, PSI=None):
        
        ### create global vectors
        self._lKE = self.da.createGlobalVec()
        self._Vol = self.da.createGlobalVec()

        ### create local vectors
        local_PSI  = self.da.createLocalVec()
        local_D = self.da.createLocalVec()

        #### load vector PSI used to compute U and V
        if PSI is None:
            self.da.globalToLocal(self.PSI, local_PSI)
        else:
            self.da.globalToLocal(PSI, local_PSI)

        self.da.globalToLocal(self.grid.D, local_D)

        #
        lKE = self.da.getVecArray(self._lKE)
        Vol = self.da.getVecArray(self._Vol)
        psi = self.da.getVecArray(local_PSI)
        D = self.da.getVecArray(local_D)

        mx, my, mz = self.da.getSizes()
        (xs, xe), (ys, ye), (zs, ze) = self.da.getRanges()

        kmask = self.grid._k_mask
        kdxu = self.grid._k_dxu
        kdyu = self.grid._k_dyu
        kdxv = self.grid._k_dxv
        kdyv = self.grid._k_dyv
        kdxt = self.grid._k_dxt
        kdyt = self.grid._k_dyt

        # Loop around volume
        for k in range(zs,ze):
            for j in range(ys, ye):
                for i in range(xs, xe): 
                    if (i==0    or j==0 or
                        i==mx-1 or j==my-1):
                        # lateral boundaries
                        lKE[i, j, k] = 0.
                        Vol[i,j,k] = 0.
                    else:
                        u = - 1. /D[i,j,kdyu] * \
                             ( 0.25*(psi[i+1,j,k]+psi[i+1,j+1,k]+psi[i,j+1,k]+psi[i,j,k]) - \
                               0.25*(psi[i+1,j-1,k]+psi[i+1,j,k]+psi[i,j,k]+psi[i,j-1,k]) )
                        v =   1. /D[i,j,kdxv] * \
                             ( 0.25*(psi[i+1,j,k]+psi[i+1,j+1,k]+psi[i,j+1,k]+psi[i,j,k]) - \
                               0.25*(psi[i,j,k]+psi[i,j+1,k]+psi[i-1,j+1,k]+psi[i-1,j,k]) )
                        Vol[i,j,k] = self.grid.dzt[k]*D[i,j,kdxt]*D[i,j,kdyt]
                        lKE[i,j,k] = 0.5 * (u**2 + v**2) *Vol[i,j,k]

        return
	
    def set_identity(self):
        ONE = self.da.createGlobalVec()
        one = self.da.getVecArray(ONE)
        (xs, xe), (ys, ye), (zs, ze) = self.da.getRanges()
        one[:] = 1.
        return ONE

