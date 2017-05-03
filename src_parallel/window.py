#!/usr/bin/python
# -*- encoding: utf8 -*-

import sys

from .grid import *

import petsc4py
#from Cython.Compiler.Main import verbose
petsc4py.init(sys.argv)
from petsc4py import PETSc

import numpy as np
from .io import read_nc_petsc, read_nc_petsc_2D
from .io import write_nc




class window():
    """  Computes a window for spectral computations
    """
    
    def __init__(self,
                 hgrid = None, vgrid=None,
                 K = 1.e-6,
                 vdom={}, hdom={},
                 ncores_x=None, ncores_y=None,
                 bdy_type_in={},
                 verbose = 1,
                 ):
        """ Window model creation
        Parameters:
        """

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
            print 'A window model object is being created'
            # print out grid parameters
            print self.grid
            # # print if a subdomain is considered
            # if self.kdown==0 or self.kup<self.grid.Nz-1:
            #     print 'Vertical subdomain: kdown=%d, kup=%d' %(self.kdown, self.kup)
            # if self.istart==0 or self.iend<self.grid.Nx-1 or self.jstart==0 or self.jend<self.grid.Ny-1:
            #     print 'Horizontal subdomain: (istart, iend) = (%d, %d), (jstart, jend) = (%d, %d)' \
            #              %(self.istart, self.iend, self.jstart, self.jend)

                
        #
        self.K = K
        self._K2 = K**2
        #
        # declare petsc vectors
        #
        # PV
        self.Q = self.da.createGlobalVec()
        # streamfunction
        self.PSI = self.da.createGlobalVec()
        # density
        #self.RHO = self.da.createGlobalVec()

        # default top and bottom boudary condition = 'N' pour Neumann. 
        # Other possibility 'D' for Direchlet
        #self.bdy_type = {'top':'D','bottom':'D'}
        #self.bdy_type.update(bdy_type_in)

        # initiate pv inversion solver
        self.wininv = wininversion(self)


    #def set_psi(self, analytical_psi=True, file_psi=None):
    #    """
    #    Set psi to a given value
    #    """
    #    if file_psi is not None:
    #        if self._verbose:
    #            print 'Set psi from file '+file_psi+' ...'
    #        read_nc_petsc(self.PSI, 'psi', file_psi, self, fillmask=0.)
    #    elif analytical_psi:
    #        self.set_psi_analytically()

    #def set_psi_analytically(self):
    #    """ Set psi analytically
    #    """
    #    psi = self.da.getVecArray(self.PSI)
    #    mx, my, mz = self.da.getSizes()
    #    (xs, xe), (ys, ye), (zs, ze) = self.da.getRanges()
    #    #
    #    if self._verbose:
    #        print 'Set psi analytically'
    #    for k in range(zs, ze):
    #        for j in range(ys, ye):
    #            for i in range(xs, xe):
    #                psi[i, j, k] = 1.


    def set_q(self, analytical_q=True, file_q=None):
        """ Set q to a given value
        """
        #
        if file_q is not None:
            if self._verbose:
                print 'Set q from file '+file_q+' ...'
            read_nc_petsc(self.Q, 'q', file_q, self, fillmask=0.)
        elif analytical_q:
            self.set_q_analytically()

    def set_q_analytically(self):
        """ Set q analytically
        """
        q = self.da.getVecArray(self.Q)
        mx, my, mz = self.da.getSizes()
        (xs, xe), (ys, ye), (zs, ze) = self.da.getRanges()
        #
        if self._verbose:
            print 'Set q analytically'
        for k in range(zs, ze):
            for j in range(ys, ye):
                for i in range(xs, xe):
                    q[i, j, k] = -self._K2


    def invert_win(self):
        """ wrapper around solver solve method
        """
        self.wininv.solve(self)





#
#==================== Parallel solver ============================================
#

class wininversion():
    """ Window inversion, parallel
    """
    
    def __init__(self, qg):
        """ Setup the PV inversion solver
        """
                
        self._verbose = qg._verbose
        
        # create the operator
        self.L = qg.da.createMat()
        #
        if self._verbose>0:
            print 'Operator L declared'

        # Fill in operator values
        self._set_L_curv(self.L, qg)

        #
        if self._verbose>0:
            print 'Operator L filled'

        # global vector for PV inversion
        self._RHS = qg.da.createGlobalVec()

        # local vectors
        #self._localRHS  = qg.da.createLocalVec()
        #self._localPSI  = qg.da.createLocalVec()

        # create solver
        self.ksp = PETSc.KSP()
        self.ksp.create(PETSc.COMM_WORLD)
        self.ksp.setOperators(self.L)
        # self.ksp.setType('cg')
        self.ksp.setType('gmres')
        # self.ksp.setType('bicg')
        self.ksp.setInitialGuessNonzero(False)
        # and incomplete Cholesky for preconditionning
        # self.ksp.getPC().setType('icc')
        # self.ksp.getPC().setType('bjacobi')
        # self.ksp.getPC().setType('asm')
        # self.ksp.getPC().setType('mg')
        # self.ksp.getPC().setType('none')
        # self.ksp.setNormType(2)
        # set tolerances
        self.ksp.setTolerances(rtol=1e-7)
        self.ksp.setTolerances(max_it=1000)
        #
        #
        for opt in sys.argv[1:]:
            PETSc.Options().setValue(opt, None)
        self.ksp.setFromOptions()
        
        if self._verbose>0:
            print 'PV inversion is set up'
            
            

    def solve(self, qg):
        """ Compute the PV inversion
        """
        # copy Q into RHS
        qg.Q.copy(self._RHS)
        # fix boundaries
        self.set_rhs_bdy(qg)
        # mask rhs 
        self.set_rhs_mask(qg)
        # actually solves the pb
        self.ksp.solve(self._RHS, qg.PSI)

        if self._verbose>1:
            print 'Inversion done'

        

    def set_rhs_bdy(self, qg):
        """
        Set South/North, East/West, Bottom/Top boundary conditions
        Set RHS along boundaries for inversion, may be an issue
        for time stepping
        :param da: abstract distributed memory object of the domain
        :param qg: qg_model instance
        :return:
        """

        rhs = qg.da.getVecArray(self._RHS)
        mx, my, mz = qg.da.getSizes()
        (xs, xe), (ys, ye), (zs, ze) = qg.da.getRanges()

        istart = qg.grid.istart
        iend = qg.grid.iend
        jstart = qg.grid.jstart
        jend = qg.grid.jend
        kdown = qg.grid.kdown
        kup = qg.grid.kup

        # lower ghost area
        if zs < kdown:
            for k in range(zs,kdown):
                for j in range(ys, ye):
                    for i in range(xs, xe):                    
                        # rhs[i,j,k]=sys.float_info.epsilon
                        rhs[i,j,k]=0.                   

        if ze > kup+1:
            for k in range(kup+1,ze):
                for j in range(ys, ye):
                    for i in range(xs, xe):
                        # rhs[i,j,k]=sys.float_info.epsilon   
                        rhs[i,j,k]= 0.

        # south bdy
        if ys <= jstart:
            #j = 0
            for k in range(zs, ze):
                for j in range(ys,min(ye,jstart+1)):
                    for i in range(xs, xe):
                        rhs[i, j, k] = 0.
        # north bdy
        if ye >= jend:
            #j = my - 1
            for k in range(zs, ze):
                for j in range(max(ys,jend),ye):
                    for i in range(xs, xe):
                        rhs[i, j, k] = 0.
        # west bdy
        if xs <= istart:
            #i = 0
            for k in range(zs, ze):
                for j in range(ys, ye):
                    for i in range(xs,min(xe,istart+1)):
                        rhs[i, j, k] = 0.
        # east bdy
        if xe >= iend:
            #i = mx - 1
            for k in range(zs, ze):
                for j in range(ys, ye):
                    for i in range(max(xs,iend),xe):
                        rhs[i, j, k] = 0.

        if self._verbose>0:
            print 'set RHS along boudaries for inversion '


    def set_rhs_mask(self, qg):
        """
        Set mask on rhs: where mask=0 (land) rhs=psi
        - param da: abstract distributed memory object of the domain
        - param qg: qg_model instance
             qg.grid.D[qg.grid._k_mask]: mask
        - self.rhs : vorticity whith boundary conditions
        return: masked rhs
        """

        rhs = qg.da.getVecArray(self._RHS)
        mask = qg.da.getVecArray(qg.grid.D)
        mx, my, mz = qg.da.getSizes()
        (xs, xe), (ys, ye), (zs, ze) = qg.da.getRanges()

        istart = qg.grid.istart
        iend = qg.grid.iend
        jstart = qg.grid.jstart
        jend = qg.grid.jend
        kdown = qg.grid.kdown
        kup = qg.grid.kup
        kmask = qg.grid._k_mask

        # interior
        for k in range(zs,ze):
            for j in range(ys, ye):
                for i in range(xs, xe):
                    if mask[i,j,kmask]==0.:
                        rhs[i, j, k] = 0.

        if self._verbose>0:
            print 'set RHS mask for inversion '


    
    def _set_L_curv(self,L, qg):
        """ Builds the laplacian operator along with boundary conditions
            Horizontally uniform grid
        """
        
        if qg._verbose>0:
            print '  ... assumes a curvilinear and/or vertically stretched grid'
        #
        mx, my, mz = qg.da.getSizes()
        #
        #D = qg.da.getVecArray(qg.grid.D)
        local_D  = qg.da.createLocalVec()
        qg.da.globalToLocal(qg.grid.D, local_D)
        D = qg.da.getVecArray(local_D)
        kmask = qg.grid._k_mask
        kdxu = qg.grid._k_dxu
        kdyu = qg.grid._k_dyu
        kdxv = qg.grid._k_dxv
        kdyv = qg.grid._k_dyv
        kdxt = qg.grid._k_dxt
        kdyt = qg.grid._k_dyt
        #

        #
        (xs, xe), (ys, ye), (zs, ze) = qg.da.getRanges()
        istart = qg.grid.istart
        iend = qg.grid.iend
        jstart = qg.grid.jstart
        jend = qg.grid.jend
        kdown = qg.grid.kdown
        kup = qg.grid.kup
        #
        L.zeroEntries()
        row = PETSc.Mat.Stencil()
        col = PETSc.Mat.Stencil()
        #
        for k in range(zs, ze):
            for j in range(ys, ye):
                for i in range(xs, xe):
                    row.index = (i,j,k)
                    row.field = 0
    
                    # masked points (land=0), L=1
                    if D[i,j,kmask]==0.:
                        L.setValueStencil(row, row, 1.)
    
                    # lateral points outside the domain: dirichlet, psi=...
                    elif (i<=istart or j<=jstart or
                        i>=iend or j>=jend):
                        L.setValueStencil(row, row, 1.0)
    
                    # points below and above the domain
                    elif (k<kdown) or (k>kup):
                        # L.setValueStencil(row, row, 0.0)
                        L.setValueStencil(row, row, 1.0)
    
                    # interior points: pv is prescribed
                    else:
                        
                        for index, value in [
                            ((i,j-1,k), 1./D[i,j,kdxt]/D[i,j,kdyt] * D[i,j-1,kdxv]/D[i,j-1,kdyv]),
                            ((i-1,j,k), 1./D[i,j,kdxt]/D[i,j,kdyt] * D[i-1,j,kdyu]/D[i-1,j,kdxu]),
                            ((i, j, k), -qg._K2 -1./D[i,j,kdxt]/D[i,j,kdyt]*( \
                                             D[i,j,kdyu]/D[i,j,kdxu] \
                                            +D[i-1,j,kdyu]/D[i-1,j,kdxu] \
                                            +D[i,j,kdxv]/D[i,j,kdyv] \
                                            +D[i,j-1,kdxv]/D[i,j-1,kdyv])),
                            ((i+1,j,k), 1./D[i,j,kdxt]/D[i,j,kdyt] * D[i,j,kdyu]/D[i,j,kdxu]),
                            ((i,j+1,k), 1./D[i,j,kdxt]/D[i,j,kdyt] * D[i,j,kdxv]/D[i,j,kdyv])
                            ]:
                            col.index = index
                            col.field = 0
                            L.setValueStencil(row, col, value)
    
                   
        L.assemble()
        return
    
    



