#!/usr/bin/python
# -*- encoding: utf8 -*-


import sys
import numpy as np

from petsc4py import PETSc

from .inout import write_nc

#
#==================== Parallel solver ============================================
#

class pvinversion():
    """ PV inversion, parallel
    """
    
    def __init__(self, da, grid, bdy_type, verbose=0, sparam=None, substract_fprime=False):
        """ Setup the PV inversion solver
        """

        self._verbose = verbose
        #
        self.bdy_type = bdy_type
        self.petscBoundaryType=False
        if ('periodic' in self.bdy_type) and (self.bdy_type['periodic']):
            self.petscBoundaryType = True
        #
        self._substract_fprime = substract_fprime
        
        # create the operator
        self.L = da.createMat()
        #
        if self._verbose>0:
            print('A PV inversion object is being created')
            print('  Operator L declared')

        # Fill in operator values
        if grid._flag_hgrid_uniform and grid._flag_vgrid_uniform:
            self._set_L(self.L, da, grid, sparam=sparam)
        else:
            self._set_L_curv(self.L, da, grid, sparam=sparam)

        #
        if self._verbose>0:
            print('  Operator L filled')

        # global vector for PV inversion
        self._RHS = da.createGlobalVec()

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
        self.ksp.setInitialGuessNonzero(True)
        # and incomplete Cholesky for preconditionning
        # self.ksp.getPC().setType('icc')
        # self.ksp.getPC().setType('bjacobi')
        # self.ksp.getPC().setType('asm')
        # self.ksp.getPC().setType('mg')
        # self.ksp.getPC().setType('none')
        # self.ksp.setNormType(2)
        # set tolerances
        self.ksp.setTolerances(rtol=1e-4)
        self.ksp.setTolerances(max_it=1000)
        #
        #
        for opt in sys.argv[1:]:
	        PETSc.Options().setValue(opt, None)
        self.ksp.setFromOptions()
        
        if self._verbose>0:
            print('  PV inversion is set up')

#
# ==================== perform inversion ===================================
#


    def solve(self, qg, Q=None, PSI=None, RHO=None):
        """ Compute the PV inversion
        """
        if Q is None or PSI is None:
            print('!Error: pvinv.solve requires Q and PSI')
            sys.exit()
        # ONE = qg.set_identity()
        # qg.pvinv.L.mult(ONE,self._RHS)
        # write_nc([self._RHS], ['id'], 'data/identity.nc', qg)
        # compute L*PSI and store in self._RHS
        #qg.pvinv.L.mult(qg.PSI,self._RHS)
        # store L*PSI in netcdf file lpsi.nc
        #write_nc([self._RHS], ['rhs'], 'output/lpsiin.nc', qg)
        # copy Q into RHS
        Q.copy(self._RHS)
        if self._substract_fprime:
            # substract f-f0 from PV
            self.substract_fprime_from_rhs(qg)
        # fix boundaries
        self.set_rhs_bdy(qg, PSI=PSI, RHO=RHO)
        # apply mask
        if qg.grid.mask:
            # mask rhs
            self.set_rhs_mask(qg)
        # store RHS in netcdf file rhs.nc
        #write_nc([self._RHS], ['rhs'], 'output/rhs.nc', qg)
        # qg.PSI.set(0)
        # actually solves the pb
        self.ksp.solve(self._RHS, qg.state.PSI)
        # compute L*PSI and store in self._RHS
        #qg.pvinv.L.mult(qg.PSI,self._RHS)
        # store L*PSI in netcdf file lpsi.nc
        #write_nc([self._RHS], ['Lpsi'], 'output/lpsiout.nc', qg)

        if self._verbose>1:
            print('Inversion done (%i iterations)' %(self.ksp.getIterationNumber()))
            
            
            
    def substract_fprime_from_rhs(self, qg):
        """
        Substract f'=f-f0 from the rhs used in PV inversion
        """
        rhs = qg.da.getVecArray(self._RHS)
        mx, my, mz = qg.da.getSizes()
        (xs, xe), (ys, ye), (zs, ze) = qg.da.getRanges()

        D = qg.da.getVecArray(qg.grid.D)

        for k in range(zs,ze):
            for j in range(ys, ye):
                for i in range(xs, xe):                    
                    rhs[i,j,k] -= D[i,j,qg.grid._k_f]  - qg.state.f0
        
        if self._verbose>1:
            print('  Substract f-f0 from pv prior to inversion')

        

    def set_rhs_bdy(self, qg, PSI=None, RHO=None):
        """
        Set South/North, East/West, Bottom/Top boundary conditions
        Set RHS along boundaries for inversion, may be an issue
        for time stepping

        :param da: abstract distributed memory object of the domain
        :param qg: qg_model instance
        :return:
        """
        
        if self._verbose>1:
            print('  Set RHS along boudaries for inversion ')

        self.set_rhs_bdy_bottom(qg, PSI=PSI, RHO=RHO)
        self.set_rhs_bdy_top(qg, PSI=PSI, RHO=RHO)
        self.set_rhs_bdy_lat(qg)

        return


    def set_rhs_bdy_bottom(self, qg, PSI=None, RHO=None):
        """
        Set bottom boundary condition

        :param PSI, RHO: Petsc vectors that will be used to compute the bdy condition
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

        # load vector used to compute boundary conditions
        if PSI is None:
            psi = qg.da.getVecArray(qg.state.PSI)
        else:
            psi = qg.da.getVecArray(PSI)
        if RHO is None:
            rho = qg.da.getVecArray(qg.state.RHO)
        else:
            rho = qg.da.getVecArray(RHO)
           
        # lower ghost area
        if zs < kdown:
            for k in range(zs,kdown):
                for j in range(ys, ye):
                    for i in range(xs, xe):                    
                        # rhs[i,j,k]=sys.float_info.epsilon
                        rhs[i,j,k]=psi[i, j, k]
        # bottom bdy
        k = kdown
        if self.bdy_type['bottom']=='N' :
            for j in range(ys, ye):
                for i in range(xs, xe):
                    rhs[i, j, k] = - qg.g*0.5*(rho[i, j, k]+rho[i, j, k+1])/(qg.rho0*qg.state.f0)
        elif self.bdy_type['bottom']=='NBG' :
            for j in range(ys, ye):
                for i in range(xs, xe):
                    rhs[i, j, k] = (psi[i,j,k+1]-psi[i,j,k])/qg.grid.dzw[k] 
        elif self.bdy_type['bottom']=='D':
            for j in range(ys, ye):
                for i in range(xs, xe):
                    rhs[i, j, k] = psi[i,j,k]

                    
        else:
            print(self.bdy_type['bottom']+" unknown bottom boundary condition")
            sys.exit()
                
        return

            
    def set_rhs_bdy_top(self, qg, PSI=None, RHO=None):
        """
        Set top boundary condition

        :param PSI, RHO: Petsc vectors that will be used to compute the bdy condition
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

        # load vector used to compute boundary conditions
        if PSI is None:
            psi = qg.da.getVecArray(qg.state.PSI)
        else:
            psi = qg.da.getVecArray(PSI)
        if RHO is None:
            rho = qg.da.getVecArray(qg.state.RHO)
        else:
            rho = qg.da.getVecArray(RHO)

        if ze > kup+1:
            for k in range(kup+1,ze):
                for j in range(ys, ye):
                    for i in range(xs, xe):
                        # rhs[i,j,k]=sys.float_info.epsilon   
                        rhs[i,j,k]= psi[i,j,k]            
        # upper bdy
        k = kup
        if self.bdy_type['top']=='N' :
            for j in range(ys, ye):
                for i in range(xs, xe):
                    rhs[i, j, k] = - qg.g*0.5*(rho[i, j, k]+rho[i, j, k-1])/(qg.rho0*qg.state.f0)
        elif self.bdy_type['top']=='NBG' :
            for j in range(ys, ye):
                for i in range(xs, xe):
                    rhs[i, j, k] = (psi[i,j,k+1]-psi[i,j,k])/qg.grid.dzw[k] 
        elif self.bdy_type['top']=='D' :
            for j in range(ys, ye):
                for i in range(xs, xe):
                    rhs[i, j, k] = psi[i,j,k]
        else:
            print(self.bdy_type['top']+" unknown top boundary condition")
            sys.exit()


    def set_rhs_bdy_lat(self, qg, PSI=None):
        """
        Set lateral boundary condition

        :param PSI: Petsc vector that will be used to compute the bdy condition
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

        # load vector used to compute boundary conditions
        if PSI is None:
            psi = qg.da.getVecArray(qg.state.PSI)
        else:
            psi = qg.da.getVecArray(PSI)

        # south bdy
        if ys <= jstart:
            #j = 0
            for k in range(zs, ze):
                for j in range(ys,min(ye,jstart+1)):
                    for i in range(xs, xe):
                        rhs[i, j, k] = psi[i, j, k]
        # north bdy
        if ye >= jend:
            #j = my - 1
            for k in range(zs, ze):
                for j in range(max(ys,jend),ye):
                    for i in range(xs, xe):
                        rhs[i, j, k] = psi[i, j, k]
        # west bdy
        if xs <= istart and self.petscBoundaryType is not 'periodic':
            #i = 0
            for k in range(zs, ze):
                for j in range(ys, ye):
                    for i in range(xs,min(xe,istart+1)):
                        rhs[i, j, k] = psi[i, j, k]
        # east bdy
        if xe >= iend and self.petscBoundaryType is not 'periodic':
            #i = mx - 1
            for k in range(zs, ze):
                for j in range(ys, ye):
                    for i in range(max(xs,iend),xe):
                        rhs[i, j, k] = psi[i, j, k]
                        
        return
    


    def set_rhs_mask(self, qg):
        """
        Set mask on rhs: where mask=0 (land) rhs=psi

        param: da abstract distributed memory object of the domain
        param: qg qg_model instance
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

        psi = qg.da.getVecArray(qg.state.PSI)
        
        # interior
        for k in range(zs,ze):
            for j in range(ys, ye):
                for i in range(xs, xe):
                    if mask[i,j,kmask]==0.:
                        rhs[i, j, k] = psi[i,j,k]

        if self._verbose>1:
            print('  Set RHS mask for inversion ')

#
# ==================== Define elliptical operators ===================================
#
    
    def _set_L(self,L, da, grid, sparam):
        """ Builds the laplacian operator along with boundary conditions
            Horizontally uniform grid
        """
        
        if self._verbose>0:
            print('  ... assumes a uniform horizontal and vertical grid')
        
        #
        mx, my, mz = da.getSizes()
        dx, dy, dz = grid.dx, grid.dy, grid.dz
        idx, idy, idz = [1.0/dl for dl in [dx, dy, dz]]
        idx2, idy2, idz2 = [1.0/dl**2 for dl in [dx, dy, dz]]
        #
        (xs, xe), (ys, ye), (zs, ze) = da.getRanges()
        #
        istart = grid.istart
        iend = grid.iend
        jstart = grid.jstart
        jend = grid.jend
        kdown = grid.kdown
        kup = grid.kup
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
    
                    # lateral points outside the domain: dirichlet, psi=...
                    if (i<=istart or j<=jstart or
                        i>=iend or j>=jend):
                        L.setValueStencil(row, row, 1.0)
    
                    # bottom bdy condition: default Neuman dpsi/dz=...
                    elif (k==kdown):
                        if self.bdy_type['bottom']=='N' or self.bdy_type['bottom']=='NBG':
                            for index, value in [
                                ((i,j,k), -idz),
                                ((i,j,k+1),  idz)
                                ]:
                                col.index = index
                                col.field = 0
                                L.setValueStencil(row, col, value)
                        elif self.bdy_type['bottom']=='D':
                            L.setValueStencil(row, row, 1.0)
                        else:
                            print('unknown bottom boundary condition')
                            sys.exit()
    
                    # top bdy condition: default Neuman dpsi/dz=...
                    elif (k==kup):
                        if self.bdy_type['top']=='N' or self.bdy_type['top']=='NBG':
                            for index, value in [
                                ((i,j,k-1), -idz),
                                ((i,j,k),  idz),
                                ]:
                                col.index = index
                                col.field = 0
                                L.setValueStencil(row, col, value)
                        elif self.bdy_type['top']=='D':
                            L.setValueStencil(row, row, 1.0)
                        else:
                            print('unknown top boundary condition')
                            sys.exit()
    
                    # points below and above the domain
                    elif (k<kdown or k>kup):
                        L.setValueStencil(row, row, 0.0)
    
                    # interior points: pv is prescribed
                    else:
                        for index, value in [
                            ((i,j,k-1), sparam[k-1]*idz2),
                            ((i,j-1,k), idy2),
                            ((i-1,j,k), idx2),
                            ((i, j, k), -2.*(idx2+idy2)-(sparam[k]*idz2+sparam[k-1]*idz2)),
                            ((i+1,j,k), idx2),
                            ((i,j+1,k), idy2),
                            ((i,j,k+1), sparam[k]*idz2)
                            ]:
                            col.index = index
                            col.field = 0
                            L.setValueStencil(row, col, value)
        L.assemble()
        return
    
    
    
    
    
    def _set_L_curv(self,L, da, grid):
        """ Builds the laplacian operator along with boundary conditions
            Horizontally uniform grid
        """
        
        if self._verbose>0:
            print('  ... assumes a curvilinear and/or vertically stretched grid')
        #
        mx, my, mz = da.getSizes()
        #
        #D = qg.da.getVecArray(qg.grid.D)
        local_D  = da.createLocalVec()
        da.globalToLocal(grid.D, local_D)
        D = da.getVecArray(local_D)
        kmask = grid._k_mask
        kdxu = grid._k_dxu
        kdyu = grid._k_dyu
        kdxv = grid._k_dxv
        kdyv = grid._k_dyv
        kdxt = grid._k_dxt
        kdyt = grid._k_dyt
        #
        idzt = 1./grid.dzt
        idzw = 1./grid.dzw
        #idz, idz2 = 1./dz, 1./dz**2
        #
        (xs, xe), (ys, ye), (zs, ze) = da.getRanges()
        istart = grid.istart
        iend = grid.iend
        jstart = grid.jstart
        jend = grid.jend
        kdown = grid.kdown
        kup = grid.kup
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
                    elif (    (i<=istart and self.BoundaryType is not 'periodic') \
                           or (i>=iend and self.BoundaryType is not 'periodic') \
                           or j<=jstart or j>=jend):
                        L.setValueStencil(row, row, 1.0)
    
                    # bottom bdy condition: default Neuman dpsi/dz=...
                    elif (k==kdown):
                        if self.bdy_type['bottom']=='N' or self.bdy_type['bottom']=='NBG':
                            for index, value in [
                                ((i,j,k), -idzw[k]),
                                ((i,j,k+1),  idzw[k])
                                ]:
                                col.index = index
                                col.field = 0
                                L.setValueStencil(row, col, value)
                        elif self.bdy_type['bottom']=='D' :
                            L.setValueStencil(row, row, 1.0)
                        else:
                            print('unknown bottom boundary condition')
                            sys.exit()
    
                    # top bdy condition: default Neuman dpsi/dz=...
                    elif (k==kup):
                        if self.bdy_type['top']=='N' or self.bdy_type['top']=='NBG':
                            for index, value in [
                                ((i,j,k-1), -idzw[k-1]),
                                ((i,j,k),  idzw[k-1]),
                                ]:
                                col.index = index
                                col.field = 0
                                L.setValueStencil(row, col, value)
                        elif self.bdy_type['top']=='D':
                            L.setValueStencil(row, row, 1.0)
                        else:
                            print('unknown top boundary condition')
                            sys.exit()
    
                    # points below and above the domain
                    elif (k<kdown) or (k>kup):
                        # L.setValueStencil(row, row, 0.0)
                        L.setValueStencil(row, row, 1.0)
    
                    # lateral points outside the domain: dirichlet, psi=...
                    # elif (i<=istart or j<=jstart or
                    #     i>=iend or j>=jend):
                    #     L.setValueStencil(row, row, 1.0)
    
                    # interior points: pv is prescribed
                    else:
                        
                        for index, value in [
    
                            ((i,j,k-1), qg._sparam[k-1]*idzt[k]*idzw[k-1]),
                            ((i,j-1,k), 1./D[i,j,kdxt]/D[i,j,kdyt] * D[i,j-1,kdxv]/D[i,j-1,kdyv]),
                            ((i-1,j,k), 1./D[i,j,kdxt]/D[i,j,kdyt] * D[i-1,j,kdyu]/D[i-1,j,kdxu]),
                            ((i, j, k), -1./D[i,j,kdxt]/D[i,j,kdyt]*( \
                                             D[i,j,kdyu]/D[i,j,kdxu] \
                                            +D[i-1,j,kdyu]/D[i-1,j,kdxu] \
                                            +D[i,j,kdxv]/D[i,j,kdyv] \
                                            +D[i,j-1,kdxv]/D[i,j-1,kdyv])
                             - (qg._sparam[k]*idzt[k]*idzw[k]+qg._sparam[k-1]*idzt[k]*idzw[k-1])),
                            ((i+1,j,k), 1./D[i,j,kdxt]/D[i,j,kdyt] * D[i,j,kdyu]/D[i,j,kdxu]),
                            ((i,j+1,k), 1./D[i,j,kdxt]/D[i,j,kdyt] * D[i,j,kdxv]/D[i,j,kdyv]),
                            ((i,j,k+1), qg._sparam[k]*idzt[k]*idzw[k])
                            ]:
                            col.index = index
                            col.field = 0
                            L.setValueStencil(row, col, value)
    
                   
        L.assemble()
        return
    
    
