#!/usr/bin/python
# -*- encoding: utf8 -*-


import sys
import numpy as np

from petsc4py import PETSc

from .inout import write_nc

#
#==================== Parallel solver ============================================
#

class omegainv():
    """ Omega equation inversion, parallel
    """
    
    def __init__(self, qg, substract_fprime):
        """ Setup the PV inversion solver
        """
                
        self._verbose = qg._verbose
        self._substract_fprime = substract_fprime
        
        # create the operator
        self.L = qg.da.createMat()
        #
        if self._verbose>0:
            print('Operator L declared')

        # Fill in operator values
        if qg.grid._flag_hgrid_uniform and qg.grid._flag_vgrid_uniform:
            self._set_L(self.L, qg)
        else:
            self._set_L_curv(self.L, qg)

        #
        if self._verbose>0:
            print('Operator L filled')

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
            print('PV inversion is set up')
            
            

    def solve(self, qg):
        """ Compute the PV inversion
        """
        # ONE = qg.set_identity()
        # self.L.mult(ONE,self._RHS)
        # write_nc([self._RHS], ['id'], 'data/identity.nc', qg)
        # compute L*PSI and store in self._RHS
        # qg.omega.L.mult(qg.PSI,self._RHS)
        # store L*PSI in netcdf file lpsi.nc
        # write_nc([self._RHS], ['rhs'], 'data/lpsiin.nc', qg)
        # Initialize  RHS
        self.set_rhs(qg)
        # fix boundaries
        self.set_rhs_bdy(qg)
        # mask rhs 
        self.set_rhs_mask(qg)
        # store RHS in netcdf file rhs.nc
        # write_nc([self._RHS], ['rhs'], 'data/rhs.nc', qg)
        # qg.PSI.set(0)
        # actually solves the pb
        self.ksp.solve(self._RHS, qg.PSI)
        # compute L*PSI and store in self._RHS
        # self.L.mult(qg.PSI,self._RHS)
        # store L*PSI in netcdf file lpsi.nc
        # write_nc([self._RHS], ['Lpsi'], 'data/lpsiout.nc', qg)

        if self._verbose>1:
            print('omega equation done')
            

    def set_uv_from_psi(self, qg, PSI=None):
        """
        Compute U & V from Psi
        U = -dPSIdy
        V =  dPSIdx
        """

         ### create global vectors
        self._U = qg.da.createGlobalVec()
        self._V = qg.da.createGlobalVec()

         ### create local vectors
        local_PSI  = qg.da.createLocalVec()
        local_D = qg.da.createLocalVec()

        #### load vector PSI used to compute U and V
        if PSI is None:
            qg.da.globalToLocal(qg.PSI, local_PSI)
        else:
            qg.da.globalToLocal(PSI, local_PSI)

        qg.da.globalToLocal(qg.grid.D, local_D)

        #
        u = qg.da.getVecArray(self._U)
        v = qg.da.getVecArray(self._V)
        psi = qg.da.getVecArray(local_PSI)
        D = qg.da.getVecArray(local_D)


        mx, my, mz = qg.da.getSizes()
        (xs, xe), (ys, ye), (zs, ze) = qg.da.getRanges()

        kmask = qg.grid._k_mask
        kdxu = qg.grid._k_dxu
        kdyu = qg.grid._k_dyu
        kdxv = qg.grid._k_dxv
        kdyv = qg.grid._k_dyv
        kdxt = qg.grid._k_dxt
        kdyt = qg.grid._k_dyt

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

    def set_rho_from_psi(self, qg, PSI=None):
        """
        Compute RHO from Psi
        rho=-rho0*f0/g dPSIdz
        """
     
        ### declare local vectors
        local_PSI  = qg.da.createLocalVec()
        local_D = qg.da.createLocalVec()
        ### declare global vectors
        self._RHO = qg.da.createGlobalVec()
        ###

        #### load vector PSI used to compute RHO
        if PSI is None:
            qg.da.globalToLocal(qg.PSI, local_PSI)
        else:
            qg.da.globalToLocal(PSI, local_PSI)

        qg.da.globalToLocal(qg.grid.D, local_D)
        #
        psi = qg.da.getVecArray(local_PSI)
        rho = qg.da.getVecArray(self._RHO)
        D = qg.da.getVecArray(local_D)

        mx, my, mz = qg.da.getSizes()
        (xs, xe), (ys, ye), (zs, ze) = qg.da.getRanges()

        kmask = qg.grid._k_mask
        kdxu = qg.grid._k_dxu
        kdyu = qg.grid._k_dyu
        kdxv = qg.grid._k_dxv
        kdyv = qg.grid._k_dyv
        kdxt = qg.grid._k_dxt
        kdyt = qg.grid._k_dyt


        # Initialize rho=-rho0*f0/g dpsidz
        for k in range(zs,ze):
            for j in range(ys, ye):
                for i in range(xs, xe): 
                    if (k==0    or k==mz-1):
                        # top and bottom boundaries
                        rho[i, j, k] = 0.
                    else:
                        rho[i,j,k] = - qg.rho0*qg.f0/qg.g/qg.grid.dzt[k] * \
                             ( 0.5*(psi[i,j,k+1]+psi[i,j,k]) - \
                               0.5*(psi[i,j,k]+psi[i,j,k-1]) )
   
    def set_Q(self,qg, U=None, V=None, RHO=None):
        """
        # qxu = g/f0/rho0 * (dudx*drhodx + dvdx*drhody) at u point
        # qyv = g/f0/rho0 * (dudy*drhodx + dvdy*drhody) at v point 
        """

        ### declare local vectors
        local_U = qg.da.createLocalVec()
        local_V = qg.da.createLocalVec()
        local_RHO = qg.da.createLocalVec()
        local_D = qg.da.createLocalVec()
        ### declare global vectors
        self._QXU = qg.da.createGlobalVec()
        self._QYV = qg.da.createGlobalVec()
        ### load vector U,V,RHO used to compute the jacobian
        qg.da.globalToLocal(qg.grid.D, local_D)
        if U is None:
            qg.da.globalToLocal(self._U, local_U)
        else:
            qg.da.globalToLocal(U, local_U)
        if V is None:
            qg.da.globalToLocal(self._V, local_V)
        else:
            qg.da.globalToLocal(V, local_V)
        if RHO is None:
            qg.da.globalToLocal(self._RHO, local_RHO)
        else:
            qg.da.globalToLocal(RHO, local_RHO)
        #
        u = qg.da.getVecArray(local_U)
        v = qg.da.getVecArray(local_V)
        rho = qg.da.getVecArray(local_RHO)
        D = qg.da.getVecArray(local_D)
        qxu = qg.da.getVecArray(self._QXU)
        qyv = qg.da.getVecArray(self._QYV)


        mx, my, mz = qg.da.getSizes()
        (xs, xe), (ys, ye), (zs, ze) = qg.da.getRanges()

        kmask = qg.grid._k_mask
        kdxu = qg.grid._k_dxu
        kdyu = qg.grid._k_dyu
        kdxv = qg.grid._k_dxv
        kdyv = qg.grid._k_dyv
        kdxt = qg.grid._k_dxt
        kdyt = qg.grid._k_dyt

        # qxu = g/f0/rho0 * (dudx*drhodx + dvdx*drhody) at u point
        # qyv = g/f0/rho0 * (dudy*drhodx + dvdy*drhody) at v point 
        for k in range(zs,ze):
            for j in range(ys, ye):
                for i in range(xs, xe): 
                    if (i==0    or j==0 or
                        i==mx-1 or j==my-1 ):
                        # lateral boundaries
                        qxu[i, j, k] = 0.
                        qyv[i, j, k] = 0. 
                    else:
                        qxu[i,j,k] = qg.g/qg.f0/qg.rho0 * ( \
                              (0.5*(u[i+1,j,k]+u[i,j,k])-0.5*(u[i,j,k]+u[i-1,j,k]))/D[i,j,kdxu] * \
                              (rho[i+1,j,k]-rho[i,j,k])/D[i,j,kdxu] + \
                              (0.5*(v[i+1,j,k]+v[i+1,j-1,k])-0.5*(v[i,j,k]+v[i,j-1,k]))/D[i,j,kdxu] * \
                              (0.25*(rho[i+1,j,k]+rho[i+1,j+1,k]+rho[i,j+1,k]+rho[i,j,k]) - \
                               0.25*(rho[i+1,j-1,k]+rho[i+1,j,k]+rho[i,j,k]+rho[i,j-1,k]) )/D[i,j,kdyu] \
                                                          )
                        qyv[i,j,k] = qg.g/qg.f0/qg.rho0 * ( \
                              (0.5*(u[i,j+1,k]+u[i-1,j+1,k])-0.5*(u[i,j,k]+u[i-1,j,k]))/D[i,j,kdyv] * \
                              (0.25*(rho[i+1,j,k]+rho[i+1,j+1,k]+rho[i,j+1,k]+rho[i,j,k]) - \
                               0.25*(rho[i,j,k]+rho[i,j+1,k]+rho[i-1,j+1,k]+rho[i-1,j,k]))/D[i,j,kdxv] + \
                              (0.5*(v[i,j,k]+v[i,j+1,k]) - 0.5*(v[i,j,k]+v[i,j-1,k]))/D[i,j,kdyv] * \
                              (rho[i,j+1,k]-rho[i,j,k])/D[i,j,kdyv] \
                                                      )

    def set_rhs(self, qg):
        """
        Compute the RHS of the omega equation i.e: 2*f0*nabla.Q with Q=-J(nabla.psi,dpsidz)
        """


        # Initialize u=-dpsidy and v=dpsidx
        self.set_uv_from_psi(qg)        
        write_nc([self._U,self._V], ['u','v'], '../output/uv.nc', qg)


        # Initialize rho=-f0/g dpsidz
        self.set_rho_from_psi(qg)
        write_nc([self._RHO], ['rho'], '../output/rho.nc', qg)

        # Initialize jacobian
        self.set_Q(qg)

        ### declare local vectors

        local_QXU = qg.da.createLocalVec()
        local_QYV = qg.da.createLocalVec()
        local_D = qg.da.createLocalVec()
        ###
        qg.da.globalToLocal(self._QXU, local_QXU)
        qg.da.globalToLocal(self._QYV, local_QYV)
        qg.da.globalToLocal(qg.grid.D, local_D)
        #
        qxu = qg.da.getVecArray(local_QXU)
        qyv = qg.da.getVecArray(local_QYV)
        D = qg.da.getVecArray(local_D)

        rhs = qg.da.getVecArray(self._RHS)

        mx, my, mz = qg.da.getSizes()
        (xs, xe), (ys, ye), (zs, ze) = qg.da.getRanges()

        kmask = qg.grid._k_mask
        kdxu = qg.grid._k_dxu
        kdyu = qg.grid._k_dyu
        kdxv = qg.grid._k_dxv
        kdyv = qg.grid._k_dyv
        kdxt = qg.grid._k_dxt
        kdyt = qg.grid._k_dyt

        for k in range(zs,ze):
            for j in range(ys, ye):
                for i in range(xs, xe):
                    if (i==0    or j==0 or
                        i==mx-1 or j==my-1 or
                        k==0    or k==mz-1 ):
                        # lateral boundaries
                        rhs[i, j, k] = 0.
                    else:

                        # qx = qxu moyenné au niveau w, qy = qyv moyenné au niveau w,
                        qxi = 0.5*(qxu[i,j,k+1]+qxu[i,j,k])
                        qxim = 0.5*(qxu[i-1,j,k+1]+qxu[i-1,j,k])  
                        qyj = 0.5*(qyv[i,j,k+1]+qyv[i,j,k])    
                        qyjm = 0.5*(qyv[i,j-1,k+1]+qyv[i,j-1,k])


                        # rhs = 2*f0* nabla.Q with 
                        # nabla:curvilinear operator =1/dx/dy(di(dy.QX) + dj(dx*QY))
                        # Q = vector(QX,QY)

                        rhs[i,j,k] = 2.*qg.f0 / D[i,j,kdxt] / D[i,j,kdyt] * ( \
                                     (D[i,j,kdyu]*qxi - D[i-1,j,kdyu]*qxim) + \
                                     (D[i,j,kdxv]*qyj - D[i,j-1,kdxv]*qyjm) )


    def set_rhs_bdy(self, qg):
        """
        Set South/North, East/West, Bottom/Top boundary conditions
        Set RHS along boundaries for inversion, may be an issue
        for time stepping
        :param da: abstract distributed memory object of the domain
        :param qg: qg_model instance
        :return:
        """
        
        if self._verbose>0:
            print('set RHS along boudaries for inversion ')

        self.set_rhs_bdy_bottom(qg)
        self.set_rhs_bdy_top(qg)
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
            psi = qg.da.getVecArray(qg.PSI)
        else:
            psi = qg.da.getVecArray(PSI)
        if RHO is None:
            rho = qg.da.getVecArray(qg.RHO)
        else:
            rho = qg.da.getVecArray(RHO)
           
        # lower ghost area
        if zs < kdown:
            for k in range(zs,kdown):
                for j in range(ys, ye):
                    for i in range(xs, xe):                    
                        # rhs[i,j,k]=sys.float_info.epsilon
                        rhs[i,j,k]=0.
        # bottom bdy
        k = kdown
        if qg.bdy_type['bottom']=='N' : 
            for j in range(ys, ye):
                for i in range(xs, xe):
                    rhs[i, j, k] = 0.
        elif qg.bdy_type['bottom']=='NBG' : 
            for j in range(ys, ye):
                for i in range(xs, xe):
                    rhs[i, j, k] = 0. 
        elif qg.bdy_type['bottom']=='D':
            for j in range(ys, ye):
                for i in range(xs, xe):
                    rhs[i, j, k] = 0.

                    
        else:
            print(qg.bdy_type['bottom']+' unknown bottom boundary condition')
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
            psi = qg.da.getVecArray(qg.PSI)
        else:
            psi = qg.da.getVecArray(PSI)
        if RHO is None:
            rho = qg.da.getVecArray(qg.RHO)
        else:
            rho = qg.da.getVecArray(RHO)

        if ze > kup+1:
            for k in range(kup+1,ze):
                for j in range(ys, ye):
                    for i in range(xs, xe):
                        # rhs[i,j,k]=sys.float_info.epsilon   
                        rhs[i,j,k]= 0.          
        # upper bdy
        k = kup
        if qg.bdy_type['top']=='N' : 
            for j in range(ys, ye):
                for i in range(xs, xe):
                    rhs[i, j, k] = 0.
        elif qg.bdy_type['top']=='NBG' : 
            for j in range(ys, ye):
                for i in range(xs, xe):
                    rhs[i, j, k] = 0. 
        elif qg.bdy_type['top']=='D' :
            for j in range(ys, ye):
                for i in range(xs, xe):
                    rhs[i, j, k] = 0.

        else:
            print(qg.bdy_type['top']+' unknown top boundary condition')
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
            psi = qg.da.getVecArray(qg.PSI)
        else:
            psi = qg.da.getVecArray(PSI)

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
        # if xs <= istart and qg.BoundaryType is not 'periodic':
            #i = 0
            for k in range(zs, ze):
                for j in range(ys, ye):
                    for i in range(xs,min(xe,istart+1)):
                        rhs[i, j, k] = 0.
        # east bdy
        if xe >= iend:
        # if xe >= iend and qg.BoundaryType is not 'periodic':
            #i = mx - 1
            for k in range(zs, ze):
                for j in range(ys, ye):
                    for i in range(max(xs,iend),xe):
                        rhs[i, j, k] = 0.
                        
        return

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

        psi = qg.da.getVecArray(qg.PSI)

        # interior
        for k in range(zs,ze):
            for j in range(ys, ye):
                for i in range(xs, xe):
                    if mask[i,j,kmask]==0.:
                        rhs[i, j, k] = 0.

        if self._verbose>0:
            print('set RHS mask for inversion ')


    
    
    def _set_L(self,L, qg):
        """ Builds the laplacian operator along with boundary conditions
            Horizontally uniform grid
        """
        
        if qg._verbose>0:
            print('  ... assumes a uniform horizontal and vertical grid')
        
        #
        mx, my, mz = qg.da.getSizes()
        dx, dy, dz = qg.grid.dx, qg.grid.dy, qg.grid.dz
        idx, idy, idz = [1.0/dl for dl in [dx, dy, dz]]
        idx2, idy2, idz2 = [1.0/dl**2 for dl in [dx, dy, dz]]
        #
        (xs, xe), (ys, ye), (zs, ze) = qg.da.getRanges()
        #
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
    
                    # lateral points outside the domain: dirichlet, psi=...
                    if (i<=istart or j<=jstart or
                        i>=iend or j>=jend):
                        L.setValueStencil(row, row, 0.0)
    
                    # bottom bdy condition: default Neuman dpsi/dz=...
                    elif (k==kdown):
                        if qg.bdy_type['bottom']=='N' :
                            L.setValueStencil(row, row, 0.0)
                        elif qg.bdy_type['bottom']=='D':
                            L.setValueStencil(row, row, 0.0)
                        else:
                            print('unknown bottom boundary condition')
                            sys.exit()
    
                    # top bdy condition: default Neuman dpsi/dz=...
                    elif (k==kup):
                        if qg.bdy_type['top']=='N' :
                            L.setValueStencil(row, row, 0.0)
                        elif qg.bdy_type['top']=='D':
                            L.setValueStencil(row, row, 0.0)
                        else:
                            print('unknown top boundary condition')
                            sys.exit()
    
                    # points below and above the domain
                    elif (k<kdown or k>kup):
                        L.setValueStencil(row, row, 0.0)
    
                    # interior points: pv is prescribed
                    else:
                        for index, value in [
                            ((i,j,k-1), qg.f0**2*idz2),
                            ((i,j-1,k), qg.N2[k]*idy2),
                            ((i-1,j,k), qg.N2[k]*idx2),
                            ((i, j, k), -2.*qg.N2[k]*(idx2+idy2)-2.*qg.f0**2*idz2),
                            ((i+1,j,k), qg.N2[k]*idx2),
                            ((i,j+1,k), qg.N2[k]*idy2),
                            ((i,j,k+1), qg.f0**2*idz2)
                            ]:
                            col.index = index
                            col.field = 0
                            L.setValueStencil(row, col, value)
        L.assemble()
        return
    
    
    
    
    
    def _set_L_curv(self,L, qg):
        """ Builds the laplacian operator along with boundary conditions
            Horizontally uniform grid
        """
        
        if qg._verbose>0:
            print('  ... assumes a curvilinear and/or vertically stretched grid')
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
        idzt = 1./qg.grid.dzt
        idzw = 1./qg.grid.dzw
        #idz, idz2 = 1./dz, 1./dz**2
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
                    elif (    (i<=istart and qg.BoundaryType is not 'periodic') \
                           or (i>=iend and qg.BoundaryType is not 'periodic') \
                           or j<=jstart or j>=jend):
                        L.setValueStencil(row, row, 1.0)              
    
                    # bottom bdy condition: default Neuman dpsi/dz=...
                    elif (k==kdown):
                        if qg.bdy_type['bottom']=='N' : 
                            L.setValueStencil(row, row, 1.0)         
                        elif qg.bdy_type['bottom']=='D' :
                            L.setValueStencil(row, row, 1.0)        
                        else:
                            print('unknown bottom boundary condition')
                            sys.exit()
    
                    # top bdy condition: default Neuman dpsi/dz=...
                    elif (k==kup):
                        if qg.bdy_type['top']=='N' : 
                            L.setValueStencil(row, row, 1.0)        
                        elif qg.bdy_type['top']=='D':
                            L.setValueStencil(row, row, 1.0)      
                        else:
                            print('unknown top boundary condition')
                            sys.exit()
    
                    # points below and above the domain
                    elif (k<kdown) or (k>kup):
                        L.setValueStencil(row, row, 1.0)       
    
                    # lateral points outside the domain: dirichlet, psi=...
                    # elif (i<=istart or j<=jstart or
                    #     i>=iend or j>=jend):
                    #     L.setValueStencil(row, row, 1.0)
    
                    # interior points: pv is prescribed
                    else:
                        
                        for index, value in [

                            ((i,j,k-1), qg.f0**2*idzt[k]*idzw[k]),
                            ((i,j-1,k), qg.N2[k]/D[i,j,kdxt]/D[i,j,kdyt] * D[i,j-1,kdxv]/D[i,j-1,kdyv]),
                            ((i-1,j,k), qg.N2[k]/D[i,j,kdxt]/D[i,j,kdyt] * D[i-1,j,kdyu]/D[i-1,j,kdxu]),
                            ((i, j, k), -qg.N2[k]/D[i,j,kdxt]/D[i,j,kdyt]*( \
                                             D[i,j,kdyu]/D[i,j,kdxu] \
                                            +D[i-1,j,kdyu]/D[i-1,j,kdxu] \
                                            +D[i,j,kdxv]/D[i,j,kdyv] \
                                            +D[i,j-1,kdxv]/D[i,j-1,kdyv])
                             - (qg.f0**2*idzt[k+1]*idzw[k]+qg.f0**2*idzt[k]*idzw[k])),
                            ((i+1,j,k), qg.N2[k]/D[i,j,kdxt]/D[i,j,kdyt] * D[i,j,kdyu]/D[i,j,kdxu]),
                            ((i,j+1,k), qg.N2[k]/D[i,j,kdxt]/D[i,j,kdyt] * D[i,j,kdxv]/D[i,j,kdyv]),
                            ((i,j,k+1), qg.f0**2*idzt[k+1]*idzw[k])   
                            
                            ]:
                            col.index = index
                            col.field = 0
                            L.setValueStencil(row, col, value)
    
                   
        L.assemble()
        return
    
    
