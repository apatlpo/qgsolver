#!/usr/bin/python
# -*- encoding: utf8 -*-


import sys
from petsc4py import PETSc
from .inout import write_nc

#
#==================== Parallel solver ============================================
#

class omegainv():
    """ Omega equation solver
    """

    def __init__(self, da, grid, bdy_type, f0, N2, verbose=0, solver='gmres', pc=None):
        """ Setup the Omega equation solver

        Parameters
        ----------
        da : petsc DMDA
            holds the petsc grid
        grid : qgsolver grid object
            grid data holder
        bdy_type : dict
            prescribe vertical and lateral boundary conditions. Examples
                bdy_type = {'bottom': 'D', 'top': 'D'}    for Dirichlet bdy conditions
                bdy_type = {'bottom': 'N', 'top': 'N'}    for Neumann bdy conditions
                bdy_type = {'bottom': 'N', 'top': 'N'}    for Neumann bdy conditions using PSI instead of RHO
                bdy_type = {'periodic': None}             for horizontal periodicity
        f0 : float
            averaged Coriolis frequency, used in operator
        N2 : ndarray
            buoyancy frequency, used in operator
        verbose : int
            degree of verbosity, 0 means no outputs
        solver : str
            petsc solver: 'gmres' (default), 'bicg', 'cg'
        pc : str, optional
            what is default?
            preconditionner: 'icc', 'bjacobi', 'asm', 'mg', 'none'

        """

        self._verbose = verbose
        #
        self.bdy_type = bdy_type
        if ('periodic' in self.bdy_type) and (self.bdy_type['periodic']):
            self.petscBoundaryType = True
        else:
            self.petscBoundaryType = False
        # physical parameters
        self.f0 = f0
        self.N2 = N2

        # create the operator
        self.L = da.createMat()
        #
        if self._verbose>0:
            print('An Omega equation inversion object is being created')
            print('  Operator L declared')

        # Fill in operator values
        if grid._flag_hgrid_uniform and grid._flag_vgrid_uniform:
            self._set_L(self.L, da, grid)
        else:
            self._set_L_curv(self.L, da, grid)

        #
        if self._verbose>0:
            print('  Operator L filled')

        # global vector for Omega equation inversion
        self._RHS = da.createGlobalVec()

        # create solver
        self.ksp = PETSc.KSP()
        self.ksp.create(PETSc.COMM_WORLD)
        self.ksp.setOperators(self.L)
        self.ksp.setType(solver)
        self.ksp.setInitialGuessNonzero(False)
        # and incomplete Cholesky for preconditionning
        if pc is not None:
            self.ksp.getPC().setType('icc')
        # set tolerances
        self.ksp.setTolerances(rtol=1e-7)
        self.ksp.setTolerances(max_it=1000)
        #
        #
        for opt in sys.argv[1:]:
            PETSc.Options().setValue(opt, None)
        self.ksp.setFromOptions()
        
        if self._verbose>0:
            print('  Omega equation inversion is set up')

#
# ==================== perform inversion ===================================
#

    def solve(self, da, grid, state, PSI=None, U=None, V=None, RHO=None):
        """ Compute the omega equation inversion
        Uses prioritarily optional Q, PSI, RHO for RHS and bdy conditions

        Parameters
        ----------
        da : petsc DMDA
            holds the petsc grid
        grid : qgsolver grid object
            grid data holder
        state : state object
            ocean state
        PSI : petsc Vec, None, optional
            streamfunction, use state.PSI if None

        Returns
        -------
        state.W:
            Put PV inversion result in state.PSI
        """

        # Initialize  RHS
        self.set_rhs(da, grid, state, PSI, U, V, RHO)

        # actually solves the pb
        self.ksp.solve(self._RHS, state.W)

        # should destroy: self._U, self._V, self._RHO

        if self._verbose>0:
            print('Omega equation solved (%i iterations)' %(self.ksp.getIterationNumber()))

#
# ==================== utils methods for inversions ===================================
#

    def set_rhs(self, da, grid, state, PSI, U, V, RHO):
        """Compute the RHS of the omega equation i.e: 2*f0*nabla.Q with Q=-J(nabla.psi,dpsidz)

        Parameters
        ----------
        da : petsc DMDA
            holds the petsc grid
        grid : qgsolver grid object
            grid data holder
        state : state object
            ocean state
        PSI : petsc Vec, None, optional
            streamfunction, use state.PSI if None
        U : petsc Vec, None, optional
            zonal velocity
        V : petsc Vec, None, optional
            meridional velocity
        RHO : petsc Vec, None, optional
            density

        """

        if PSI is None:
            PSI = state.PSI

        # get u/v
        if U is None or V is None:
            # Initialize u=-dpsidy and v=dpsidx
            self.set_uv_from_psi(da, grid, PSI)
        else:
            self._U = U
            self._V = V
        write_nc([self._U,self._V], ['u','v'], '../output/uv.nc', da, grid)

        # get rho
        if RHO is None:
            # Initialize rho=-f0/g dpsidz
            self.set_rho_from_psi(da, grid, PSI)
        else:
            self._RHO = RHO
        write_nc([self._RHO], ['rho'], '../output/rho.nc', da, grid)

        # Initialize jacobian
        self.set_Q(da, grid)

        # compute Q vector divergence
        self.compute_divQ(da, grid)

        # fix boundaries
        self.set_rhs_bdy(da, grid)

        # mask rhs
        self.set_rhs_mask(da, grid)

    def set_uv_from_psi(self, da, grid, PSI):
        """ Compute U & V from Psi:
            U = -dPSIdy    V =  dPSIdx

        Parameters
        ----------
        da : petsc DMDA
            holds the petsc grid
        grid : qgsolver grid object
            grid data holder
        PSI : petsc Vec, None, optional
            streamfunction, use state.PSI if None

        """

        # create global vectors
        self._U = da.createGlobalVec()
        self._V = da.createGlobalVec()

        # create local vectors
        local_PSI  = da.createLocalVec()
        local_D = da.createLocalVec()

        # load vector PSI used to compute U and V
        da.globalToLocal(PSI, local_PSI)
        da.globalToLocal(qg.grid.D, local_D)

        #
        u = da.getVecArray(self._U)
        v = da.getVecArray(self._V)
        psi = da.getVecArray(local_PSI)
        D = da.getVecArray(local_D)

        mx, my, mz = da.getSizes()
        (xs, xe), (ys, ye), (zs, ze) = da.getRanges()

        kdyu = grid._k_dyu
        kdxv = grid._k_dxv

        # set u=-dpsidy and v=dpsidx
        for k in range(zs,ze):
            for j in range(ys, ye):
                for i in range(xs, xe): 
                    if (i==0    or j==0 or i==mx-1 or j==my-1):
                        # lateral boundaries
                        u[i, j, k] = 0.
                        v[i, j, k] = 0.
                    else:
                        u[i,j,k] = - 1. /D[i,j,kdyu] * \
                             ( 0.25*(psi[i+1,j,k]+psi[i+1,j+1,k]+psi[i,j+1,k]+psi[i,j,k]) -
                               0.25*(psi[i+1,j-1,k]+psi[i+1,j,k]+psi[i,j,k]+psi[i,j-1,k]) )
                        v[i,j,k] =   1. /D[i,j,kdxv] * \
                             ( 0.25*(psi[i+1,j,k]+psi[i+1,j+1,k]+psi[i,j+1,k]+psi[i,j,k]) -
                               0.25*(psi[i,j,k]+psi[i,j+1,k]+psi[i-1,j+1,k]+psi[i-1,j,k]) )

    def set_rho_from_psi(self, da, grid, PSI):
        """ Compute RHO from Psi
            rho=-rho0*f0/g dPSIdz

        Parameters
        ----------
        da : petsc DMDA
            holds the petsc grid
        grid : qgsolver grid object
            grid data holder
        rho0 : float
            1000kg/m^3
        g : float
            9.81 m/s^2
        PSI : petsc Vec, None, optional
            streamfunction, use state.PSI if None

        """
     
        # declare local vectors
        local_PSI  = da.createLocalVec()

        # declare global vectors
        self._RHO = da.createGlobalVec()

        # load vector PSI used to compute RHO
        da.globalToLocal(PSI, local_PSI)

        #
        psi = da.getVecArray(local_PSI)
        rho = da.getVecArray(self._RHO)

        mx, my, mz = da.getSizes()
        (xs, xe), (ys, ye), (zs, ze) = da.getRanges()

        # set rho=-rho0*f0/g dpsidz
        for k in range(zs,ze):
            for j in range(ys, ye):
                for i in range(xs, xe): 
                    if (k==0    or k==mz-1):
                        # top and bottom boundaries
                        rho[i, j, k] = 0.
                    else:
                        rho[i,j,k] = - self.rho0*self.f0/self.g/grid.dzt[k] * \
                             ( 0.5*(psi[i,j,k+1]+psi[i,j,k]) - \
                               0.5*(psi[i,j,k]+psi[i,j,k-1]) )
   
    def set_Q(self, da, grid, U=None, V=None, RHO=None):
        """ Compute Q vector
            qxu = g/f0/rho0 * (dudx*drhodx + dvdx*drhody) at u point
            qyv = g/f0/rho0 * (dudy*drhodx + dvdy*drhody) at v point

        Parameters
        ----------
        da : petsc DMDA
            holds the petsc grid
        grid : qgsolver grid object
            grid data holder
        U : petsc Vec, None, optional
            zonal velocity
        V : petsc Vec, None, optional
            meridional velocity
        RHO : petsc Vec, None, optional
            density

        """

        # declare local vectors
        local_U = da.createLocalVec()
        local_V = da.createLocalVec()
        local_RHO = da.createLocalVec()
        local_D = da.createLocalVec()

        # declare global vectors
        self._QXU = da.createGlobalVec()
        self._QYV = da.createGlobalVec()

        # load vector U,V,RHO used to compute the jacobian
        da.globalToLocal(grid.D, local_D)
        if U is None:
            da.globalToLocal(self._U, local_U)
        else:
            da.globalToLocal(U, local_U)
        if V is None:
            da.globalToLocal(self._V, local_V)
        else:
            da.globalToLocal(V, local_V)
        if RHO is None:
            da.globalToLocal(self._RHO, local_RHO)
        else:
            da.globalToLocal(RHO, local_RHO)
        #
        u = da.getVecArray(local_U)
        v = da.getVecArray(local_V)
        rho = da.getVecArray(local_RHO)
        D = da.getVecArray(local_D)
        qxu = da.getVecArray(self._QXU)
        qyv = da.getVecArray(self._QYV)

        mx, my, mz = da.getSizes()
        (xs, xe), (ys, ye), (zs, ze) = da.getRanges()

        kdxu = grid._k_dxu
        kdyu = grid._k_dyu
        kdxv = grid._k_dxv
        kdyv = grid._k_dyv

        # qxu = g/f0/rho0 * (dudx*drhodx + dvdx*drhody) at u point
        # qyv = g/f0/rho0 * (dudy*drhodx + dvdy*drhody) at v point 
        for k in range(zs,ze):
            for j in range(ys, ye):
                for i in range(xs, xe): 
                    if (i==0    or j==0 or i==mx-1 or j==my-1 ):
                        # lateral boundaries
                        qxu[i, j, k] = 0.
                        qyv[i, j, k] = 0. 
                    else:
                        qxu[i,j,k] = self.g/self.f0/self.rho0 * (
                              (0.5*(u[i+1,j,k]+u[i,j,k])-0.5*(u[i,j,k]+u[i-1,j,k]))/D[i,j,kdxu] *
                              (rho[i+1,j,k]-rho[i,j,k])/D[i,j,kdxu] +
                              (0.5*(v[i+1,j,k]+v[i+1,j-1,k])-0.5*(v[i,j,k]+v[i,j-1,k]))/D[i,j,kdxu] *
                              (0.25*(rho[i+1,j,k]+rho[i+1,j+1,k]+rho[i,j+1,k]+rho[i,j,k]) -
                               0.25*(rho[i+1,j-1,k]+rho[i+1,j,k]+rho[i,j,k]+rho[i,j-1,k]) )/D[i,j,kdyu])
                        qyv[i,j,k] = self.g/self.f0/self.rho0 * (
                              (0.5*(u[i,j+1,k]+u[i-1,j+1,k])-0.5*(u[i,j,k]+u[i-1,j,k]))/D[i,j,kdyv] *
                              (0.25*(rho[i+1,j,k]+rho[i+1,j+1,k]+rho[i,j+1,k]+rho[i,j,k]) -
                               0.25*(rho[i,j,k]+rho[i,j+1,k]+rho[i-1,j+1,k]+rho[i-1,j,k]))/D[i,j,kdxv] +
                              (0.5*(v[i,j,k]+v[i,j+1,k]) - 0.5*(v[i,j,k]+v[i,j-1,k]))/D[i,j,kdyv] *
                              (rho[i,j+1,k]-rho[i,j,k])/D[i,j,kdyv] )

    def compute_divQ(self, da, grid):
        """ Compute Q vector divergence
        
        Parameters
        ----------
        da : petsc DMDA
            holds the petsc grid
        grid : qgsolver grid object
            grid data holder
        
        """
        # declare local vectors
        local_QXU = da.createLocalVec()
        local_QYV = da.createLocalVec()
        local_D = da.createLocalVec()
        #
        da.globalToLocal(self._QXU, local_QXU)
        da.globalToLocal(self._QYV, local_QYV)
        da.globalToLocal(grid.D, local_D)
        #
        qxu = da.getVecArray(local_QXU)
        qyv = da.getVecArray(local_QYV)
        D = da.getVecArray(local_D)

        rhs = da.getVecArray(self._RHS)

        mx, my, mz = da.getSizes()
        (xs, xe), (ys, ye), (zs, ze) = da.getRanges()

        kdyu = grid._k_dyu
        kdxv = grid._k_dxv
        kdxt = grid._k_dxt
        kdyt = grid._k_dyt

        for k in range(zs,ze):
            for j in range(ys, ye):
                for i in range(xs, xe):
                    if (i==0    or j==0 or i==mx-1 or j==my-1 or k==0    or k==mz-1 ):
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
                        rhs[i,j,k] = 2.*self.f0 / D[i,j,kdxt] / D[i,j,kdyt] * (
                                     (D[i,j,kdyu]*qxi - D[i-1,j,kdyu]*qxim) +
                                     (D[i,j,kdxv]*qyj - D[i,j-1,kdxv]*qyjm) )

    def _set_rhs_bdy(self, da, grid):
        """
        Set South/North, East/West, Bottom/Top boundary conditions
        Set RHS along boundaries for inversion, may be an issue
        for time stepping
        Mainly a wrapper
        """
        
        if self._verbose>0:
            print('  set RHS along boudaries for inversion ')

        self._set_rhs_bdy_bottom(da, grid)
        self._set_rhs_bdy_top(da, grid)
        self._set_rhs_bdy_lat(da, grid)

    def _set_rhs_bdy_bottom(self, da, grid, W=None):
        """ Set bottom boundary condition
        """

        rhs = da.getVecArray(self._RHS)
        (xs, xe), (ys, ye), (zs, ze) = da.getRanges()

        kdown = grid.kdown

        # load vector used to compute boundary conditions
        if W is None:
            w = da.getVecArray(qg.W)
        else:
            w = da.getVecArray(W)
           
        # lower ghost area
        if zs < kdown:
            for k in range(zs,kdown):
                for j in range(ys, ye):
                    for i in range(xs, xe):                    
                        # rhs[i,j,k]=sys.float_info.epsilon
                        rhs[i,j,k] = w[i, j, k]
        # bottom bdy
        k = kdown
        if self.bdy_type['bottom']=='N' :
            for j in range(ys, ye):
                for i in range(xs, xe):
                    rhs[i, j, k] = w[i, j, k]
        elif self.bdy_type['bottom']=='NBG' :
            for j in range(ys, ye):
                for i in range(xs, xe):
                    rhs[i, j, k] = w[i, j, k]
        elif self.bdy_type['bottom']=='D':
            for j in range(ys, ye):
                for i in range(xs, xe):
                    rhs[i, j, k] = w[i, j, k]
        else:
            print(self.bdy_type['bottom']+' unknown bottom boundary condition')
            sys.exit()

    def _set_rhs_bdy_top(self, da, grid, W=None):
        """ Set top boundary condition
        """
        
        rhs = da.getVecArray(self._RHS)
        (xs, xe), (ys, ye), (zs, ze) = da.getRanges()

        kup = grid.kup

        # load vector used to compute boundary conditions
        if W is None:
            w = da.getVecArray(qg.W)
        else:
            w = da.getVecArray(W)

        if ze > kup+1:
            for k in range(kup+1,ze):
                for j in range(ys, ye):
                    for i in range(xs, xe):
                        # rhs[i,j,k]=sys.float_info.epsilon   
                        rhs[i,j,k] = w[i, j, k]       
        # upper bdy
        k = kup
        if self.bdy_type['top']=='N' :
            for j in range(ys, ye):
                for i in range(xs, xe):
                    rhs[i, j, k] = w[i, j, k]
        elif self.bdy_type['top']=='NBG' :
            for j in range(ys, ye):
                for i in range(xs, xe):
                    rhs[i, j, k] = w[i, j, k]
        elif self.bdy_type['top']=='D' :
            for j in range(ys, ye):
                for i in range(xs, xe):
                    rhs[i, j, k] = w[i, j, k]

        else:
            print(self.bdy_type['top']+' unknown top boundary condition')
            sys.exit()

    def _set_rhs_bdy_lat(self, da, grid, W=None):
        """ Set lateral boundary condition
        """
        
        rhs = da.getVecArray(self._RHS)
        (xs, xe), (ys, ye), (zs, ze) = qg.da.getRanges()

        istart = qg.grid.istart
        iend = qg.grid.iend
        jstart = qg.grid.jstart
        jend = qg.grid.jend

        # load vector used to compute boundary conditions
        if W is None:
            w = da.getVecArray(qg.W)
        else:
            w = da.getVecArray(W)

        # south bdy
        if ys <= jstart:
            #j = 0
            for k in range(zs, ze):
                for j in range(ys,min(ye,jstart+1)):
                    for i in range(xs, xe):
                        rhs[i, j, k] = w[i, j, k]
        # north bdy
        if ye >= jend:
            #j = my - 1
            for k in range(zs, ze):
                for j in range(max(ys,jend),ye):
                    for i in range(xs, xe):
                        rhs[i, j, k] = w[i, j, k]
        # west bdy
        #if xs <= istart:
        if xs <= istart and self.petscBoundaryType is not 'periodic':
            #i = 0
            for k in range(zs, ze):
                for j in range(ys, ye):
                    for i in range(xs,min(xe,istart+1)):
                        rhs[i, j, k] = w[i, j, k]
        # east bdy
        #if xe >= iend:
        if xe >= iend and self.petscBoundaryType is not 'periodic':
            #i = mx - 1
            for k in range(zs, ze):
                for j in range(ys, ye):
                    for i in range(max(xs,iend),xe):
                        rhs[i, j, k] = w[i, j, k]

    def _set_rhs_mask(self, da, grid):
        """ Set mask on rhs: where mask=0 (land) rhs=psi
        """

        rhs = da.getVecArray(self._RHS)
        mask = da.getVecArray(grid.D)
        (xs, xe), (ys, ye), (zs, ze) = da.getRanges()

        kmask = grid._k_mask

        w = da.getVecArray(qg.W)

        # interior
        for k in range(zs,ze):
            for j in range(ys, ye):
                for i in range(xs, xe):
                    if mask[i,j,kmask]==0.:
                        rhs[i, j, k] = w[i, j, k]

        if self._verbose>0:
            print('  set RHS mask for inversion ')

#
# ==================== Define elliptical operators ===================================
#

    def _set_L(self, L, grid, da):
        """ Builds the laplacian operator along with boundary conditions
            Horizontally uniform grid

        Parameters
        ----------
        L : petsc Mat
            potential vorticity operator
        da : petsc DMDA
            holds the petsc grid
        grid : qgsolver grid object
            grid data holder

        """
        
        if self._verbose>0:
            print('  ... assumes a uniform horizontal and vertical grid')
        
        #
        dx, dy, dz = grid.dx, grid.dy, grid.dz
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

                    # lateral points outside the domain: dirichlet, w=...
                    if ( (i<=istart and self.petscBoundaryType is not 'periodic')
                           or (i>=iend and self.petscBoundaryType is not 'periodic')
                           or j<=jstart or j>=jend):
                        L.setValueStencil(row, row, 1.0)

                    # bottom bdy condition: default Neuman dw/dz=...
                    elif (k==kdown):
                        if self.bdy_type['bottom']=='N' :
                            L.setValueStencil(row, row, 0.0)
                        elif self.bdy_type['bottom']=='D':
                            L.setValueStencil(row, row, 0.0)
                        else:
                            print('unknown bottom boundary condition')
                            sys.exit()
    
                    # top bdy condition: default Neuman dw/dz=...
                    elif (k==kup):
                        if self.bdy_type['top']=='N' :
                            L.setValueStencil(row, row, 0.0)
                        elif self.bdy_type['top']=='D':
                            L.setValueStencil(row, row, 0.0)
                        else:
                            print('unknown top boundary condition')
                            sys.exit()
    
                    # points below and above the domain
                    elif (k<kdown or k>kup):
                        L.setValueStencil(row, row, 0.0)
    
                    # interior points: Q div is prescribed
                    else:
                        for index, value in [
                                ((i,j,k-1), self.f0**2*idz2),
                                ((i,j-1,k), self.N2[k]*idy2),
                                ((i-1,j,k), self.N2[k]*idx2),
                                ((i, j, k), -2.*self.N2[k]*(idx2+idy2)-2.*self.f0**2*idz2),
                                ((i+1,j,k), self.N2[k]*idx2),
                                ((i,j+1,k), self.N2[k]*idy2),
                                ((i,j,k+1), self.f0**2*idz2)]:
                            col.index = index
                            col.field = 0
                            L.setValueStencil(row, col, value)
        L.assemble()

    def _set_L_curv(self,L, da, grid):
        """ Builds the laplacian operator along with boundary conditions
            Curvilinear grid

        Parameters
        ----------
        L : petsc Mat
            potential vorticity operator
        da : petsc DMDA
            holds the petsc grid
        grid : qgsolver grid object
            grid data holder

        """
        
        if self._verbose>0:
            print('  ... assumes a curvilinear and/or vertically stretched grid')
        #
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
        (xs, xe), (ys, ye), (zs, ze) = qg.da.getRanges()
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
                    elif ( (i<=istart and self.petscBoundaryType is not 'periodic')
                           or (i>=iend and self.petscBoundaryType is not 'periodic')
                           or j<=jstart or j>=jend):
                        L.setValueStencil(row, row, 1.0)              
    
                    # bottom bdy condition: default Neuman dpsi/dz=...
                    elif (k==kdown):
                        if self.bdy_type['bottom']=='N' :
                            L.setValueStencil(row, row, 1.0)         
                        elif self.bdy_type['bottom']=='D' :
                            L.setValueStencil(row, row, 1.0)        
                        else:
                            print('unknown bottom boundary condition')
                            sys.exit()
    
                    # top bdy condition: default Neuman dpsi/dz=...
                    elif (k==kup):
                        if self.bdy_type['top']=='N' :
                            L.setValueStencil(row, row, 1.0)        
                        elif self.bdy_type['top']=='D':
                            L.setValueStencil(row, row, 1.0)      
                        else:
                            print('unknown top boundary condition')
                            sys.exit()

                    # points below and above the domain
                    elif (k<kdown) or (k>kup):
                        L.setValueStencil(row, row, 1.0)       

                    # interior points: Q div is prescribed
                    else:
                        
                        for index, value in [
                                ((i,j,k-1), self.f0**2*idzt[k]*idzw[k]),
                                ((i,j-1,k), self.N2[k]/D[i,j,kdxt]/D[i,j,kdyt] * D[i,j-1,kdxv]/D[i,j-1,kdyv]),
                                ((i-1,j,k), self.N2[k]/D[i,j,kdxt]/D[i,j,kdyt] * D[i-1,j,kdyu]/D[i-1,j,kdxu]),
                                ((i, j, k), -self.N2[k]/D[i,j,kdxt]/D[i,j,kdyt]*( \
                                                 D[i,j,kdyu]/D[i,j,kdxu] \
                                                +D[i-1,j,kdyu]/D[i-1,j,kdxu] \
                                                +D[i,j,kdxv]/D[i,j,kdyv] \
                                                +D[i,j-1,kdxv]/D[i,j-1,kdyv])
                                 - (self.f0**2*idzt[k+1]*idzw[k]+self.f0**2*idzt[k]*idzw[k])),
                                ((i+1,j,k), self.N2[k]/D[i,j,kdxt]/D[i,j,kdyt] * D[i,j,kdyu]/D[i,j,kdxu]),
                                ((i,j+1,k), self.N2[k]/D[i,j,kdxt]/D[i,j,kdyt] * D[i,j,kdxv]/D[i,j,kdyv]),
                                ((i,j,k+1), self.f0**2*idzt[k+1]*idzw[k])]:
                            col.index = index
                            col.field = 0
                            L.setValueStencil(row, col, value)
        L.assemble()

