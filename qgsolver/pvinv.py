#!/usr/bin/python
# -*- encoding: utf8 -*-


import sys
from petsc4py import PETSc
from .utils import g, rho0

#
#==================== PV inversion solver object ============================================
#


class pvinversion():
    ''' PV inversion solver
    '''
    
    def __init__(self, da, grid, bdy_type, sparam, verbose=0, solver='gmres', pc=None):
        ''' Setup the PV inversion solver

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
        sparam : ndarray
            numpy array containing f^2/N^2
        verbose : int, optional
            degree of verbosity, 0 means no outputs
        solver : str, optional
            petsc solver: 'gmres' (default), 'bicg', 'cg'
        pc : str, optional
            what is default?
            preconditionner: 'icc', 'bjacobi', 'asm', 'mg', 'none'

        '''

        self._verbose = verbose
        #
        self.bdy_type = bdy_type
        if ('periodic' in self.bdy_type) and (self.bdy_type['periodic']):
            self.petscBoundaryType = 'periodic'
        else:
            self.petscBoundaryType = None

        # create the operator
        self.L = da.createMat()
        #
        if self._verbose>0:
            print('A PV inversion object is being created')
            print('  Operator L declared')

        # Fill in operator values
        if grid._flag_hgrid_uniform and grid._flag_vgrid_uniform:
            self._set_L(self.L, da, grid, sparam)
        else:
            self._set_L_curv(self.L, da, grid, sparam)

        #
        if self._verbose>0:
            print('  Operator L filled')

        # global vector for PV inversion
        self._RHS = da.createGlobalVec()

        # create solver
        self.ksp = PETSc.KSP()
        self.ksp.create(PETSc.COMM_WORLD)
        self.ksp.setOperators(self.L)
        self.ksp.setType(solver)
        self.ksp.setInitialGuessNonzero(True)
        if pc is not None:
            self.ksp.getPC().setType(pc)
        # set tolerances
        self.ksp.setTolerances(rtol=1e-4)
        self.ksp.setTolerances(max_it=100)
        # tests:
        #self.ksp.setPCSide(PETSc.PC.Side.R)
        #self.ksp.setInitialGuessNonzero(False)
        #self.ksp.setTolerances(atol=1e-1)
        #self.ksp.setTolerances(max_it=100)
        #
        for opt in sys.argv[1:]:
            PETSc.Options().setValue(opt, None)
        self.ksp.setFromOptions()
        
        if self._verbose>0:
            print('  PV inversion is set up')

#
# ==================== perform inversion ===================================
#

    def solve(self, da, grid, state, Q=None, PSI=None, RHO=None, \
              bstate=None, addback_bstate=True, topdown_rho=False, numit=False):
        ''' Compute the PV inversion
        Uses prioritarily optional Q, PSI, RHO for RHS and bdy conditions

        Parameters
        ----------
        da : petsc DMDA
            holds the petsc grid
        grid : qgsolver grid object
            grid data holder
        state : state object
            ocean state
        Q : petsc Vec, None, optional
            potential vorticity, use state.Q if None
        PSI : petsc Vec, None, optional
            streamfunction, use state.PSI if None
        RHO : petsc Vec, None, optional
            density, use state.RHO if None
        bstate : state object, None, optional
            background state that will be added in advective terms
        addback_bstate : boolean
            if True, add background state back to output variables ()
        topdown_rho : boolean
            if True, indicates that RHO used for top down boundary conditions
            is contained in state.Q at indices kdown and kup
        numit : boolean
            if True, returns the number of iterations

        Returns
        -------
        state.PSI:
            Put PV inversion result in state.PSI
        '''
        if Q is None and not hasattr(state,'Q'):
            print('!Error: pvinv.solve requires state.Q or Q')
            sys.exit()
        elif Q is None:
            Q = state.Q
        if PSI is None and not hasattr(state, 'PSI'):
            print('!Error: pvinv.solve requires state.PSI or PSI')
            sys.exit()
        elif PSI is None:
            PSI = state.PSI
        if RHO is None and not hasattr(state, 'RHO'):
            print('!Error: pvinv.solve requires state.RHO or RHO')
            sys.exit()
        elif RHO is None:
            if topdown_rho:
                RHO = Q
            else:
                RHO = state.RHO
        #
        # substract background fields
        if bstate is not None:
            Q += - bstate.Q
            PSI += - bstate.PSI
            RHO += - bstate.RHO
        #
        # copy Q into RHS
        Q.copy(self._RHS)
        # fix boundaries
        self.set_rhs_bdy(da, grid, state, PSI, RHO, topdown_rho)
        # apply mask
        if grid.mask:
            # mask rhs
            self.set_rhs_mask(da, grid, PSI)
        # actually solves the pb
        self.ksp.solve(self._RHS, state.PSI)
        # add back background state
        if bstate is not None and addback_bstate:
            if self._verbose>1:
                print('add back bstate after pv inversion')
            Q += bstate.Q
            state.PSI += bstate.PSI
            RHO += bstate.RHO

        if self._verbose>1:
            print('Inversion done (%i iterations)' %(self.ksp.getIterationNumber()), flush=True)
        if numit:
            return self.ksp.getIterationNumber()


#
# ==================== utils methods for inversions ===================================
#

    def q_from_psi(self, Q, PSI):
        ''' Compute PV from a streamfunction
        
        Parameters
        ----------
        Q : petsc Vec
            output vector where data is stored
        PSI : petsc Vec
            input streamfunction used for the computation of PV
        '''
        self.L.mult(PSI, Q)
        # should fix boundary conditions
        

    def set_rhs_bdy(self, da, grid, state, PSI, RHO, topdown_rho):
        ''' Set South/North, East/West, Bottom/Top boundary conditions
        Set RHS along boundaries for inversion, may be an issue
        for time stepping

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
        RHO : petsc Vec, None, optional
            density, use state.RHO if None
        '''
        
        if self._verbose>1:
            print('  Set RHS along boudaries for inversion ')

        self.set_rhs_bdy_bottom(da, grid, state, PSI, RHO)
        self.set_rhs_bdy_top(da, grid, state, PSI, RHO, topdown_rho)
        self.set_rhs_bdy_lat(da, grid, PSI)

    def set_rhs_bdy_bottom(self, da, grid, state, PSI, RHO):
        ''' Set bottom boundary condition

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
        RHO : petsc Vec, None, optional
            density, use state.RHO if None
        '''

        rhs = da.getVecArray(self._RHS)
        (xs, xe), (ys, ye), (zs, ze) = da.getRanges()

        kdown = grid.kdown

        # load vector used to compute boundary conditions
        psi = da.getVecArray(PSI)
        rho = da.getVecArray(RHO)
           
        # lower ghost area
        if zs < kdown:
            for k in range(zs,kdown):
                for j in range(ys, ye):
                    for i in range(xs, xe):                    
                        # rhs[i,j,k]=sys.float_info.epsilon
                        rhs[i,j,k]=psi[i, j, k]
        # bottom bdy
        k = kdown
        if self.bdy_type['bottom'] == 'N':
            for j in range(ys, ye):
                for i in range(xs, xe):
                    rhs[i, j, k] = - g*rho[i, j, k]/(rho0*state.f0)
        elif self.bdy_type['bottom'] == 'NBG':
            for j in range(ys, ye):
                for i in range(xs, xe):
                    rhs[i, j, k] = (psi[i,j,k+1]-psi[i,j,k])/grid.dzw[k]
        elif self.bdy_type['bottom'] == 'D':
            for j in range(ys, ye):
                for i in range(xs, xe):
                    rhs[i, j, k] = psi[i,j,k]
        else:
            print(self.bdy_type['bottom']+" unknown bottom boundary condition")
            sys.exit()
                
        return

    def set_rhs_bdy_top(self, da, grid, state, PSI, RHO, topdown_rho):
        '''Set top boundary condition

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
        RHO : petsc Vec, None, optional
            density, use state.RHO if None
        '''
        
        rhs = da.getVecArray(self._RHS)
        (xs, xe), (ys, ye), (zs, ze) = da.getRanges()

        kup = grid.kup

        # load vector used to compute boundary conditions
        psi = da.getVecArray(PSI)
        rho = da.getVecArray(RHO)

        if ze > kup+1:
            for k in range(kup+1,ze):
                for j in range(ys, ye):
                    for i in range(xs, xe):
                        # rhs[i,j,k]=sys.float_info.epsilon   
                        rhs[i,j,k]= psi[i,j,k]            
        # upper bdy
        k = kup
        if self.bdy_type['top'] == 'N':
            if topdown_rho:
                krho = k
            else:
                krho = k-1
            for j in range(ys, ye):
                for i in range(xs, xe):
                    rhs[i, j, k] = - g*rho[i, j, krho]/(rho0*state.f0)
        elif self.bdy_type['top'] == 'NBG':
            for j in range(ys, ye):
                for i in range(xs, xe):
                    rhs[i, j, k] = (psi[i,j,k]-psi[i,j,k-1])/grid.dzw[k-1]
        elif self.bdy_type['top'] == 'D':
            for j in range(ys, ye):
                for i in range(xs, xe):
                    rhs[i, j, k] = psi[i,j,k]
        else:
            print(self.bdy_type['top']+" unknown top boundary condition")
            sys.exit()

    def set_rhs_bdy_lat(self, da, grid, PSI):
        '''Set lateral boundary condition

        Parameters
        ----------
        da : petsc DMDA
            holds the petsc grid
        grid : qgsolver grid object
            grid data holder
        PSI : petsc Vec, None, optional
            streamfunction, use state.PSI if None
        '''
        
        rhs = da.getVecArray(self._RHS)
        (xs, xe), (ys, ye), (zs, ze) = da.getRanges()

        istart = grid.istart
        iend = grid.iend
        jstart = grid.jstart
        jend = grid.jend

        # load vector used to compute boundary conditions
        psi = da.getVecArray(PSI)

        # south bdy
        if ys <= jstart and self.petscBoundaryType is not 'periodic':
            #j = 0
            for k in range(zs, ze):
                for j in range(ys,min(ye,jstart+1)):
                    for i in range(xs, xe):
                        rhs[i, j, k] = psi[i, j, k]
        # north bdy
        if ye >= jend and self.petscBoundaryType is not 'periodic':
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

    def set_rhs_mask(self, da, grid, PSI):
        '''Set mask on rhs: where mask=0 (land) rhs=psi

        Parameters
        ----------
        da : petsc DMDA
            holds the petsc grid
        grid : qgsolver grid object
            grid data holder
        PSI : petsc Vec
            streamfunction used over masked areas
        '''

        rhs = da.getVecArray(self._RHS)
        mask = da.getVecArray(grid.D)
        (xs, xe), (ys, ye), (zs, ze) = da.getRanges()

        kmask = grid._k_mask

        psi = da.getVecArray(PSI)
        
        # interior
        for k in range(zs,ze):
            for j in range(ys, ye):
                for i in range(xs, xe):
                    if mask[i,j,kmask] == 0.:
                        rhs[i, j, k] = psi[i,j,k]

        if self._verbose>1:
            print('  Set RHS mask for inversion ')

#
# ==================== Define elliptical operators ===================================
#
    
    def _set_L(self,L, da, grid, sparam):
        ''' Builds the laplacian operator along with boundary conditions
            Horizontally uniform grid

        Parameters
        ----------
        L : petsc Mat
            potential vorticity operator
        da : petsc DMDA
            holds the petsc grid
        grid : qgsolver grid object
            grid data holder
        sparam : ndarray
            f0^2/N^2/ array
        '''
        
        if self._verbose>0:
            print('  ... assumes a uniform horizontal and vertical grid')
        
        #
        dx, dy, dz = grid.dx, grid.dy, grid.dz
        idx, idy, idz = [1.0/dl for dl in [dx, dy, dz]]
        idx2, idy2, idz2 = [1.0/dl**2 for dl in [dx, dy, dz]]
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

                    # lateral points outside the domain: dirichlet, psi=...
                    if (i <= istart or i >= iend or j <= jstart or j >= jend) \
                            and self.petscBoundaryType is not 'periodic':
                        L.setValueStencil(row, row, 1.0)
    
                    # bottom bdy condition: default Neuman dpsi/dz=...
                    elif k == kdown:
                        if self.bdy_type['bottom'] == 'N' or self.bdy_type['bottom'] == 'NBG':
                            for index, value in [
                                    ((i,j,k), -idz),
                                    ((i,j,k+1),  idz)]:
                                col.index = index
                                col.field = 0
                                L.setValueStencil(row, col, value)
                        elif self.bdy_type['bottom'] == 'D':
                            L.setValueStencil(row, row, 1.0)
                        else:
                            print('unknown bottom boundary condition')
                            sys.exit()
    
                    # top bdy condition: default Neuman dpsi/dz=...
                    elif k == kup:
                        if self.bdy_type['top'] == 'N' or self.bdy_type['top'] == 'NBG':
                            for index, value in [
                                    ((i,j,k-1), -idz),
                                    ((i,j,k),  idz)]:
                                col.index = index
                                col.field = 0
                                L.setValueStencil(row, col, value)
                        elif self.bdy_type['top'] == 'D':
                            L.setValueStencil(row, row, 1.0)
                        else:
                            print('unknown top boundary condition')
                            sys.exit()
    
                    # points below and above the domain
                    elif k < kdown or k > kup:
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
                                ((i,j,k+1), sparam[k]*idz2)]:
                            col.index = index
                            col.field = 0
                            L.setValueStencil(row, col, value)
        L.assemble()
        return

    def _set_L_curv(self,L, da, grid, sparam):
        ''' Builds the laplacian operator along with boundary conditions
            Horizontally uniform grid

        Parameters
        ----------
        L : petsc Mat
            potential vorticity operator
        da : petsc DMDA
            holds the petsc grid
        grid : qgsolver grid object
            grid data holder
        sparam : ndarray
            f0^2/N^2/ array
        '''
        
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

                    # row index
                    row.index = (i, j, k)
                    row.field = 0
    
                    # masked points (land=0), L=1
                    if D[i,j,kmask] == 0.:
                        L.setValueStencil(row, row, 1.)
   
                    # domain edges 
                    elif (i <= istart or i >= iend or j <= jstart or j >= jend) \
                            and self.petscBoundaryType is not 'periodic':
                        L.setValueStencil(row, row, 1.0)
    
                    # bottom bdy condition: default Neuman dpsi/dz=...
                    elif k == kdown:
                        if self.bdy_type['bottom'] == 'N' or self.bdy_type['bottom'] == 'NBG':
                            for index, value in [
                                    ((i,j,k), -idzw[k]),
                                    ((i,j,k+1),  idzw[k])]:
                                col.index = index
                                col.field = 0
                                L.setValueStencil(row, col, value)
                        elif self.bdy_type['bottom']=='D' :
                            L.setValueStencil(row, row, 1.0)
                        else:
                            print('unknown bottom boundary condition')
                            sys.exit()
    
                    # top bdy condition: default Neuman dpsi/dz=...
                    elif k == kup:
                        if self.bdy_type['top'] == 'N' or self.bdy_type['top'] == 'NBG':
                            for index, value in [
                                    ((i,j,k-1), -idzw[k-1]),
                                    ((i,j,k),  idzw[k-1])]:
                                col.index = index
                                col.field = 0
                                L.setValueStencil(row, col, value)
                        elif self.bdy_type['top']=='D':
                            L.setValueStencil(row, row, 1.0)
                        else:
                            print('unknown top boundary condition')
                            sys.exit()
    
                    # points below and above the domain
                    elif k < kdown or k > kup:
                        L.setValueStencil(row, row, 1.0)

                    # interior points: pv is prescribed
                    else:
                        
                        for index, value in [
                                ((i,j,k-1), sparam[k-1]*idzt[k]*idzw[k-1]),
                                ((i,j-1,k), 1./D[i,j,kdxt]/D[i,j,kdyt] * D[i,j-1,kdxv]/D[i,j-1,kdyv]),
                                ((i-1,j,k), 1./D[i,j,kdxt]/D[i,j,kdyt] * D[i-1,j,kdyu]/D[i-1,j,kdxu]),
                                ((i, j, k), -1./D[i,j,kdxt]/D[i,j,kdyt]*(
                                                 D[i,j,kdyu]/D[i,j,kdxu]
                                                +D[i-1,j,kdyu]/D[i-1,j,kdxu]
                                                +D[i,j,kdxv]/D[i,j,kdyv]
                                                +D[i,j-1,kdxv]/D[i,j-1,kdyv])
                                            - (sparam[k]*idzt[k]*idzw[k]+sparam[k-1]*idzt[k]*idzw[k-1])),
                                ((i+1,j,k), 1./D[i,j,kdxt]/D[i,j,kdyt] * D[i,j,kdyu]/D[i,j,kdxu]),
                                ((i,j+1,k), 1./D[i,j,kdxt]/D[i,j,kdyt] * D[i,j,kdxv]/D[i,j,kdyv]),
                                ((i,j,k+1), sparam[k]*idzt[k]*idzw[k])]:
                            col.index = index
                            col.field = 0
                            L.setValueStencil(row, col, value)
        #
        L.assemble()
