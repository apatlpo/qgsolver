#!/usr/bin/python
# -*- encoding: utf8 -*-

import sys
import petsc4py
#from Cython.Compiler.Main import verbose
petsc4py.init(sys.argv)
#
from .state import *
from .pvinv import *
from .omegainv import *
from .timestepper import *
from .inout import write_nc


class qg_model():
    '''
    QG model
    '''


#
#==================== Object init ============================================
#

    def __init__(self,
                 ncores_x=None, ncores_y=None,
                 hgrid = None, vgrid=None,
                 vdom={}, hdom={}, mask=False,
                 boundary_types={},
                 N2 = 1.e-3, f0 = 7.e-5,
                 f0N2_file = None,
                 dt = None, K = 1.e2,
                 verbose = 1,
                 flag_pvinv=True,
                 flag_omega=False,
                 **kwargs
                 ):
        '''
        QG model initializer

        Parameters
        ----------
        ncores_x : int
            number of MPI tilings in x direction
        ncores_y : int
            number of MPI tilings in y direction
        hgrid : dict or str
            defines horizontal grid choice, see grid doc
        vgrid : dict or str
            defines vertical grid choice, see grid doc
        hdom : dict
            select subdomain, see grid doc
        vdom : dict
            select subdomain, see grid doc
        boundary_types : dict
            may be used to turn on periodic boundary conditions {'periodic'}
        N2 : float, optional
            Brunt Vaisala frequency, default=1.e-3
        f0 : float, optional
            Coriolis frequency, default=7e-5
        f0N2_file : str
            netcdf file containing N2 and f0
        dt : float, optional
            time step
        K : float, optional
            dissipation coefficient, default = 1.e2
        verbose : int, optional
            degree of verbosity, 0 means no outputs
        flag_pvinv: boolean, optional
            turn on setup of PV inversion solver, default is True
        flag_omega: boolean, optional
            turn on setup of omega equation inversion solver, default is False
        '''

        #
        # Build grid object
        #
        self.grid = grid(hgrid, vgrid, hdom, vdom, mask=mask, verbose=verbose)

        # set boundary conditions
        if ('periodic' in boundary_types.keys()) and (boundary_types['periodic']):
            self.petscBoundaryType = 'periodic'
        else:
            self.petscBoundaryType = None
        # default top and bottom boudary condition = 'N_PSI', i.e. Neumann using PSI to compute boundary values.
        # Other possibility 'D' for Dirichlet
        self.bdy_type = {'top':'N_PSI','bottom':'N_PSI'}
        self.bdy_type.update(boundary_types)


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
            self.grid.load_metric_terms(self.da)

        if self.grid.mask:
            # initialize mask
            self.grid.load_mask(self.grid.hgrid_file, self.da)

        #
        if self._verbose>0:
            # print out grid parameters
            print(self.grid)
            # periodicity
            if self._verbose and self.petscBoundaryType is 'periodic':
                print('Petsc boundaries are periodic')

        #
        # create an ocean state
        #
        self.state = state(self.da, self.grid, N2=N2, f0=f0, f0N2_file=f0N2_file, verbose=self._verbose)

        #
        # Init solvers
        #

        # initiate pv inversion solver
        if flag_pvinv:
            self.pvinv = pvinversion(self.da, self.grid, self.bdy_type, sparam=self.state._sparam,
                                     verbose=self._verbose, **kwargs)

        # initiate omega inversion
        if flag_omega:
            self.W = self.da.createGlobalVec()
            self.omegainv = omegainv(self.da, self.grid, self.bdy_type, self.state.f0, self.state.N2,
                                     verbose=self._verbose)

        # initiate time stepper
        if dt is not None:
            self.tstepper = time_stepper(self.da, self.grid, dt, K, self.petscBoundaryType, verbose=self._verbose)



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
                                      stencil_width=2, boundary_type=self.petscBoundaryType)
        # http://lists.mcs.anl.gov/pipermail/petsc-dev/2016-April/018889.html

        self.comm = self.da.getComm()
        self.rank = self.comm.getRank()

#
# ==================== Wrappers to set values of critical variables ===================================
#

    def __getitem__(self,key):
        if hasattr(self.state,key):
            return getattr(self.state,key)
        else:
            print(key+' not in qg.state')

#
#==================== Wrappers to set values of critical variables ===================================
#

    def set_psi(self, **kwargs):
        ''' Set psi, wrapper around state.set_psi
        '''
        self.state.set_psi(self.da, self.grid, **kwargs)


    def set_q(self, **kwargs):
        ''' Set q, wrapper around state.set_q
        '''
        self.state.set_q(self.da, self.grid, **kwargs)

    def set_rho(self, **kwargs):
        ''' Set rho, wrapper around state.set_rho
        '''
        self.state.set_rho(self.da, self.grid, **kwargs)

    def set_bstate(self,**kwargs):
        ''' Set background state
        '''
        if self._verbose:
            print('Set background state:')
        bstate = add(self.state, self.state, da=self.da, a1=0., a2=0.)
        bstate.set_psi(self.da, self.grid, **kwargs)
        bstate.set_rho(self.da, self.grid, **kwargs)
        bstate.set_q(self.da, self.grid, **kwargs) # read from file
        #self.pvinv.q_from_psi(bstate.Q,bstate.PSI)
        return bstate

    def set_w(self, **kwargs):
        ''' Set w, wrapper around state.set_w
        '''
        self.state.set_w(self.da, self.grid, **kwargs)

#
#==================== useful wrappers for solvers ============================================
#

    def invert_pv(self, **kwargs):
        ''' wrapper around pv inversion solver pvinv.solve
        '''
        if hasattr(self,'state'):
            if hasattr(self.state,'RHO'):
                RHO = self.state.RHO
            else:
                RHO = None
            dkwargs = {'Q': self.state.Q, 'PSI': self.state.PSI, 'RHO': RHO,
                       'bstate': None, 'addback_bstate': True}
            dkwargs.update(kwargs)
            self.pvinv.solve(self.da, self.grid, self.state, **dkwargs)
        else:
            print('!Error qg.inver_pv requires qg.state (with Q/PSI and RHO depending on bdy conditions)')

    def invert_omega(self):
        ''' wrapper around solver solve method omegainv.solve
        '''
        self.omegainv.solve(self)

    def tstep(self, nt=1, rho_sb=True, bstate=None):
        ''' Time step wrapper tstepper.go
        '''
        self.tstepper.go(nt, self.da, self.grid, self.state, self.pvinv, rho_sb, bstate=bstate)


#
# ==================== IO ============================================
#

    def write_state(self,v=['PSI','Q'], vname=['psi','q'], filename='output.nc', append=False):
        ''' Outputs state to a netcdf file

        Parameters
        ----------
        v : list of str
            List of variables to output (must be contained in state object)
        vname : list of str
            list of the names used in netcdf files
        filename : str
            netcdf output filename
        create : boolean, optional
            if true creates a new file, append otherwise (default is True)
        '''
        V=[]
        for vv in v:
            if hasattr(self.state,vv):
                V.append(getattr(self.state,vv))
            else:
                print('Warning: variable '+vv+' not present in state vector and thus not outputted')
        write_nc(V, vname, filename, self.da, self.grid, append=append)

#
#==================== utils ============================================
#
    def compute_CFL(self, PSI=None):
        ''' Compute CFL = max (u*dt/dx)

        Parameters
        ----------
        PSI: petsc Vec, optional
            PSI vector used for velocity computation

        Returns
        -------
        CFL: float
            CFL number
        '''

        # compute U from psi
        self.state.get_uv(self.da, self.grid, PSI=PSI)

        # compute abs(u*dt/dx)
        self._compute_dudx(PSI=PSI)

        CFL=self.state._U.max()[1]*self.tstepper.dt
        self.state._U.destroy()
        return CFL

    def _compute_dudx(self, PSI=None):
        ''' Compute abs(u*dt/dx)
        '''
        # get u
        u = self.da.getVecArray(self.state._U)
        # get dx
        D = self.da.getVecArray(self.grid.D)

        kdxu = self.grid._k_dxu
        (xs, xe), (ys, ye), (zs, ze) = self.da.getRanges()

        for k in range(zs,ze):
            u[:,:,k] = u[:,:,k]/D[:,:,kdxu]

    def compute_KE(self, PSI=None):
        ''' Compute the domain averaged kinetic energy, wrapper around state.compute_KE

        Parameters
        ----------
        PSI: petsc Vec, optional
            PSI vector used for velocity computation

        Returns
        -------
        KE: float
            Kinetic energy in m/s
        '''
        return self.state.compute_KE(self.da, self.grid, PSI=PSI)
