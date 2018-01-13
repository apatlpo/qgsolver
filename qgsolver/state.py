#!/usr/bin/python
# -*- encoding: utf8 -*-

from .grid import *
from .inout import read_nc, read_nc_petsc
from .utils import g, rho0

class state():
    '''
    Ocean state variable holder
    '''

    def __init__(self, da, grid, N2=1.e-3, f0=7e-5, f0N2_file=None, verbose=0):
        """Declare Petsc vectors

        Parameters
        ----------
        da : petsc DMDA
            holds the petsc grid
        grid : qgsolver grid object
            grid data holder
        N2 : float, optional
            Brunt Vaisala frequency, default=1.e-3
        f0 : float, optional,
            Coriolis frequency, default=7.e-5
        f0N2_file : str, optional
            netcdf file containing N2 and f0, default is None
        verbose : int, optional
            degree of verbosity, 0 means no outputs

        """

        self._verbose = verbose

        # PV
        self.Q = da.createGlobalVec()
        # streamfunction
        self.PSI = da.createGlobalVec()
        # density
        self.RHO = da.createGlobalVec()

        #
        # vertical stratification and Coriolis
        #
        # N2 is at w points (cell faces), N2[k] is between q[k] and q[k+1]
        if f0N2_file is not None:
            if self._verbose>0:
                print('  Reads N2, f0 and f from '+f0N2_file)
            #
            self.N2 = read_nc('N2', f0N2_file, grid)
            self.f0 = read_nc('f0', f0N2_file, grid)
            grid.load_coriolis_parameter(f0N2_file, da)
            #
            self._compute_sparam()
        elif (N2 is not None) :
            if self._verbose>0:
                print('  Set N2 from user prescribed value = '+str(N2)+' 1/s^2')
                print('  Sets f0 to %.3e' % f0)
            #
            self.N2 = N2*np.ones(grid.Nz)
            self.f0 = f0
            #
            self._compute_sparam()

    def __str__(self):
        """ Plot useful statistics about the state
        """
        out = '  State vector desribed by: \n' \
              + '    N max = %.2e 1/s  \n'%(np.sqrt(self.N2)) \
              + '    f0 = %.2e \n' %(self.f0)
        return out

    def __getitem__(self,key):
        if hasattr(self,key):
            return getattr(self,key)
        else:
            print(key+' not in state')

    def __setitem__(self, key, item):
        if hasattr(self, key):
            return setattr(self, key, item)
        else:
            print(key + ' not in state')

    def __mul__(self,other):
        return None 

    def __add__(self,other):
        return None

    def __sub__(self,other):
        return None

#
# ==================== fill state variables ============================================
#

    def _compute_sparam(self):
        """ Compute f^2/N^2
        """
        self._sparam = self.f0**2 /self.N2

    def set_psi(self, da, grid, analytical_psi=True, psi0=0., file=None, **kwargs):
        """ Set psi (streamfunction)

        Parameters
        ----------
        da : petsc DMDA
            holds the petsc grid
        grid : qgsolver grid object
            grid data holder
        analytical_psi : boolean, optional
            True set psi analytically, default is True
        file : str, optional
            filename where psi can be found

        """
        if file is not None:
            if self._verbose:
                print('  Set psi from file ' + file + ' ...')
            read_nc_petsc(self.PSI, 'psi', file, da, grid, fillmask=0.)
        elif analytical_psi:
            self.set_psi_analytically(da, psi0)
        else:
            print('Error: you need to provide arguments (filename or turn on analytical definition) ' \
                  +'in order to set psi')
            sys.exit()

    def set_psi_analytically(self, da, psi0):
        """ Set psi analytically

        """
        psi = da.getVecArray(self.PSI)
        (xs, xe), (ys, ye), (zs, ze) = da.getRanges()
        #
        if self._verbose:
            print('  Set psi analytically')
        for k in range(zs, ze):
            for j in range(ys, ye):
                for i in range(xs, xe):
                    psi[i, j, k] = psi0

    def set_q(self, da, grid, analytical_q=True, q0=1.e-5, beta=0., file=None, **kwargs):
        """ Set q (PV)

        Parameters
        ----------
        da : petsc DMDA
            holds the petsc grid
        grid : qgsolver grid object
            grid data holder
        analytical_psi : boolean, optional
            True set psi analytically, default is True
        file : str, optional
            filename where q can be found

        """
        #
        if file is not None:
            if self._verbose:
                print('  Set q from file ' + file + ' ...')
            read_nc_petsc(self.Q, 'q', file, da, grid, fillmask=0.)
        elif analytical_q:
            self.set_q_analytically(da, grid, q0, beta)

    def set_q_analytically(self,da, grid, q0, beta):
        """ Set q analytically
        """
        q = da.getVecArray(self.Q)
        mx, my, mz = da.getSizes()
        (xs, xe), (ys, ye), (zs, ze) = da.getRanges()
        #
        if self._verbose:
            print('  Set q analytically')
        for k in range(zs, ze):
            for j in range(ys, ye):
                for i in range(xs, xe):
                    q[i, j, k] = q0 * np.exp(-((i / float(mx - 1) - 0.5) ** 2
                                                  + (j / float(my - 1) - 0.5) ** 2) / 0.1 ** 2)
                    q[i, j, k] *= np.sin(i / float(mx - 1) * np.pi)
                    q[i, j, k] *= np.sin(2 * j / float(my - 1) * np.pi)
                    q[i, j, k] += beta*grid.dy*(j-my/2.)

    def set_rho(self, da, grid, analytical_rho=True, rho0=0., file=None, **kwargs):
        """ Set rho (density)

        Parameters
        ----------
        da : petsc DMDA
            holds the petsc grid
        grid : qgsolver grid object
            grid data holder
        analytical_psi : boolean, optional
            True set psi analytically, default is True
        file : str, optional
            filename where rho can be found

        """
        #
        if file is not None:
            if self._verbose:
                print('  Set rho from file ' + file + ' ...')
            read_nc_petsc(self.RHO, 'rho', file, da, grid, fillmask=0.)
        elif analytical_rho:
            self.set_rho_analytically(da, rho0)

    def set_rho_analytically(self, da, rho0):
        """ Set rho analytically
        """
        rho = da.getVecArray(self.RHO)
        (xs, xe), (ys, ye), (zs, ze) = da.getRanges()
        #
        if self._verbose:
            print('  Set rho analytically')
        for k in range(zs, ze):
            for j in range(ys, ye):
                for i in range(xs, xe):
                    rho[i, j, k] = rho0

    def set_w(self, da, grid, analytical_w=True, file=None, **kwargs):
        """ Set w

        Parameters
        ----------
        da : petsc DMDA
            holds the petsc grid
        grid : qgsolver grid object
            grid data holder
        analytical_psi : boolean, optional
            True set psi analytically, default is True
        file : str, optional
            filename where w can be found

        """
        #
        if file is not None:
            if self._verbose:
                print('Set w from file ' + file + ' ...')
            read_nc_petsc(self.W, 'w', file, da, grid, fillmask=0.)
        elif analytical_w:
            self.set_w_analytically(da)

    def set_w_analytically(self,da):
        """ Set w analytically
        """
        w = da.getVecArray(self.W)
        mx, my, mz = da.getSizes()
        (xs, xe), (ys, ye), (zs, ze) = da.getRanges()
        #
        if self._verbose:
            print('Set w analytically to zero')
        for k in range(zs, ze):
            for j in range(ys, ye):
                for i in range(xs, xe):
                    w[i, j, k] = 0.


#
# ==================== utils ============================================
#

    def update_rho(self, da, grid, PSI=None, RHO=None):
        """ update rho from psi
        """

        if PSI is None:
            PSI = self.PSI
        if RHO is None:
            RHO = self.RHO
        psi = da.getVecArray(PSI)
        rho = da.getVecArray(RHO)

        #
        idzw = 1. / grid.dzw
        #
        (xs, xe), (ys, ye), (zs, ze) = da.getRanges()
        kdown = grid.kdown
        kup = grid.kup

        for k in range(kdown + 1, kup):
            for j in range(ys, ye):
                for i in range(xs, xe):
                    rho[i, j, k] = -rho0 * self.f0 /g * \
                                   0.5 * ((psi[i, j, k + 1] - psi[i, j, k]) * idzw[k] \
                                          + (psi[i, j, k] - psi[i, j, k - 1]) * idzw[k - 1])
        # extrapolate top and bottom
        k = kdown
        for j in range(ys, ye):
            for i in range(xs, xe):
                rho[i, j, k] = -rho0 * self.f0 /g * (psi[i, j, k + 1] - psi[i, j, k]) * idzw[k]
        k = kup
        for j in range(ys, ye):
            for i in range(xs, xe):
                rho[i, j, k] = -rho0 * self.f0 / g * (psi[i, j, k] - psi[i, j, k - 1]) * idzw[k - 1]
        return

    def get_uv(self, da, grid, PSI=None):
        """ Compute horizontal velocities
        Compute U & V from Psi
        U = -dPSIdy
        V =  dPSIdx
        """

        ### create global vectors
        self._U = da.createGlobalVec()
        self._V = da.createGlobalVec()

        ### create local vectors
        local_PSI = da.createLocalVec()
        local_D = da.createLocalVec()

        #### load vector PSI used to compute U and V
        if PSI is None:
            da.globalToLocal(self.PSI, local_PSI)
        else:
            da.globalToLocal(PSI, local_PSI)

        da.globalToLocal(grid.D, local_D)

        #
        u = da.getVecArray(self._U)
        v = da.getVecArray(self._V)
        psi = da.getVecArray(local_PSI)
        D = da.getVecArray(local_D)

        mx, my, mz = da.getSizes()
        (xs, xe), (ys, ye), (zs, ze) = da.getRanges()

        kdyu = grid._k_dyu
        kdxv = grid._k_dxv

        # Initialize u=-dpsidy and v=dpsidx

        for k in range(zs, ze):
            for j in range(ys, ye):
                for i in range(xs, xe):
                    if (i == 0 or j == 0 or
                                i == mx - 1 or j == my - 1):
                        # lateral boundaries
                        u[i, j, k] = 0.
                        v[i, j, k] = 0.
                    else:
                        u[i, j, k] = - 1. / D[i, j, kdyu] * \
                                     (0.25 * (
                                     psi[i + 1, j, k] + psi[i + 1, j + 1, k] + psi[i, j + 1, k] + psi[i, j, k]) - \
                                      0.25 * (
                                      psi[i + 1, j - 1, k] + psi[i + 1, j, k] + psi[i, j, k] + psi[i, j - 1, k]))
                        v[i, j, k] = 1. / D[i, j, kdxv] * \
                                     (0.25 * (
                                     psi[i + 1, j, k] + psi[i + 1, j + 1, k] + psi[i, j + 1, k] + psi[i, j, k]) - \
                                      0.25 * (
                                      psi[i, j, k] + psi[i, j + 1, k] + psi[i - 1, j + 1, k] + psi[i - 1, j, k]))

    def compute_KE(self, da, grid, PSI=None):
        """ Compute kinetic energy = 0.5 * sum (u**2+v**2)
        """

        # compute local kinetic energy
        self._compute_local_KE(da, grid, PSI=PSI)

        # average spatially
        KE = self._lKE.sum()
        Vol = self._Vol.sum()
        self._lKE.destroy()
        self._Vol.destroy()

        return KE / Vol

    def _compute_local_KE(self, da, grid, PSI=None):

        # create global vectors
        self._lKE = da.createGlobalVec()
        self._Vol = da.createGlobalVec()

        # create local vectors
        local_PSI = da.createLocalVec()
        local_D = da.createLocalVec()

        # load vector PSI used to compute U and V
        if PSI is None:
            da.globalToLocal(self.PSI, local_PSI)
        else:
            da.globalToLocal(PSI, local_PSI)

        da.globalToLocal(grid.D, local_D)

        #
        lKE = da.getVecArray(self._lKE)
        Vol = da.getVecArray(self._Vol)
        psi = da.getVecArray(local_PSI)
        D = da.getVecArray(local_D)

        mx, my, mz = da.getSizes()
        (xs, xe), (ys, ye), (zs, ze) = da.getRanges()

        kdyu = grid._k_dyu
        kdxv = grid._k_dxv
        kdxt = grid._k_dxt
        kdyt = grid._k_dyt

        # Loop around volume
        for k in range(zs, ze):
            for j in range(ys, ye):
                for i in range(xs, xe):
                    if (i == 0 or j == 0 or
                                i == mx - 1 or j == my - 1):
                        # lateral boundaries
                        lKE[i, j, k] = 0.
                        Vol[i, j, k] = 0.
                    else:
                        u = - 1. / D[i, j, kdyu] * \
                            (0.25 * (psi[i + 1, j, k] + psi[i + 1, j + 1, k] + psi[i, j + 1, k] + psi[i, j, k]) - \
                             0.25 * (psi[i + 1, j - 1, k] + psi[i + 1, j, k] + psi[i, j, k] + psi[i, j - 1, k]))
                        v = 1. / D[i, j, kdxv] * \
                            (0.25 * (psi[i + 1, j, k] + psi[i + 1, j + 1, k] + psi[i, j + 1, k] + psi[i, j, k]) - \
                             0.25 * (psi[i, j, k] + psi[i, j + 1, k] + psi[i - 1, j + 1, k] + psi[i - 1, j, k]))
                        Vol[i, j, k] = self.grid.dzt[k] * D[i, j, kdxt] * D[i, j, kdyt]
                        lKE[i, j, k] = 0.5 * (u ** 2 + v ** 2) * Vol[i, j, k]

def add(state1, state2, da=None, a1=1., a2=1., a3=0.):
    """ add fields of two states: a1*state1 + a2*state2 + a3

    Parameters
    ----------
    state1 : qgsolver state
    state2 : qgsolver state
    da : None or petsc DMDA
        if None state1 is updated; otherwise a new state is created
    a1 : float, optional
        default value = 1.
    a2 : float, optional
        default value = 1.
    a2 : float, optional
        default value = 0.

    """
    if state1 is not None and state2 is not None:
        if da is not None:
            new_state = state(da,None,N2=None,verbose=state1._verbose)
        else:
            new_state = state1
        new_state['Q'] = a1 * state1['Q'] + a2 * state2['Q'] + a3
        new_state['PSI'] = a1 * state1['Q'] + a2 * state2['Q'] + a3
        new_state['RHO'] = a1 * state1['Q'] + a2 * state2['Q'] + a3

        # we use state1 parameters
        new_state.N2 = state1.N2
        new_state.f0 = state1.f0
        new_state._compute_sparam()
    else:
        new_state = None

    if da is not None:
        return new_state

