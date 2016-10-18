#!/usr/bin/python
# -*- encoding: utf8 -*-

""" Test the basic features of the library:
Setup of uniform grid
PV inversion of an analytical PV distribution
"""

#import qgsolver.qg as qg
from qgsolver.qg import qg_model
from qgsolver.io import write_nc


def uniform_grid_runs():
    ''' Tests with uniform grid, closed domains
    '''
    qg = qg_model(hgrid = {'Nx':150, 'Ny':100}, vgrid = {'Nz':3 },
            K = 0.e0, dt = 0.5*86400.e0)
    qg.case='uniform'
    #
    qg.set_q()
    qg.invert_pv()
    write_nc([qg.PSI, qg.Q], ['psi', 'q'], 'output.nc', qg)
    #
    test=2
    if test==0:
        # one time step and store
        qg.tstep(1)
        write_nc([qg.PSI, qg.Q], ['psi', 'q'], 'output.nc', qg, create=False)
    elif test==1:
        # write/read/write
        qg.tstep(1)
        write_nc([qg.PSI, qg.Q], ['psi', 'q'], 'output1.nc', qg, create=True)
        qg.set_q(file_q='output.nc')
        qg.tstep(1)
        write_nc([qg.PSI, qg.Q], ['psi', 'q'], 'output1.nc', qg, create=False)
    else:
        while qg.tstepper.t/86400. < 200 :
            qg.tstep(1)
            write_nc([qg.PSI, qg.Q], ['psi', 'q'], 'output.nc', qg, create=False)
    
    return qg


def curvilinear_runs():
    ''' Tests with curvilinear grid
    ''' 
    
    qg = qg_model(hgrid = 'curv_metrics.nc', vgrid = 'curv_metrics.nc',
                  f0N2_file = 'curv_pv.nc',
                  K = 1.e3, dt = 0.5*86400.e0)
    qg.case='curv'
    #
    qg.set_q(file_q='curv_pv.nc')
    qg.invert_pv()
    write_nc([qg.PSI, qg.Q], ['psi', 'q'], 'output.nc', qg)
    
    test=1
    if test==0:
        # one time step and store
        qg.tstep(1)
        write_nc([qg.PSI, qg.Q], ['psi', 'q'], 'output.nc', qg, create=False)
    elif test==1:
        while qg.tstepper.t/86400. < 200 :
            qg.tstep(1)
            write_nc([qg.PSI, qg.Q], ['psi', 'q'], 'output.nc', qg, create=False)

    return qg


def roms_input_runs():
    ''' Tests with curvilinear grid
    '''

    hgrid = {'Lx':(512-1)*2.e3, 'Ly':(1440-1)*2.e3, 'H':4.e3, \
             'Nx':512, 'Ny':1440, 'Nz':100}
    vgrid = 'data/jet_cfg1_wp5_2km_k1e7_TSUP5_2000a3000j_zlvl_pv.nc'
    qg = qg_model(hgrid = hgrid, vgrid = vgrid, f0N2_file = vgrid, K = 1.e0, dt = 0.5*86400.e0)
    qg.case='roms'

    qg.set_q(file_q = vgrid)
    qg.set_psi(file_psi = vgrid)
    qg.set_rho(file_rho = vgrid)
    qg.invert_pv()
    write_nc([qg.PSI, qg.Q], ['psi', 'q'], 'output.nc', qg)

    return qg


if __name__ == "__main__":
    
    #qg = uniform_grid_runs()
    
    qg = roms_input_runs()
    

    if qg._verbose:
        print 'Test done \n'
