#!/usr/bin/python
# -*- encoding: utf8 -*-

"""
Test the basic features of the library:
Setup of uniform grid
PV inversion of an analytical PV distribution
"""

import time
import sys
from qgsolver.qg import qg_model
from qgsolver.io import write_nc


def uniform_grid_runs():
    """
    Tests with uniform grid, closed domains
    """
    qg = qg_model(hgrid = {'Nx0':150, 'Ny0':100}, vgrid = {'Nz0':3 },
            K = 0.e0, dt = 0.5*86400.e0)
    qg.case='uniform'
    #
    qg.set_q()
    qg.invert_pv()
    write_nc([qg.PSI, qg.Q], ['psi', 'q'], 'data/output.nc', qg)
    #
    test=2
    if test==0:
        # one time step and store
        qg.tstep(1)
        write_nc([qg.PSI, qg.Q], ['psi', 'q'], 'data/output.nc', qg, create=False)
    elif test==1:
        # write/read/write
        qg.tstep(1)
        write_nc([qg.PSI, qg.Q], ['psi', 'q'], 'data/output1.nc', qg, create=True)
        qg.set_q(file_q='data/output.nc')
        qg.tstep(1)
        write_nc([qg.PSI, qg.Q], ['psi', 'q'], 'data/output1.nc', qg, create=False)
    else:
        while qg.tstepper.t/86400. < 200 :
            qg.tstep(1)
            write_nc([qg.PSI, qg.Q], ['psi', 'q'], 'data/output.nc', qg, create=False)
    
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
    write_nc([qg.PSI, qg.Q], ['psi', 'q'], 'data/output.nc', qg)
    
    test=1
    if test==0:
        # one time step and store
        qg.tstep(1)
        write_nc([qg.PSI, qg.Q], ['psi', 'q'], 'data/output.nc', qg, create=False)
    elif test==1:
        while qg.tstepper.t/86400. < 200 :
            qg.tstep(1)
            write_nc([qg.PSI, qg.Q], ['psi', 'q'], 'data/output.nc', qg, create=False)

    return qg


def roms_input_runs():
    ''' Tests with curvilinear grid
    '''

    start_time = time.time()
    cur_time = start_time

    # MPI decomposition of the domain
    # case must be defined before ncores for run_caparmor.py
    casename='roms'
    ncores_x= 2
    ncores_y=4
    
    # vertical subdomain
    # vdom = {'kdown': 0, 'kup': 49, 'k0': 0 }
    vdom = {'kdown': 25, 'kup': 35, 'k0': 15 }

    # horizontal subdomain
    # hdom = {'istart': 0, 'iend': 255, 'i0': 0, 'jstart':0, 'jend':721,  'j0': 0}
    hdom = {'istart': 50, 'iend': 200, 'i0': 40, 'jstart':100, 'jend':600,  'j0': 90}


    # hgrid = {'Lx':(512-1)*2.e3, 'Ly':(1440-1)*2.e3, 'H':4.e3, \
    #          'Nx':512, 'Ny':1440, 'Nz':100}
    hgrid = {'Lx':(256-1)*4.e3, 'Ly':(720-1)*4.e3, 'Nx0':256, 'Ny0':722}
    # vgrid = 'data/jet_cfg1_wp5_2km_k1e7_TSUP5_2000a3000j_zlvl_pv.nc'
    vgrid = 'data/jet_cfg1_wp5_4km_k3.2e8_0a1500j_zlvl_pv.nc'
    qg = qg_model(hgrid = hgrid, vgrid = vgrid, f0N2_file = vgrid, K = 1.e0, dt = 0.5*86400.e0, 
                  vdom=vdom, hdom=hdom, ncores_x=ncores_x, ncores_y=ncores_y)
    qg.case=casename


    if qg.rank == 0: print '----------------------------------------------------'
    if qg.rank == 0: print 'Elapsed time for qg_model ',str(time.time() - cur_time)
    cur_time = time.time()

    qg.set_q(file_q = vgrid)
    if qg.rank == 0: print '----------------------------------------------------'
    if qg.rank == 0: print 'Elapsed time for set_q ',str(time.time() - cur_time)
    cur_time = time.time()

    qg.set_psi(file_psi = vgrid)
    if qg.rank == 0: print '----------------------------------------------------'
    if qg.rank == 0: print 'Elapsed time for set_psi ',str(time.time() - cur_time)
    cur_time = time.time()

    qg.set_rho(file_rho = vgrid)
    if qg.rank == 0: print '----------------------------------------------------'
    if qg.rank == 0: print 'Elapsed time for set_rho ',str(time.time() - cur_time)
    cur_time = time.time()

    qg.invert_pv()
    if qg.rank == 0: print '----------------------------------------------------'
    if qg.rank == 0: print 'Elapsed time for invert_pv ',str(time.time() - cur_time)
    cur_time = time.time()

    write_nc([qg.PSI, qg.Q], ['psi', 'q'], 'data/output.nc', qg)
    if qg.rank == 0: print '----------------------------------------------------'
    if qg.rank == 0: print 'Elapsed time for write_nc ',str(time.time() - cur_time)
    cur_time = time.time()

    if qg.rank == 0: print '----------------------------------------------------'
    if qg.rank == 0: print 'Elapsed time  ',str(cur_time - start_time)

    return qg


def nemo_input_runs():
    ''' Tests with curvilinear grid
    '''

    start_time = time.time()
    cur_time = start_time

    # MPI decomposition of the domain
    # case must be defined before ncores for run_caparmor.py
    casename = 'nemo'
    ncores_x = 2
    ncores_y = 4

    # vertical subdomain
    vdom = {'kdown': 150, 'kup': 250, 'k0': 150 }

    # horizontal subdomain
    hdom = {'istart': 230, 'iend': 450, 'i0': 230,'jstart': 250, 'jend': 500,  'j0': 250 }


    hgrid = 'data/nemo_metrics.nc'
    vgrid = 'data/nemo_metrics.nc'
    file_q = 'data/nemo_pv.nc'
    file_psi = 'data/nemo_psi.nc'
    file_rho = 'data/nemo_rho.nc'
    qg = qg_model(hgrid=hgrid, vgrid=vgrid, f0N2_file=file_q, K=1.e0, dt=0.5 * 86400.e0,
                  vdom=vdom, hdom=hdom, ncores_x=ncores_x, ncores_y=ncores_y)
    qg.case=casename

    if qg.rank == 0: print '----------------------------------------------------'
    if qg.rank == 0: print 'Elapsed time for qg_model ', str(time.time() - cur_time)
    cur_time = time.time()

    qg.set_q(file_q=file_q)
    if qg.rank == 0: print '----------------------------------------------------'
    if qg.rank == 0: print 'Elapsed time for set_q ', str(time.time() - cur_time)
    cur_time = time.time()

    qg.set_psi(file_psi=file_psi)
    if qg.rank == 0: print '----------------------------------------------------'
    if qg.rank == 0: print 'Elapsed time for set_psi ', str(time.time() - cur_time)
    cur_time = time.time()

    qg.set_rho(file_rho=file_rho)
    if qg.rank == 0: print '----------------------------------------------------'
    if qg.rank == 0: print 'Elapsed time for set_rho ', str(time.time() - cur_time)
    cur_time = time.time()

    qg.invert_pv()
    if qg.rank == 0: print '----------------------------------------------------'
    if qg.rank == 0: print 'Elapsed time for invert_pv ', str(time.time() - cur_time)
    cur_time = time.time()

    write_nc([qg.PSI, qg.Q], ['psi', 'q'], 'data/output.nc', qg)
    if qg.rank == 0: print '----------------------------------------------------'
    if qg.rank == 0: print 'Elapsed time for write_nc ', str(time.time() - cur_time)
    cur_time = time.time()

    if qg.rank == 0: print '----------------------------------------------------'
    if qg.rank == 0: print 'Elapsed time  ', str(cur_time - start_time)

    return qg

def test_L():

    start_time = time.time()
    cur_time = start_time

    # hgrid = {'Lx':(512-1)*2.e3, 'Ly':(1440-1)*2.e3, 'H':4.e3, \
    #          'Nx':512, 'Ny':1440, 'Nz':100}
    hgrid = {'Lx':(256-1)*4.e3, 'Ly':(720-1)*4.e3, 'H':4.e3,
             'Nx':256, 'Ny':720, 'Nz':50}
    # vgrid = 'data/jet_cfg1_wp5_2km_k1e7_TSUP5_2000a3000j_zlvl_pv.nc'
    vgrid = 'data/jet_cfg1_wp5_4km_k3.2e8_0a1500j_zlvl_pv.nc'
    qg = qg_model(hgrid = hgrid, vgrid = vgrid, f0N2_file = vgrid, K = 1.e0, dt = 0.5*86400.e0)
    qg.case='roms'
    if qg.rank == 0: print '----------------------------------------------------'
    if qg.rank == 0: print 'Elapsed time for qg_model ',str(time.time() - cur_time)
    cur_time = time.time()

    qg.set_psi(file_psi = vgrid)
    if qg.rank == 0: print '----------------------------------------------------'
    if qg.rank == 0: print 'Elapsed time for set_psi ',str(time.time() - cur_time)
    cur_time = time.time()

    qg.pvinv.L.mult(qg.PSI,qg.Q)
    qg.invert_pv()
    if qg.rank == 0: print '----------------------------------------------------'
    if qg.rank == 0: print 'Elapsed time for invert_pv ',str(time.time() - cur_time)

    write_nc([qg.PSI, qg.Q], ['psi', 'q'], 'data/Lpsi_invPV.nc', qg)

if __name__ == "__main__":
    
    #qg = uniform_grid_runs()
    
    # qg = roms_input_runs()
    qg = nemo_input_runs()
    # qg = test_L()
    

    if qg._verbose:
        print 'Test done \n'
