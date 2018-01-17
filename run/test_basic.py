#!/usr/bin/python
# -*- encoding: utf8 -*-

"""
Test the basic features of the library:
Setup of uniform grid
PV inversion of an analytical PV distribution
"""

import time
import sys

#try:
#    from qgsolver.qg import qg_model
#    from qgsolver.io import write_nc
#except:
#    print 'qgsolver not yet in path'
#    sys.exit

from qgsolver.qg import qg_model
from qgsolver.io import write_nc


#
#==================== Synthetic curvilinear case ============================================
#



def curvilinear_runs(ncores_x=8, ncores_y=8, ping_mpi_cfg=False):
    ''' Tests with curvilinear grid
    ''' 
    
    if ping_mpi_cfg:
        # escape before computing
        return ncores_x, ncores_y
    
    else:
        # proceeds with computations            
    
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


#
#==================== NEMO case ============================================
#

def nemo_input_runs(ncores_x=2, ncores_y=6, ping_mpi_cfg=False):
    ''' Tests with curvilinear grid
    '''


    # LMX domain: Nx=1032, Ny=756, Nz=300

    # vertical subdomain
    # vdom = {'kdown': 0, 'kup': 50-1, 'k0': 200 }    # linux
    vdom = {'kdown': 0, 'kup': 50-1, 'k0': 98 }     # linux with mask
    # vdom = {'kdown': 0, 'kup': 100-1, 'k0': 115 }     # Datarmor
    # horizontal subdomain
    # hdom = {'istart': 0, 'iend': 100-1, 'i0': 450,'jstart': 0, 'jend': 100-1,  'j0': 300}   # linux
    hdom = {'istart': 0, 'iend': 100-1, 'i0': 410,'jstart': 0, 'jend': 100-1,  'j0': 590}   # linux with mask
    # hdom = {'istart': 0, 'iend': 270-1, 'i0': 135,'jstart': 0, 'jend': 384-1,  'j0': 165}     # medium datarmor
    # hdom = {'istart': 0, 'iend': 672-1, 'i0': 230,'jstart': 0, 'jend': 256-1,  'j0': 200}   # large datarmor
    # 448=8x56
    # 512=8x64
    
    # set tiling
    ncores_x=2; ncores_y=4 #   linux
    # ncores_x=2; ncores_y=8   #   medium datarmor
    # ncores_x=21; ncores_y=8  # large datarmor

    
    if ping_mpi_cfg:
        # escape before computing
        return ncores_x, ncores_y
    
    else:
        # proceeds with computations
        start_time = time.time()
        cur_time = start_time
    
        # MPI decomposition of the domain
        # case must be defined before ncores for run_caparmor.py
        casename = 'nemo'
    
        # Top and Bottom boundary condition type: 'N' for Neumann, 'D' for Dirichlet
        bdy_type = {'top':'N', 'bottom':'N'}
    
        datapath = 'data/'
        # datapath = '/home7/pharos/othr/NATL60/DIAG_DIMUP/qgsolver/mode2_fcTrue_fvertTrue/'
        # datapath = '/home7/pharos/othr/NATL60/DIAG_DIMUP/qgsolver/mode2_fcTrue_fvertTrue_new/'
        hgrid = datapath+'nemo_metrics.nc'
        vgrid = datapath+'nemo_metrics.nc'
        file_q = datapath+'nemo_pv.nc'
        file_psi = datapath+'nemo_psi.nc'
        file_rho = datapath+'nemo_rho.nc'
        qg = qg_model(hgrid=hgrid, vgrid=vgrid, f0N2_file=file_q, K=1.e0, dt=0.5 * 86400.e0,
                      vdom=vdom, hdom=hdom, ncores_x=ncores_x, ncores_y=ncores_y, 
                      bdy_type_in=bdy_type, substract_fprime=True)
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
   
        write_nc([qg.PSI, qg.Q], ['psi', 'q'], 'data/input.nc', qg)

        if qg.rank == 0: print '----------------------------------------------------'
        if qg.rank == 0: print 'Elapsed time for write_nc ', str(time.time() - cur_time)
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



def main(ping_mpi_cfg=False):    
    
    # qg = uniform_grid_runs(ping_mpi_cfg=ping_mpi_cfg)
    #
    #qg = roms_input_runs(ping_mpi_cfg=ping_mpi_cfg)
    #
    qg = nemo_input_runs(ping_mpi_cfg=ping_mpi_cfg)
    #
    # qg = test_L(ping_mpi_cfg=ping_mpi_cfg)
    
    if ping_mpi_cfg:
        return qg[0], qg[1]
    elif qg._verbose:
        print 'Test done \n'


if __name__ == "__main__":
    main()
    
    
