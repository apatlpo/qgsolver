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
#==================== Uniform case ============================================
#

def uniform_grid_runs(ncores_x=16, ncores_y=16, ping_mpi_cfg=False):
    """
    Tests with uniform grid, closed domains
    """
    
    start_time = time.time()
    cur_time = start_time
    
    # MPI decomposition of the domain
    # case must be defined before ncores for run_caparmor.py
    casename='uniform'
    
    # grid
    # nhoes: 512³,  512 procs, 1000³, 4096 procs 
    # LMX case: Nx=1032=2³x3x43, Ny=756=2²x3³x7, Nz=300
    #hgrid = {'Nx':1032, 'Ny':756}
    #vgrid = {'Nz':300}
    #
    
    #
    #ncores_x=8; ncores_y=2
    #hgrid = {'Nx':256, 'Ny':256}
    #vgrid = {'Nz':50}           
    
    # :
    ncores_x=16; ncores_y=16
    hgrid = {'Nx':512, 'Ny':512}
    #vgrid = {'Nz':100}       
    vgrid = {'Nz':200}
    
    # requires 16x8:
    #ncores_x=16; ncores_y=8
    #hgrid = {'Nx':512, 'Ny':512}
    #vgrid = {'Nz':300}    
    
    # requires 16x16
    #ncores_x=16; ncores_y=16
    #hgrid = {'Nx':1024, 'Ny':512}
    #vgrid = {'Nz':300}
    # crashes with io
    
    # no io
    # 512 x 512 x 100 on 8x8: 160 s, 136 iter
    # 512 x 512 x 300 on 8x8: crash, out of memory
    #     -------     on 16x8: 237 s, 144 iter
    # 1024 x 512 x 300 on 16 x 8: crash, out of memory
    #     -------     on 16x16: 379s s, 236 iter

    
    
    if ping_mpi_cfg:
        # escape before computing
        return ncores_x, ncores_y
    
    else:
        # proceeds with computations        
        qg = qg_model(hgrid = hgrid, vgrid = vgrid,
                K = 0.e0, dt = 0.5*86400.e0,
                ncores_x=ncores_x, ncores_y=ncores_y)
        qg.case='uniform'
        #
        qg.set_q()
        qg.set_rho()
        qg.set_psi()    
        qg.invert_pv()
        write_nc([qg.PSI, qg.Q], ['psi', 'q'], 'data/output.nc', qg)
        
        if qg._verbose>0:
            print '----------------------------------------------------'
            print 'Elapsed time for all ',str(time.time() - cur_time)
        
        #
        test=-1
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
        elif test==2:
            while qg.tstepper.t/86400. < 200 :
                qg.tstep(1)
                write_nc([qg.PSI, qg.Q], ['psi', 'q'], 'data/output.nc', qg, create=False)
        
        return qg








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
#==================== ROMS case ============================================
#



def roms_input_runs(ncores_x=2, ncores_y=4, ping_mpi_cfg=False):
    ''' Tests with roms configuration (spatially uniform, vertically stretched)
    '''

    if ping_mpi_cfg:
        # escape before computing
        return ncores_x, ncores_y
    
    else:
        # proceeds with computations
    
        start_time = time.time()
        cur_time = start_time
    
        # MPI decomposition of the domain
        # case must be defined before ncores for run_caparmor.py
        casename='roms'
        
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





#
#==================== NEMO case ============================================
#

def nemo_input_runs(ncores_x=2, ncores_y=6, ping_mpi_cfg=False):
    ''' Tests with curvilinear grid
    '''


    # LMX domain: Nx=1032, Ny=756, Nz=300

    # vertical subdomain
    vdom = {'kdown': 0, 'kup': 100-1, 'k0': 115 }

    # horizontal subdomain
    hdom = {'istart': 0, 'iend': 270-1, 'i0': 135,'jstart': 0, 'jend': 384-1,  'j0': 165}
    #hdom = {'istart': 0, 'iend': 448-1, 'i0': 230,'jstart': 0, 'jend': 256-1,  'j0': 200}
    # 448=8x56
    # 512=8x64
    
    # set tiling
    ncores_x=2; ncores_y=16
    ncores_x=6; ncores_y=16

    
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
    
        hgrid = 'data/nemo_metrics.nc'
        vgrid = 'data/nemo_metrics.nc'
        file_q = 'data/nemo_pv.nc'
        file_psi = 'data/nemo_psi.nc'
        file_rho = 'data/nemo_rho.nc'
        qg = qg_model(hgrid=hgrid, vgrid=vgrid, f0N2_file=file_q, K=1.e0, dt=0.5 * 86400.e0,
                      vdom=vdom, hdom=hdom, ncores_x=ncores_x, ncores_y=ncores_y, 
                      bdy_type=bdy_type, substract_fprime=True)
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

        #write_nc(qg.PSI, 'psi', 'data/input.nc', qg)
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





def test_L(ncores_x=2, ncores_y=4, ping_mpi_cfg=False):
    
    if ping_mpi_cfg:
        # escape before computing
        return ncores_x, ncores_y
    
    else:
        # proceeds with computations
        
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




def main(ping_mpi_cfg=False):    
    
    qg = uniform_grid_runs(ping_mpi_cfg=ping_mpi_cfg)
    #
    #qg = roms_input_runs(ping_mpi_cfg=ping_mpi_cfg)
    #
    #qg = nemo_input_runs(ping_mpi_cfg=ping_mpi_cfg)
    #
    # qg = test_L(ping_mpi_cfg=ping_mpi_cfg)
    
    if ping_mpi_cfg:
        return qg[0], qg[1]
    elif qg._verbose:
        print 'Test done \n'


if __name__ == "__main__":
    main()
    
    
