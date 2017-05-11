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
        
        # Top and Bottom boundary condition type: 'N' for Neumann, 'D' for Dirichlet
        #bdy_type = {'top':'N', 'bottom':'N'}
        bdy_type = {'top':'N', 'bottom':'N', 'periodic':True}

        # vertical subdomain
        #vdom = {'kdown': 0, 'kup': 49, 'k0': 0 }
        #vdom = {'kdown': 0, 'kup': 30, 'k0': 0 }
        #vdom = {'Nz': 30}
        vdom = {'Nz': 50}

        # horizontal subdomain
        # hdom = {'istart': 0, 'iend': 255, 'i0': 0, 'jstart':0, 'jend':721,  'j0': 0}
        hdom = {'Nx': 256, 'j0':300, 'Ny':200}
        #hdom = {'Nx': 256, 'Ny': 720}


        datapath = 'input/'
        hgrid = datapath+'roms_metrics.nc'
        vgrid = datapath+'roms_metrics.nc'
        file_q = datapath+'roms_pv.nc'
        file_psi = datapath+'roms_psi.nc'
        file_rho = datapath+'roms_rho.nc'
        qg = qg_model(hgrid=hgrid, vgrid=vgrid, f0N2_file=file_q, K=1.e0, dt=0.5 * 86400.e0,
                      vdom=vdom, hdom=hdom, ncores_x=ncores_x, ncores_y=ncores_y, 
                      bdy_type_in=bdy_type, substract_fprime=True, verbose=1)        
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
   
        write_nc([qg.PSI, qg.Q], ['psi', 'q'], 'output/input.nc', qg)

        if qg.rank == 0: print '----------------------------------------------------'
        if qg.rank == 0: print 'Elapsed time for write_nc ', str(time.time() - cur_time)
        cur_time = time.time()
    
        qg.invert_pv()
        if qg.rank == 0: print '----------------------------------------------------'
        if qg.rank == 0: print 'Elapsed time for invert_pv ', str(time.time() - cur_time)
        cur_time = time.time()

        write_nc([qg.PSI, qg.Q], ['psi', 'q'], 'output/output.nc', qg)
        if qg.rank == 0: print '----------------------------------------------------'
        if qg.rank == 0: print 'Elapsed time for write_nc ', str(time.time() - cur_time)
        cur_time = time.time()

        #
        test=0
        if test==0:
            # one time step and store
            qg.tstep(1)
            write_nc([qg.PSI, qg.Q], ['psi', 'q'], 'output/output.nc', qg, create=False)
        elif test==1:
            # write/read/write
            qg.tstep(1)
            write_nc([qg.PSI, qg.Q], ['psi', 'q'], 'output/output1.nc', qg, create=True)
            qg.set_q(file_q='data/output.nc')
            qg.tstep(1)
            write_nc([qg.PSI, qg.Q], ['psi', 'q'], 'output/output1.nc', qg, create=False)
        elif test==2:
            while qg.tstepper.t/86400. < 200 :
                qg.tstep(1)
                write_nc([qg.PSI, qg.Q], ['psi', 'q'], 'output/output.nc', qg, create=False)
             

        if qg.rank == 0: print '----------------------------------------------------'
        if qg.rank == 0: print 'Elapsed time  ', str(cur_time - start_time)
    
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
    
    #
    qg = roms_input_runs(ping_mpi_cfg=ping_mpi_cfg)
    #
    #qg = nemo_input_runs(ping_mpi_cfg=ping_mpi_cfg)
    #
    
    if ping_mpi_cfg:
        return qg[0], qg[1]
    elif qg._verbose:
        print 'Test done \n'


if __name__ == "__main__":
    main()
    
    
