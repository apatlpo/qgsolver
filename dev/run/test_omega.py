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
from qgsolver.inout import write_nc

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
        write_nc([qg.PSI, qg.Q], ['psi', 'q'], '../output/output.nc', qg)
        
        if qg._verbose>0:
            print '----------------------------------------------------'
            print 'Elapsed time for all ',str(time.time() - cur_time)
        
        #
        test=-1
        if test==0:
            # one time step and store
            qg.tstep(1)
            write_nc([qg.PSI, qg.Q], ['psi', 'q'], '../output/output.nc', qg, create=False)
        elif test==1:
            # write/read/write
            qg.tstep(1)
            write_nc([qg.PSI, qg.Q], ['psi', 'q'], '../output/output1.nc', qg, create=True)
            qg.set_q(file_q='../output/output.nc')
            qg.tstep(1)
            write_nc([qg.PSI, qg.Q], ['psi', 'q'], '../output/output1.nc', qg, create=False)
        elif test==2:
            while qg.tstepper.t/86400. < 200 :
                qg.tstep(1)
                write_nc([qg.PSI, qg.Q], ['psi', 'q'], '../output/output.nc', qg, create=False)
        
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
        
        # Top and Bottom boundary condition type: 'N' for Neumann, 'D' for Dirichlet
        bdy_type = {'top':'D', 'bottom':'D', 'periodic':True}
        #bdy_type = {'top':'D', 'bottom':'D'}

        # vertical subdomain
        #vdom = {'kdown': 0, 'kup': 49, 'k0': 0 }
        # vdom = {'kdown': 0, 'kup': 30, 'k0': 10 }
        #vdom = {'Nz': 30}
        vdom = {'Nz': 50}

        # horizontal subdomain
        # hdom = {'istart': 0, 'iend': 255, 'i0': 0, 'jstart':0, 'jend':721,  'j0': 0}
        # hdom = {'Nx': 256, 'j0':300, 'Ny':200}
        hdom = {'Nx': 256, 'Ny': 720}
        # 256 = 2^8
        # 720 = 2^4 x 3^2 x 5

        datapath = '../input/'
        outdir = '../output/'
        hgrid = datapath+'roms_metrics.nc'
        vgrid = datapath+'roms_metrics.nc'
        file_psi = datapath+'roms_psi.nc'
        file_q = datapath+'roms_pv.nc'
        qg = qg_model(hgrid=hgrid, vgrid=vgrid, f0N2_file=file_q, K=20.e0,
                      vdom=vdom, hdom=hdom, ncores_x=ncores_x, ncores_y=ncores_y, 
                      bdy_type_in=bdy_type, substract_fprime=True, verbose=1, 
                      flag_pvinv=False, flag_omega=True)        
        qg.case=casename
    
        # prepare inversion
        qg.set_psi(file_psi=file_psi)
        qg.set_w()
        write_nc([qg.PSI, qg.W], ['psi', 'w'], outdir+'input.nc', qg)
    
        # invert and store
        qg.invert_omega()
        write_nc([qg.W], ['w'], outdir+'output.nc', qg)
    
        return qg





#
#==================== NEMO case ============================================
#

def nemo_input_runs(ncores_x=2, ncores_y=6, ping_mpi_cfg=False):
    ''' Tests with curvilinear grid
    '''


    # LMX domain: Nx=1032, Ny=756, Nz=300

    # vertical subdomain
    # vdom = {'kdown': 0, 'kup': 50-1, 'k0': 200 }    # small without mask
    vdom = {'kdown': 0, 'kup': 160-1, 'k0': 140 }    # medium without mask
    # vdom = {'kdown': 0, 'kup': 50-1, 'k0': 98 }     # linux with mask
    # vdom = {'kdown': 0, 'kup': 100-1, 'k0': 115 }     # Datarmor
    # horizontal subdomain
    # hdom = {'istart': 0, 'iend': 100-1, 'i0': 450,'jstart': 0, 'jend': 100-1,  'j0': 300}   # small without mask
    hdom = {'istart': 0, 'iend': 300-1, 'i0': 200,'jstart': 0, 'jend': 300-1,  'j0': 200}   # medium without mask
    # hdom = {'istart': 0, 'iend': 100-1, 'i0': 410,'jstart': 0, 'jend': 100-1,  'j0': 590}   # linux with mask
    # hdom = {'istart': 0, 'iend': 270-1, 'i0': 135,'jstart': 0, 'jend': 384-1,  'j0': 165}     # medium datarmor
    # hdom = {'istart': 0, 'iend': 672-1, 'i0': 230,'jstart': 0, 'jend': 256-1,  'j0': 200}   # large datarmor
    # 448=8x56
    # 512=8x64
    
    # set tiling
    # ncores_x=2; ncores_y=4 #   small witout mask
    ncores_x=5; ncores_y=5 #   medium without mask
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
    
        datapath = '../input/'
        # datapath = '/home7/pharos/othr/NATL60/DIAG_DIMUP/qgsolver/mode2_fcTrue_fvertTrue/'
        # datapath = '/home7/pharos/othr/NATL60/DIAG_DIMUP/qgsolver/mode2_fcTrue_fvertTrue_new/'
        hgrid = datapath+'nemo_metrics.nc'
        vgrid = datapath+'nemo_metrics.nc'
        file_q = datapath+'nemo_pv.nc'
        file_psi = datapath+'nemo_psi.nc'
        file_rho = datapath+'nemo_rho.nc'
        qg = qg_model(hgrid=hgrid, vgrid=vgrid, f0N2_file=file_q, K=1.e0, dt=0.5 * 86400.e0,
                      vdom=vdom, hdom=hdom, ncores_x=ncores_x, ncores_y=ncores_y, 
                      bdy_type_in=bdy_type, substract_fprime=True,
                      flag_pvinv=False, flag_omega=True)
        qg.case=casename
    
        if qg.rank == 0: print '----------------------------------------------------'
        if qg.rank == 0: print 'Elapsed time for qg_model ', str(time.time() - cur_time)
        cur_time = time.time()
    
        qg.set_psi(file_psi=file_psi)
        if qg.rank == 0: print '----------------------------------------------------'
        if qg.rank == 0: print 'Elapsed time for set_psi ', str(time.time() - cur_time)
        cur_time = time.time()
        
        #qg.set_w(file_psi=file_psi)
        qg.set_w()
        if qg.rank == 0: print '----------------------------------------------------'
        if qg.rank == 0: print 'Elapsed time for set_w ', str(time.time() - cur_time)
        cur_time = time.time()        
        
        write_nc([qg.PSI, qg.W], ['psi', 'w'], '../output/input.nc', qg)
        if qg.rank == 0: print '----------------------------------------------------'
        if qg.rank == 0: print 'Elapsed time for write_nc ', str(time.time() - cur_time)
        cur_time = time.time()
    
        qg.invert_omega()
        if qg.rank == 0: print '----------------------------------------------------'
        if qg.rank == 0: print 'Elapsed time for invert_pv ', str(time.time() - cur_time)
        cur_time = time.time()

        write_nc([qg.W], ['w'], '../output/output.nc', qg)
        if qg.rank == 0: print '----------------------------------------------------'
        if qg.rank == 0: print 'Elapsed time for write_nc ', str(time.time() - cur_time)
        cur_time = time.time()
    
        if qg.rank == 0: print '----------------------------------------------------'
        if qg.rank == 0: print 'Elapsed time  ', str(cur_time - start_time)
    
        return qg






def main(ping_mpi_cfg=False):    
    
    # qg = uniform_grid_runs(ping_mpi_cfg=ping_mpi_cfg)
    #
    qg = roms_input_runs(ping_mpi_cfg=ping_mpi_cfg)
    #
    #qg = nemo_input_runs(ping_mpi_cfg=ping_mpi_cfg)
    #
    # qg = test_L(ping_mpi_cfg=ping_mpi_cfg)
    
    if ping_mpi_cfg:
        return qg[0], qg[1]
    elif qg._verbose:
        print 'All done \n'


if __name__ == "__main__":
    main()
    
    
