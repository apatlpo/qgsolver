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


def main(ping_mpi_cfg=False):    
    
    qg = uniform_grid_runs(ping_mpi_cfg=ping_mpi_cfg)
    
    if ping_mpi_cfg:
        return qg[0], qg[1]
    elif qg._verbose:
        print 'Test done \n'


if __name__ == "__main__":
    main()
    
    
