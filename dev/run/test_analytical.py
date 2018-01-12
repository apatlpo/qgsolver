#!/usr/bin/python
# -*- encoding: utf8 -*-

"""
Test the basic features of the library:
Setup of uniform grid
PV inversion of an analytical PV distribution
Time stepping

mpirun -n 4 python test_analytical.py
"""

import time
import sys

sys.path.append('../../')
from qgsolver.qg import qg_model
from qgsolver.state import add

#
#==================== analytical / uniform grid case =========================================
#

def uniform_grid_runs(ncores_x=16, ncores_y=16, ping_mpi_cfg=False):
    """
    Tests with uniform grid, closed domains
    """
    
    start_time = time.time()
    cur_time = start_time

    # grid
    # nhoes: 512³,  512 procs, 1000³, 4096 procs 
    # LMX case: Nx=1032=2³x3x43, Ny=756=2²x3³x7, Nz=300
    #hgrid = {'Nx':1032, 'Ny':756}
    #vgrid = {'Nz':300}
    #
    
    #
    ncores_x=2; ncores_y=2
    hgrid = {'Nx':256, 'Ny':256}
    vgrid = {'Nz':5}
    
    # :
    #ncores_x=16; ncores_y=16
    #hgrid = {'Nx':512, 'Ny':512}
    #vgrid = {'Nz':100}       
    #vgrid = {'Nz':200}
    
    # requires 16x8:
    #ncores_x=16; ncores_y=8
    #hgrid = {'Nx':512, 'Ny':512}
    #vgrid = {'Nz':300}    
    
    # requires 16x16
    #ncores_x=16; ncores_y=16
    #hgrid = {'Nx':1024, 'Ny':512}
    #vgrid = {'Nz':300}

    if ping_mpi_cfg:
        # escape before computing
        return ncores_x, ncores_y
    else:

        # proceeds with computations
        qg = qg_model(hgrid = hgrid, vgrid = vgrid, boundary_types={'periodic': True}, 
                      K = 0.e0, dt = 0.5*86400.e0,
                      ncores_x=ncores_x, ncores_y=ncores_y, verbose=1)

        # pv inversion
        qg.set_q()
        qg.invert_pv()
        qg.write_state(filename='data/output.nc')
        #
        if True:
            bstate = qg.set_bstate(psi0=0., q0=0., beta=1.e-11, rho0=0.)
            #
            # bstate=None # turns off background state
            if False:
                # add the background state to qg.state, debug
                add(qg.state,bstate,da=None)
                qg.write_state(filename='data/output.nc', append=True)
                #
                qg.invert_pv(bstate=bstate)
                #qg.invert_pv(bstate=bstate, addback_bstate=False) # test
        
        #
        test=1
        if test==0:
            # one time step and store
            qg.tstep(1, bstate=bstate)
            qg.write_state(filename='data/output.nc', append=True)
        elif test==1:
            while qg.tstepper.t/86400. < 200 :
                qg.tstep(2)
                qg.write_state(filename='data/output.nc', append=True)
    
        if qg._verbose>0:
            print('----------------------------------------------------')
            print('Elapsed time for all ',str(time.time() - cur_time))
    
        return qg

#
#==================== main wrappers =========================================
#


def main(ping_mpi_cfg=False):    
    
    qg = uniform_grid_runs(ping_mpi_cfg=ping_mpi_cfg)
    
    if ping_mpi_cfg:
        return qg[0], qg[1]
    elif qg._verbose>0:
        print('Test analytical done \n')

if __name__ == "__main__":
    main()
    
