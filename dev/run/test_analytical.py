#!/usr/bin/python
# -*- encoding: utf8 -*-

"""
Test the basic features of the library:
Setup of uniform grid
PV inversion of an analytical PV distribution
Time stepping
"""

import time
import sys

sys.path.append('../../')
from qgsolver.qg import qg_model
from qgsolver.inout import write_nc


#
#==================== Uniform grid case ============================================
#

def uniform_grid_runs():
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

    # proceeds with computations
    #qg = qg_model(hgrid = hgrid, vgrid = vgrid,
    #              K = 0.e0, dt = 0.5*86400.e0,
    #              ncores_x=ncores_x, ncores_y=ncores_y)
    qg = qg_model(hgrid = hgrid, vgrid = vgrid,
                  K = 0.e0, dt = None,
                  ncores_x=ncores_x, ncores_y=ncores_y)
    #
    qg.set_q()
    qg.invert_pv()
    write_nc([qg['PSI'], qg['Q']], ['psi', 'q'], 'data/output.nc', qg)

    # load background PV

    #
    test=-1
    if test==0:
        # one time step and store
        qg.tstep(1)
        write_nc([qg['PSI'], qg['Q']], ['psi', 'q'], 'data/output.nc', qg, create=False)
    elif test==1:
        # write/read/write
        qg.tstep(1)
        write_nc([qg['PSI'], qg['Q']], ['psi', 'q'], 'data/output.nc', qg, create=False)
        #qg.set_q(file_q='data/output.nc')
        qg.tstep(1)
        write_nc([qg['PSI'], qg['Q']], ['psi', 'q'], 'data/output.nc', qg, create=False)
    elif test==2:
        while qg.tstepper.t/86400. < 200 :
            qg.tstep(1)
            write_nc([qg.state.PSI, qg.state.Q], ['psi', 'q'], 'data/output.nc', qg, create=False)

    if qg._verbose>0:
        print('----------------------------------------------------')
        print('Elapsed time for all ',str(time.time() - cur_time))

    return qg


if __name__ == "__main__":

    qg = uniform_grid_runs()

    if qg._verbose>0:
        print('Test done \n')

    
