#!/usr/bin/python
# -*- encoding: utf8 -*-

"""
Test the basic features of the library:
Setup of uniform grid
PV inversion of an analytical PV distribution
Time stepping

mpirun -n 4 python uniform.py
or
mpirun -n 4 python uniform.py -mf -ksp_view -ksp_monitor -ksp_converged_reason
"""

import time
import sys

sys.path.append('../')
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

    if ping_mpi_cfg:
        # escape before computing
        return ncores_x, ncores_y
    else:

        # proceeds with computations
        qg = qg_model(hgrid = hgrid, vgrid = vgrid, boundary_types={'periodic': True},
                      ncores_x=ncores_x, ncores_y=ncores_y, verbose=1)

        # pv inversion
        qg.set_q()
        PSI = qg.state.PSI
        Q = qg.state.Q
        qg.write_state(filename='data/output.nc') # 1

        # axpy(self, alpha, Vec x)
        PSI.axpy(2.,Q)   # PSI = PSI + 2*Q = 2Q
        qg.write_state(filename='data/output.nc',append=True) # 2
        PSI.axpy(2.,Q)   # PSI = PSI + 2*Q = 4*Q
        qg.write_state(filename='data/output.nc',append=True) # 3
        # Computes y = alpha x + y
        # http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Vec/VecAXPY.html

        # waxpy(self, alpha, Vec x, Vec y)
        PSI.waxpy(2.,Q,Q)   # PSI = 2*Q + Q = 3*Q
        qg.write_state(filename='data/output.nc',append=True) #4
        # Computes w = alpha x + y.
        # http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Vec/VecWAXPY.html#VecWAXPY

        # copy
        Q.copy(PSI)   # PSI = Q  !! i.e. reverse assignement
        qg.write_state(filename='data/output.nc',append=True) #4



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
        print('Test vector oprations done \n')

if __name__ == "__main__":
    main()

