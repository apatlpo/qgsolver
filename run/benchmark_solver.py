#!/usr/bin/python
# -*- encoding: utf8 -*-

"""
Test the basic features of the library:
Setup of uniform grid
PV inversion of an analytical PV distribution
Time stepping

mpirun -n 4 python benchmark_solver.py
or
mpirun -n 4 python benchmark_solver.py -mf -ksp_view -ksp_monitor -ksp_converged_reason
"""

import time
import sys

sys.path.append('../')
from qgsolver.qg import qg_model
from qgsolver.state import add


from timeit import default_timer as timer

class benchmark(object):

    def __enter__(self):
        self.start = timer()
        return self

    def __exit__(self, *args):
        self.delt = timer() - self.start

#
#==================== analytical / uniform grid case =========================================
#

def uniform_grid_runs(ncores_x=16, ncores_y=16, ping_mpi_cfg=False):
    """
    Tests with uniform grid, closed domains
    """
    
    start_time = time.time()
    cur_time = start_time

    #
    ncores_x=2; ncores_y=2
    hgrid = {'Nx':128, 'Ny':128}
    vgrid = {'Nz':5}

    if ping_mpi_cfg:
        # escape before computing
        return ncores_x, ncores_y
    else:

        with benchmark() as b:

            # proceeds with computations
            qg = qg_model(hgrid = hgrid, vgrid = vgrid, boundary_types={'periodic': True},
                          K = 0.e0, dt = 0.5*86400.e0,
                          ncores_x=ncores_x, ncores_y=ncores_y, verbose=1,
                          solver='lgmres')
            # default is right preconditioning for: fgmres, qcg
            # default is left  preconditioning for: gmres, bicg
            # bcgsl
            # restart: can you see with one time step effect?
        #
        delt_init = b.delt

        # pv inversion in order to compute psi
        qg.set_q()
        qg.invert_pv()
        #
        bstate = qg.set_bstate(psi0=0., q0=0., beta=1.e-11)

        with benchmark() as b:
            #
            # 100 time step and store
            #qg.tstep(5, bstate=bstate)
            pass
        #
        delt_compute = b.delt

        if qg._verbose>0:
            print('----------------------------------------------------')
            print('Elapsed time for all %.1f s \n' %(time.time() - cur_time))

        verbose = qg._verbose
        del qg

        return delt_init, delt_compute, verbose

#
#==================== main wrappers =========================================
#


def main(ping_mpi_cfg=False):    

    #
    delt_init, delt_compute, verbose = uniform_grid_runs(ping_mpi_cfg=ping_mpi_cfg)

    if not ping_mpi_cfg and verbose:
        print('delt_init = %.1f s' %delt_init)
        print('delt_compute = %.1f s\n' %delt_compute)
        print('Benchmark test done\n')

if __name__ == "__main__":
    main()

