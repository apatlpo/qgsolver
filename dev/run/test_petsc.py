#!/usr/bin/python
# -*- encoding: utf8 -*-

"""
Test petsc4py basic
"""

import sys

sys.path.append('../../')
from qgsolver.qg import qg_model

#
#==================== analytical / uniform grid case =========================================
#

def uniform_grid_runs(ncores_x=16, ncores_y=16, ping_mpi_cfg=False):
    """
    Tests with uniform grid, closed domains
    """

    #
    ncores_x=2; ncores_y=2
    hgrid = {'Nx':256, 'Ny':256}
    vgrid = {'Nz':5}

    if ping_mpi_cfg:
        # escape before computing
        return ncores_x, ncores_y
    else:

        # proceeds with computations
        qg = qg_model(hgrid = hgrid, vgrid = vgrid, boundary_types={'periodic': None}, 
                      K = 0.e0, dt = 0.5*86400.e0,
                      ncores_x=ncores_x, ncores_y=ncores_y, verbose=1)
    
        # init q and psi
        qg.set_q()
        qg.state.PSI = qg.state.Q*2.
        qg.write_state(filename='data/output.nc') # 0
        # 1
        PSI = qg.state.PSI
        PSI.set(4.e-6) # correctly set qg.state.PSI to 1.e-6
        qg.write_state(filename='data/output.nc', append=True) # 1
        # 2
        Q = qg.state.Q
        Q = Q*2. # qq.state.Q is un changed, a new vector has been created
        qg.write_state(filename='data/output.nc', append=True) # 2
        # 3
        Q = qg.state.Q
        qg.state.Q = Q*2. # correctly multiply qg.state.Q by 2, Q holds old qg.state.Q
        qg.write_state(filename='data/output.nc', append=True) # 3
        # 4
        qg.state.Q = Q # correctly reset qg.state.Q to its old value
        qg.write_state(filename='data/output.nc', append=True) # 4
        # 5
        qg.state.Q = qg.state.Q/2. # qg.state.Q is divided by 2 and becomes a different object
        qg.write_state(filename='data/output.nc', append=True) # 5
        # 6
        qg.state.Q = Q # verifies that Q holds old value (2,4)
        qg.write_state(filename='data/output.nc', append=True) # 6

        return qg

def test_scope(state,q):
    pass

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
    
