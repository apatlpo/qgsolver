#!/usr/bin/python
# -*- encoding: utf8 -*-

"""
Create a window for spectral computations in the presence of a mask
"""

import time
import sys

from qgsolver.window import window
from qgsolver.inout import write_nc
from petsc4py import PETSc

#
#==================== Window case ============================================
#

def window_runs(ncores_x=2, ncores_y=6, ping_mpi_cfg=False):
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

    # proceeds with computations
    start_time = time.time()
    cur_time = start_time

    # MPI decomposition of the domain
    # case must be defined before ncores for run_caparmor.py
    casename = 'window'

    datapath = 'data/'
    hgrid = datapath+'window_metrics.nc'
    vgrid = datapath+'window_metrics.nc'
    win = window(hgrid=hgrid, vgrid=vgrid, K=1.e-5, mask3D=True,
                vdom=vdom, hdom=hdom, ncores_x=ncores_x, ncores_y=ncores_y)
    win.case=casename

    if win.rank == 0: print('----------------------------------------------------')
    if win.rank == 0: print('Elapsed time for window_model ', str(time.time() - cur_time))
    cur_time = time.time()

    win.set_q()
    if win.rank == 0: print('----------------------------------------------------')
    if win.rank == 0: print('Elapsed time for set_q (//RHS)', str(time.time() - cur_time))
    cur_time = time.time()

    win.invert_win()
    if win.rank == 0: print('----------------------------------------------------')
    if win.rank == 0: print('Elapsed time for invert_win ', str(time.time() - cur_time))
    cur_time = time.time()

    # normalize by its maximum value
    Wmax = win.PSI.norm(norm_type=PETSc.NormType.INF)
    if win.rank == 0: print('----------------------------------------------------')
    if win.rank == 0: print('Infinite norm of window is =%f' %Wmax
    win.PSI.scale(1./Wmax)
    if win.rank == 0: print('Window is normalized by this infinite norm')


    # output the window
    write_nc([win.PSI], ['window'], 'data/output.nc', win)
    if win.rank == 0: print('----------------------------------------------------')
    if win.rank == 0: print('Elapsed time for write_nc ', str(time.time() - cur_time))
    cur_time = time.time()

    if win.rank == 0: print('----------------------------------------------------')
    if win.rank == 0: print('Elapsed time  ', str(cur_time - start_time))

    return win


def main(ping_mpi_cfg=False):

    win = window_runs(ping_mpi_cfg=ping_mpi_cfg)

    if ping_mpi_cfg:
        return win[0], win[1]
    elif win._verbose:
        print('Test done \n')


if __name__ == "__main__":
    main()
