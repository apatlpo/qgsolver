#!/usr/bin/python
# -*- encoding: utf8 -*-
"""
Create a window for spectral computations in the presence of a mask

The window \psi satisfies:

\psi = \psi_0 / max(\psi_0)

where:

\Delta \psi_0 - k^2 \psi_0 = 1
\psi_0 = 0 on frontiers and land

Run with:
mpirun -n 2 python window_uniform.py
mpirun -n 2 python window_uniform.py -mf -ksp_view -ksp_monitor -ksp_converged_reason

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
    """ Tests with curvilinear grid
    """


    # set tiling
    ncores_x, ncores_y= 2, 1 #   linux
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
    casename = "window_uniform"

    #datapath = "data/"

    # from notebook:
    hgrid = {"Lx": (200-1)*11100, "Nx":200,
             "Ly": (200-1)*9618,"Ny":200}

    vgrid = {"H": 1000, "Nz":10}
    # horizontal subdomain
    hdom = {"istart": 0, "iend": 100-1, "i0": 0,
            "jstart": 0, "jend": 100-1,  "j0":0,
            }
    # vertical subdomain
    vdom = {"kdown": 0, "kup": 10-1, "k0":0}
    win = window(hgrid=hgrid, vgrid=vgrid,
                 K=1.e-6,
                 mask3D=False,
                 mask_file="mask.nc",
                 vdom=vdom, hdom=hdom,
                 ncores_x=ncores_x, ncores_y=ncores_y,
                 )
    win.case=casename

    if win.rank == 0: print("----------------------------------------------------")
    if win.rank == 0: print("Elapsed time for window_model ", str(time.time() - cur_time))
    cur_time = time.time()

    win.set_q()
    if win.rank == 0: print("----------------------------------------------------")
    if win.rank == 0: print("Elapsed time for set_q (//RHS)", str(time.time() - cur_time))
    cur_time = time.time()

    win.invert_win()
    if win.rank == 0: print("----------------------------------------------------")
    if win.rank == 0: print("Elapsed time for invert_win ", str(time.time() - cur_time))
    cur_time = time.time()

    # normalize by its maximum value
    Wmax = win.PSI.norm(norm_type=PETSc.NormType.INF)
    if win.rank == 0: print("----------------------------------------------------")
    if win.rank == 0: print("Infinite norm of window is =%f" %Wmax)
    win.PSI.scale(1./Wmax)
    if win.rank == 0: print("Window is normalized by this infinite norm")

    # output the window
    write_nc([win.PSI], ["window"], "output.nc", win.da, win.grid)
    if win.rank == 0: print("----------------------------------------------------")
    if win.rank == 0: print("Elapsed time for write_nc ", str(time.time() - cur_time))
    cur_time = time.time()

    if win.rank == 0: print("----------------------------------------------------")
    if win.rank == 0: print("Elapsed time  ", str(cur_time - start_time))

    return win

def main(ping_mpi_cfg=False):

    win = window_runs(ping_mpi_cfg=ping_mpi_cfg)

    if ping_mpi_cfg:
        return win[0], win[1]
    elif win._verbose:
        print("Test done \n")


if __name__ == "__main__":
    main()
