#!/usr/bin/python
# -*- encoding: utf8 -*-

""" Test the basic features of the library:
Setup of uniform grid
PV inversion of an analytical PV distribution
"""

#import qgsolver.qg as qg
from qgsolver.qg import qg
from qgsolver.io import write_nc

if __name__ == "__main__":
    
    qg = qg(grid_uniform_input = {'Nx':100, 'Ny':150, 'Nz':3 },
            K = 1.e1, dt = 86400.e1)
    #
    qg.set_q()
    qg.invert_pv()
    write_nc([qg.PSI, qg.Q], ['q', 'psi'], 'output.nc', qg)
    #
    qg.tstep(1)
    write_nc([qg.PSI, qg.Q], ['q', 'psi'], 'output.nc', qg, create=False)

    if qg._verbose:
        print 'Test done \n'
