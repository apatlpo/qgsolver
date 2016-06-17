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
    
    qg = qg(grid_uniform_input = {'Nx':150, 'Ny':100, 'Nz':3 },
            K = 0.e0, dt = 5.*86400.e0)
    #
    qg.set_q()
    qg.invert_pv()
    write_nc([qg.PSI, qg.Q], ['psi', 'q'], 'output.nc', qg)
    #
    if False:
        qg.tstep(1)
        write_nc([qg.PSI, qg.Q], ['psi', 'q'], 'output.nc', qg, create=False)
    else:
        while qg.tstepper.t/86400. < 100 :
            qg.tstep(1)
            write_nc([qg.PSI, qg.Q], ['psi', 'q'], 'output.nc', qg, create=False)

    if qg._verbose:
        print 'Test done \n'
