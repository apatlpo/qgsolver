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
    
    #qg = qg()
    #qg = qg(grid_uniform_input = {'Lx':1.e8})
    qg = qg(grid_uniform_input = {'Nx':100, 'Ny':150, 'Nz':10 })
    qg.set_q()
    qg.pvinv.solve(qg)
    write_nc(qg.PSI,'psi','psi.nc', qg)

    if qg._verbose:
        #print '%e \n' % qg.grid.H
        #print dir(qg)
        pass

