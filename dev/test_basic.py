#!/usr/bin/python
# -*- encoding: utf8 -*-

""" Test the basic features of the library:
Setup of uniform grid
PV inversion of an analytical PV distribution
"""

#import qgsolver.qg as qg
from qgsolver.qg import qg

if __name__ == "__main__":
    #qg = qg()
    #qg = qg(grid_uniform_input = {'Lx':1.e8})
    qg = qg(grid_uniform_input = {'Nx':10, 'Ny':15, 'Nz':2 })
    print qg.grid.H
    qg.set_psi()
    print dir(qg)

