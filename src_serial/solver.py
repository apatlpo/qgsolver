#!/usr/bin/python
# -*- encoding: utf8 -*-


import sys
import petsc4py
petsc4py.init(sys.argv)
from petsc4py import PETSc

#from .grid import *


#
#==================== Serial solver ============================================
#

class pvinversion_serial():
    """ PV inversion, serial
    """
    
    def __init__(self, grid):
        
        # get petsc options from command line
        OptDB = PETSc.Options()

        # determine the tile decomposition        
        n  = OptDB.getInt('n', 16)
        #nx = OptDB.getInt('nx', n)
        #ny = OptDB.getInt('ny', n)
        #nz = OptDB.getInt('nz', n)
        
        #kplt = OptDB.getInt('kplt', nz//2)
        
        
        ### setup the solver
        
        #da = PETSc.DMDA().create([nx, ny, nz], stencil_width=2)
        da = PETSc.DMDA().create([grid.Nx, grid.Ny, grid.Nz], stencil_width=2)
        comm = da.getComm()
        rank = comm.getRank()
        
        # create the operator
        L = da.createMat()
        print 'Operator L created'
#         set_L()
#         
#         # create solver
#         ksp = PETSc.KSP()
#         ksp.create(PETSc.COMM_WORLD)
#         ksp.setOperators(L)
#         # use conjugate gradients
#         ksp.setType('cg')
#         #ksp.setType('gmres')
#         # and incomplete Cholesky for preconditionning
#         #ksp.getPC().setType('icc')
#         # set tolerances
#         #ksp.setTolerances(rtol=1e-10) # nope
#         ksp.setFromOptions()
#         
#         
#         ### setup time stepping
#         
#         # 4 steps explicit RungeKutta
#         a = [1./6., 1./3., 1./3., 1./6.]
#         b = [0.5, 0.5, 1.]
#         
#         # vector containing PV
#         Q = da.createGlobalVec()
#         Q0 = da.createGlobalVec()
#         Q1 = da.createGlobalVec()
#         # RHS
#         dQ = da.createGlobalVec()
#         # PV inv RHS vector
#         Qinv = da.createGlobalVec()
#         # vector containing the streamfunction
#         PSI = da.createGlobalVec()
#         U = da.createGlobalVec()
#         V = da.createGlobalVec()
#         # local vectors
#         localQ  = da.createLocalVec()
#         localdQ  = da.createLocalVec()
#         localPSI  = da.createLocalVec()
#         localU  = da.createLocalVec()
#         localV  = da.createLocalVec()
#         # for plotting purposes
#         PSIn = da.createNaturalVec()
#         Qn = da.createNaturalVec()
#         
#         # init PV
#         init_q()
#         set_q_bdy()
#         # compute corresponding PSI
#         Q.copy(Qinv) # copies Q into Qinv
#         set_qinv_bdy()
#         ksp.solve(Qinv, PSI)



class time_stepper_serial():
    """ Time stepper, serial
    """
    
    def __init__(self):
        pass
    




#
#==================== Parallel solver ==========================================
#

class pvinversion_parallel():
    """ PV inversion, parallel with petsc4py
    """
    
    def __init__(self):
        pass



class time_stepper_parallel():
    """ Time stepper, parallel with petsc4py
    """
    
    def __init__(self):
        pass
    
    