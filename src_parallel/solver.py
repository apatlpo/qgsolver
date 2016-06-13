#!/usr/bin/python
# -*- encoding: utf8 -*-


import sys
#import petsc4py
#from petsc4py import PETSc

#from .grid import *
from .set_L import *

#
#==================== Serial solver ============================================
#

class pvinversion():
    """ PV inversion, parallel
    """
    
    def __init__(self, qg):
        
        # get petsc options from command line
        #OptDB = PETSc.Options()

        # determine the tile decomposition        
        #n  = OptDB.getInt('n', 16)
        #nx = OptDB.getInt('nx', n)
        #ny = OptDB.getInt('ny', n)
        #nz = OptDB.getInt('nz', n)
        
        #kplt = OptDB.getInt('kplt', nz//2)
        
        
        ### setup the solver
        
        #da = PETSc.DMDA().create([nx, ny, nz], stencil_width=2)
        #da = PETSc.DMDA().create([grid.Nx, grid.Ny, grid.Nz], stencil_width=2)
        #comm = da.getComm()
        #rank = comm.getRank()
        
        # create the operator
        #print dir(qg)
        self.L = qg.da.createMat()
        #
        if qg._verbose>0:
            print 'Operator L declared \n'

        # Fill in operator values
        self.L = set_L(self.L, qg.da)
        #
        if qg._verbose>0:
            print 'Operator L filled \n'

        # local vectors
        self._localQ  = qg.da.createLocalVec()
        self._localPSI  = qg.da.createLocalVec()


        # create solver
        self.ksp = PETSc.KSP()
        self.ksp.create(PETSc.COMM_WORLD)
        self.ksp.setOperators(self.L)
        # use conjugate gradients
        self.ksp.setType('cg')
        #self.ksp.setType('gmres')
        # and incomplete Cholesky for preconditionning
        #self.ksp.getPC().setType('icc')
        # set tolerances
        #self.ksp.setTolerances(rtol=1e-10) # nope
        self.ksp.setFromOptions()
         
         
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

    def solve(self, qg):
        """ Compute the PV inversion
        """
        self.ksp.solve(qg.Q, qg.PSI)
        if qg._verbose>0:
            print 'Inversion done'




#
#==================== Time stepper ============================================
#


class time_stepper():
    """ Time stepper, parallel with petsc4py
    """
    
    def __init__(self):
        pass
    
    