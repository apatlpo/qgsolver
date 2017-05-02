#!/usr/bin/python
# -*- encoding: utf8 -*-

"""
Test solve a 3D poisson equation: 
Laplacian psi = sin(pi x) sin(pi y) sin (pi z) for 0<x<1, 0<y<1, 0<z<1
with Dirichlet boundary conditions: psi=0
The exact solution is psi = - sin(pi x) sin(pi y) sin (pi z) / 3 pi^2
"""

import time
import sys

import petsc4py
petsc4py.init(sys.argv)
from petsc4py import PETSc

import numpy as np
from netCDF4 import Dataset

class minipb(object):

    def __init__(self,N, ncores_x, ncores_y):
        
        self.N = N
        self.da = PETSc.DMDA().create(sizes = [N, N, N],
                                      proc_sizes = [ncores_x,ncores_y,1],
                                      stencil_width = 2)
        self.comm = self.da.getComm()
        self.rank = self.comm.getRank()
        if self.rank==0:
            print 'The 3D grid is tiled according to (nproc_x, nproc_y, nproc_z) : '\
                    +str(self.da.proc_sizes) 
        
        # init useful vectors
        self.PSI = self.da.createGlobalVec()
        self.RHS = self.da.createGlobalVec()
        
        # fill right-hand-side
        self.fill_RHS()
        
        # create and fill the operator
        self.L = self.da.createMat()
        self.set_L()
        
        # create the solver
        self.ksp = PETSc.KSP()
        self.ksp.create(PETSc.COMM_WORLD)
        self.ksp.setOperators(self.L)
        self.ksp.setType('bicg')
        self.ksp.setInitialGuessNonzero(True)
        self.ksp.setTolerances(max_it=1000)
        #
        for opt in sys.argv[1:]:
            PETSc.Options().setValue(opt, None)
        self.ksp.setFromOptions()       
        
        pass
    
    
    def fill_RHS(self):
        rhs = self.da.getVecArray(self.RHS)
        mx, my, mz = self.da.getSizes()
        (xs, xe), (ys, ye), (zs, ze) = self.da.getRanges()
        #
        def f(i,j,k): return np.sin(np.pi*i/mx)*np.sin(np.pi*j/my)*np.sin(np.pi*k/mz)
        #
        for k in range(zs,ze):
            for j in range(ys, ye):
                for i in range(xs, xe):
                    if (i>0 and j>0 and k>0 and i<mx-1 and j<my-1 and k<mz-1):
                        rhs[i, j, k] = f(i,j,k)
                    else:
                        rhs[i, j, k] = 0.
        

    def set_L(self):
        #
        self.L.zeroEntries()
        #
        mx, my, mz = self.da.getSizes()
        dx, dy, dz = 1./mx, 1./my, 1./mz
        #idx, idy, idz = [1.0/dl for dl in [dx, dy, dz]]
        idx2, idy2, idz2 = [1.0/dl**2 for dl in [dx, dy, dz]]
        #
        (xs, xe), (ys, ye), (zs, ze) = self.da.getRanges()
        row = PETSc.Mat.Stencil()
        col = PETSc.Mat.Stencil()
        for k in range(zs, ze):
            for j in range(ys, ye):
                for i in range(xs, xe):
                    row.index = (i,j,k)
                    row.field = 0
                    # bdy condition
                    if (i==0 or j==0 or k==0 or
                        i==mx-1 or j==my-1 or k==mz-1):
                        self.L.setValueStencil(row, row, 1.0)
                    # interior points
                    else:
                        for index, value in [
                        ((i,j,k-1), idz2),
                        ((i,j-1,k), idy2),
                        ((i-1,j,k), idx2),
                        ((i, j, k), -2.*(idx2+idy2+idz2)),
                        ((i+1,j,k), idx2),
                        ((i,j+1,k), idy2),
                        ((i,j,k+1), idz2)
                        ]:
                            col.index = index
                            col.field = 0
                            self.L.setValueStencil(row, col, value)
        self.L.assemble()
    
    
    def write(self, V, vname, filename):
        """ Write a variable to a netcdf file
        """
        # process rank
        rank = self.rank
        #
        if rank == 0:
            ### create a netcdf file to store QG pv for inversion
            rootgrp = Dataset(filename, 'w',
                              format='NETCDF4_CLASSIC', clobber=True)
            # create dimensions
            rootgrp.createDimension('x', self.N)
            rootgrp.createDimension('y', self.N)
            rootgrp.createDimension('z', self.N)
            #
            dtype='f8'
            nc_V = rootgrp.createVariable(vname,dtype,('z','y','x',))
        #
        Vn = self.da.createNaturalVec()
        self.da.globalToNatural(V, Vn)
        scatter, Vn0 = PETSc.Scatter.toZero(Vn)
        scatter.scatter(Vn, Vn0, False, PETSc.Scatter.Mode.FORWARD)
        #    
        if rank == 0:
            Vf = Vn0[...].reshape(self.da.sizes[::-1], order='c')
            nc_V[:] = Vf[...]
        self.comm.barrier()
        #    
        if rank == 0:
            rootgrp.close()



if __name__ == "__main__":
    
    pb = minipb(64, 4,2)
    pb.write(pb.RHS, 'rhs', 'RHS.nc')
    #
    pb.ksp.solve(pb.RHS, pb.PSI)
    pb.write(pb.PSI, 'psi', 'PSI.nc')
    #
    pb.L.mult(pb.PSI, pb.RHS)
    pb.write(pb.RHS, 'Lpsi', 'Lpsi.nc')
    if pb.rank==0: print 'done'
        
    

    
