#!/usr/bin/python
# -*- encoding: utf8 -*-

"""
Create an operator L=I and multiply 
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
        
        # fill psi
        self.fill_PSI()
        
        # create and fill the operator
        self.L = self.da.createMat()
        self.set_L()
            
    
    def fill_PSI(self):
        psi = self.da.getVecArray(self.PSI)
        mx, my, mz = self.da.getSizes()
        (xs, xe), (ys, ye), (zs, ze) = self.da.getRanges()
        #
        #def f(i,j,k): return np.sin(np.pi*i/mx)*np.sin(np.pi*j/my)*np.sin(np.pi*k/mz)
        #
        for k in range(zs,ze):
            for j in range(ys, ye):
                for i in range(xs, xe):
                    #psi[i, j, k] = f(i,j,k)
                    psi[i, j, k] = 1.


    def set_L(self):
        #
        self.L.zeroEntries()
        #
        (xs, xe), (ys, ye), (zs, ze) = self.da.getRanges()
        row = PETSc.Mat.Stencil()
        for k in range(zs, ze):
            for j in range(ys, ye):
                for i in range(xs, xe):
                    row.index = (i,j,k)
                    row.field = 0
                    self.L.setValueStencil(row, row, 1.0)
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
    pb.write(pb.PSI, 'psi', 'PSI.nc')
    #
    pb.L.mult(pb.PSI, pb.RHS)
    pb.write(pb.RHS, 'Lpsi', 'Lpsi.nc')
    if pb.rank==0: print 'done'
        
    

    
