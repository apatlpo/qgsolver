#!/usr/bin/python
# -*- encoding: utf8 -*-

from petsc4py import PETSc

def set_L(L, qg):
    """ Builds the laplacian operator along with boundary conditions """
    #
    mx, my, mz = qg.da.getSizes()
    dx, dy, dz = qg.grid.dx, qg.grid.dy, qg.grid.dz
    idx, idy, idz = [1.0/dl for dl in [dx, dy, dz]]
    idx2, idy2, idz2 = [1.0/dl**2 for dl in [dx, dy, dz]]
    #print idx, idy, idz 
    #print idx2, idy2, idz2   
    #
    (xs, xe), (ys, ye), (zs, ze) = qg.da.getRanges()
    #
    L.zeroEntries()
    row = PETSc.Mat.Stencil()
    col = PETSc.Mat.Stencil()
    #
    for k in range(zs, ze):
        for j in range(ys, ye):
            for i in range(xs, xe):
                row.index = (i,j,k)
                row.field = 0
                if (k==0):
                    # bottom bdy condition: Neuman dpsi/dz=0
                    for index, value in [
                        ((i,j,k), -idz),
                        ((i,j,k+1),  idz)
                        ]:
                        col.index = index
                        col.field = 0
                        L.setValueStencil(row, col, value)
                elif (k==mz-1):
                    # top bdy condition: Neuman dpsi/dz=0
                    for index, value in [
                        ((i,j,k-1), -idz),
                        ((i,j,k),  idz),
                        ]:
                        col.index = index
                        col.field = 0
                        L.setValueStencil(row, col, value)
                elif (i==0    or j==0 or
                      i==mx-1 or j==my-1):
                    # Dirichlet lateral bdy condition: psi=0
                    L.setValueStencil(row, row, 1.0)
                else:
                    # interior pv
                    for index, value in [
                        ((i,j,k-1), qg._sparam[k]*idz2),
                        ((i,j-1,k), idy2),
                        ((i-1,j,k), idx2),
                        ((i, j, k), -2.*(idx2+idy2)-(qg._sparam[k]*idz2+qg._sparam[k+1]*idz2)),
                        ((i+1,j,k), idx2),
                        ((i,j+1,k), idy2),
                        ((i,j,k+1), qg._sparam[k+1]*idz2)
                        ]:
                        col.index = index
                        col.field = 0
                        L.setValueStencil(row, col, value)
    L.assemble()
    return