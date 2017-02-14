#!/usr/bin/python
# -*- encoding: utf8 -*-

from petsc4py import PETSc
import sys


def set_L(L, qg):
    """ Builds the laplacian operator along with boundary conditions
        Horizontally uniform grid
    """
    
    if qg._verbose>0:
        print '  ... assumes a uniform horizontal and vertical grid'
    
    #
    mx, my, mz = qg.da.getSizes()
    dx, dy, dz = qg.grid.dx, qg.grid.dy, qg.grid.dz
    idx, idy, idz = [1.0/dl for dl in [dx, dy, dz]]
    idx2, idy2, idz2 = [1.0/dl**2 for dl in [dx, dy, dz]]
    #
    (xs, xe), (ys, ye), (zs, ze) = qg.da.getRanges()
    #
    istart = qg.grid.istart
    iend = qg.grid.iend
    jstart = qg.grid.jstart
    jend = qg.grid.jend
    kdown = qg.grid.kdown
    kup = qg.grid.kup
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

                # lateral points outside the domain: dirichlet, psi=...
                if (i<=istart or j<=jstart or
                    i>=iend or j>=jend):
                    L.setValueStencil(row, row, 1.0)

                # bottom bdy condition: default Neuman dpsi/dz=...
                elif (k==kdown):
                    if qg.bdy_type['bottom']=='N' :
                        for index, value in [
                            ((i,j,k), -idz),
                            ((i,j,k+1),  idz)
                            ]:
                            col.index = index
                            col.field = 0
                            L.setValueStencil(row, col, value)
                    elif qg.bdy_type['bottom']=='D':
                        L.setValueStencil(row, row, 1.0)
                    else:
                        print "unknown bottom boundary condition"
                        sys.exit()

                # top bdy condition: default Neuman dpsi/dz=...
                elif (k==kup):
                    if qg.bdy_type['top']=='N' :
                        for index, value in [
                            ((i,j,k-1), -idz),
                            ((i,j,k),  idz),
                            ]:
                            col.index = index
                            col.field = 0
                            L.setValueStencil(row, col, value)
                    elif qg.bdy_type['top']=='D':
                        L.setValueStencil(row, row, 1.0)
                    else:
                        print "unknown bottom boundary condition"
                        sys.exit()

                # points below and above the domain
                elif (k<kdown or k>kup):
                    L.setValueStencil(row, row, 0.0)

                # interior points: pv is prescribed
                else:
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





def set_L_curv(L, qg):
    """ Builds the laplacian operator along with boundary conditions
        Horizontally uniform grid
    """
    
    if qg._verbose>0:
        print '  ... assumes a curvilinear and/or vertically stretched grid'
    #
    mx, my, mz = qg.da.getSizes()
    #
    #D = qg.da.getVecArray(qg.grid.D)
    local_D  = qg.da.createLocalVec()
    qg.da.globalToLocal(qg.grid.D, local_D)
    D = qg.da.getVecArray(local_D)
    kdx = qg.grid._k_dx
    kdy = qg.grid._k_dy
    #
    idzc = 1./qg.grid.dzc
    idzf = 1./qg.grid.dzf
    #idz, idz2 = 1./dz, 1./dz**2
    #
    (xs, xe), (ys, ye), (zs, ze) = qg.da.getRanges()
    istart = qg.grid.istart
    iend = qg.grid.iend
    jstart = qg.grid.jstart
    jend = qg.grid.jend
    kdown = qg.grid.kdown
    kup = qg.grid.kup
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

                # lateral points outside the domain: dirichlet, psi=...
                if (i<=istart or j<=jstart or
                    i>=iend or j>=jend):
                    L.setValueStencil(row, row, 1.0)

                # bottom bdy condition: default Neuman dpsi/dz=...
                elif (k==kdown):
                    if qg.bdy_type['bottom']=='N' : 
                        for index, value in [
                            ((i,j,k), -idzc[k]),
                            ((i,j,k+1),  idzc[k])
                            ]:
                            col.index = index
                            col.field = 0
                            L.setValueStencil(row, col, value)
                    elif qg.bdy_type['bottom']=='D' :
                        L.setValueStencil(row, row, 1.0)
                    else:
                        print "unknown bottom boundary condition"
                        sys.exit()

                # top bdy condition: default Neuman dpsi/dz=...
                elif (k==kup):
                    if qg.bdy_type['top']=='N' : 
                        for index, value in [
                            ((i,j,k-1), -idzc[k-1]),
                            ((i,j,k),  idzc[k-1]),
                            ]:
                            col.index = index
                            col.field = 0
                            L.setValueStencil(row, col, value)
                    elif qg.bdy_type['top']=='D':
                        L.setValueStencil(row, row, 1.0)
                    else:
                        print "unknown top boundary condition"
                        sys.exit()

                # points below and above the domain
                elif (k<kdown) or (k>kup):
                    L.setValueStencil(row, row, 0.0)
                    # L.setValueStencil(row, row, 1.0)
                    
                # interior points: pv is prescribed
                else:
                    
                    for index, value in [
                        ((i,j,k-1), qg._sparam[k]*idzc[k-1]*idzf[k]),
                        ((i,j-1,k), 1./D[i,j,kdx]/D[i,j,kdy] \
                         * (D[i,j-1,kdx]+D[i,j,kdx])/(D[i,j-1,kdy]+D[i,j,kdy])),
                        ((i-1,j,k), 1./D[i,j,kdx]/D[i,j,kdy] \
                         * (D[i-1,j,kdy]+D[i,j,kdy])/(D[i-1,j,kdx]+D[i,j,kdx])),
                        ((i, j, k), -1./D[i,j,kdx]/D[i,j,kdy]*( \
                                         (D[i,j-1,kdx]+D[i,j,kdx])/(D[i,j-1,kdy]+D[i,j,kdy]) \
                                        +(D[i-1,j,kdy]+D[i,j,kdy])/(D[i-1,j,kdx]+D[i,j,kdx]) \
                                        +(D[i,j,kdy]+D[i+1,j,kdy])/(D[i+1,j,kdx]+D[i,j,kdx]) \
                                        +(D[i,j,kdx]+D[i,j+1,kdx])/(D[i,j,kdy]+D[i,j+1,kdy]) \
                                                               ) \
                         - (qg._sparam[k]*idzc[k-1]*idzf[k]+qg._sparam[k+1]*idzc[k]*idzf[k])),
                        ((i+1,j,k), 1./D[i,j,kdx]/D[i,j,kdy] \
                         * (D[i,j,kdy]+D[i+1,j,kdy])/(D[i+1,j,kdx]+D[i,j,kdx])),
                        ((i,j+1,k), 1./D[i,j,kdx]/D[i,j,kdy] \
                         * (D[i,j,kdx]+D[i,j+1,kdx])/(D[i,j,kdy]+D[i,j+1,kdy])),
                        ((i,j,k+1), qg._sparam[k+1]*idzc[k]*idzf[k])
                        ]:
                        col.index = index
                        col.field = 0
                        L.setValueStencil(row, col, value)

               
    L.assemble()
    return


