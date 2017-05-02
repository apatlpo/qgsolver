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
                        print "unknown top boundary condition"
                        sys.exit()

                # points below and above the domain
                elif (k<kdown or k>kup):
                    L.setValueStencil(row, row, 0.0)

                # interior points: pv is prescribed
                else:
                    for index, value in [
                        ((i,j,k-1), qg._sparam[k-1]*idz2),
                        ((i,j-1,k), idy2),
                        ((i-1,j,k), idx2),
                        ((i, j, k), -2.*(idx2+idy2)-(qg._sparam[k]*idz2+qg._sparam[k-1]*idz2)),
                        ((i+1,j,k), idx2),
                        ((i,j+1,k), idy2),
                        ((i,j,k+1), qg._sparam[k]*idz2)
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
    kmask = qg.grid._k_mask
    kdxu = qg.grid._k_dxu
    kdyu = qg.grid._k_dyu
    kdxv = qg.grid._k_dxv
    kdyv = qg.grid._k_dyv
    kdxt = qg.grid._k_dxt
    kdyt = qg.grid._k_dyt
    #
    idzt = 1./qg.grid.dzt
    idzw = 1./qg.grid.dzw
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

                # masked points (land=0), L=1
                if D[i,j,kmask]==0.:
                    L.setValueStencil(row, row, 1.)

                # lateral points outside the domain: dirichlet, psi=...
                elif (i<=istart or j<=jstart or
                    i>=iend or j>=jend):
                    L.setValueStencil(row, row, 1.0)

                # bottom bdy condition: default Neuman dpsi/dz=...
                elif (k==kdown):
                    if qg.bdy_type['bottom']=='N' : 
                        for index, value in [
                            ((i,j,k), -idzw[k]),
                            ((i,j,k+1),  idzw[k])
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
                            ((i,j,k-1), -idzw[k-1]),
                            ((i,j,k),  idzw[k-1]),
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
                    # L.setValueStencil(row, row, 0.0)
                    L.setValueStencil(row, row, 1.0)

                # lateral points outside the domain: dirichlet, psi=...
                # elif (i<=istart or j<=jstart or
                #     i>=iend or j>=jend):
                #     L.setValueStencil(row, row, 1.0)

                # interior points: pv is prescribed
                else:
                    
                    for index, value in [

                        ((i,j,k-1), qg._sparam[k-1]*idzt[k]*idzw[k-1]),
                        ((i,j-1,k), 1./D[i,j,kdxt]/D[i,j,kdyt] * D[i,j-1,kdxv]/D[i,j-1,kdyv]),
                        ((i-1,j,k), 1./D[i,j,kdxt]/D[i,j,kdyt] * D[i-1,j,kdyu]/D[i-1,j,kdxu]),
                        ((i, j, k), -1./D[i,j,kdxt]/D[i,j,kdyt]*( \
                                         D[i,j,kdyu]/D[i,j,kdxu] \
                                        +D[i-1,j,kdyu]/D[i-1,j,kdxu] \
                                        +D[i,j,kdxv]/D[i,j,kdyv] \
                                        +D[i,j-1,kdxv]/D[i,j-1,kdyv])
                         - (qg._sparam[k]*idzt[k]*idzw[k]+qg._sparam[k-1]*idzt[k]*idzw[k-1])),
                        ((i+1,j,k), 1./D[i,j,kdxt]/D[i,j,kdyt] * D[i,j,kdyu]/D[i,j,kdxu]),
                        ((i,j+1,k), 1./D[i,j,kdxt]/D[i,j,kdyt] * D[i,j,kdxv]/D[i,j,kdyv]),
                        ((i,j,k+1), qg._sparam[k]*idzt[k]*idzw[k])
                        ]:
                        col.index = index
                        col.field = 0
                        L.setValueStencil(row, col, value)

               
    L.assemble()
    return
