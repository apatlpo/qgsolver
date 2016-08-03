#!/usr/bin/python
# -*- encoding: utf8 -*-

from petsc4py import PETSc

import numpy as np
from netCDF4 import Dataset


def write_nc(V, vname, filename, qg, create=True):    
    """ Write a variable to a netcdf file
    Parameters:
        V list of petsc vectors
        vname list of corresponding names
        filename
        qg object
    """

    # number of variables to be stored
    Nv=len(vname)
    # process rank
    rank = qg.rank

    if rank == 0 and create:

        ### create a netcdf file to store QG pv for inversion
        rootgrp = Dataset(filename, 'w',
                          format='NETCDF4_CLASSIC', clobber=True)

        # create dimensions
        rootgrp.createDimension('x', qg.grid.Nx)
        rootgrp.createDimension('y', qg.grid.Ny)
        rootgrp.createDimension('z', qg.grid.Nz)
        rootgrp.createDimension('t', None)
        
        # create variables
        dtype='f8'
        nc_x = rootgrp.createVariable('x',dtype,('x'))
        nc_y = rootgrp.createVariable('y',dtype,('y'))
        nc_z = rootgrp.createVariable('z',dtype,('z'))
        #x,y,z=qg.grid.get_xyz()
        nc_x[:], nc_y[:], nc_z[:] = qg.grid.get_xyz()
        # 3D variables
        nc_V=[]
        for name in vname:
            nc_V.append(rootgrp.createVariable(name,dtype,('t','z','y','x',)))
    
    elif rank == 0:
        ### open netcdf file
        rootgrp = Dataset(filename, 'a',
                          format='NETCDF4_CLASSIC')
        # 3D variables
        nc_V=[]
        for name in vname:
            nc_V.append(rootgrp.variables[name])
        

    # loop around variables now and store them
    Vn = qg.da.createNaturalVec()
    for i in xrange(Nv):    
        qg.da.globalToNatural(V[i], Vn)
        scatter, Vn0 = PETSc.Scatter.toZero(Vn)
        scatter.scatter(Vn, Vn0, False, PETSc.Scatter.Mode.FORWARD)
        if rank == 0:
            Vf = Vn0[...].reshape(qg.da.sizes[::-1], order='c')
            if create:
                nc_V[i][:] = Vf[np.newaxis,...]
            else:
                if i==0: it=nc_V[i].shape[0]
                nc_V[i][it,...] = Vf[:]
        qg.comm.barrier()
      
    if rank == 0:
        # close the netcdf file
        rootgrp.close()
        


def read_nc_petsc(V, vname, filename, qg):    
    """ Read a variable from a netcdf file and stores it in a petsc Vector
    Parameters:
        V one(!) petsc vector
        vname corresponding name in netcdf file
        filename
        qg object
    """
    v = qg.da.getVecArray(V)
    mx, my, mz = qg.da.getSizes()
    (xs, xe), (ys, ye), (zs, ze) = qg.da.getRanges()
    rootgrp = Dataset(filename, 'r')
    for k in range(zs, ze):
        for j in range(ys, ye):
            for i in range(xs, xe):
                #v[i, j, k] = rootgrp.variables['q'][-1,k,j,i]
                # line above does not work for early versions of netcdf4 python library
                # print netCDF4.__version__  1.1.1 has a bug and one cannot call -1 for last index:
                # https://github.com/Unidata/netcdf4-python/issues/306
                v[i, j, k] = rootgrp.variables['q'][rootgrp.variables['q'].shape[0]-1,k,j,i]
    rootgrp.close()
    qg.comm.barrier()
    #if qg.rank ==0: print 'Variable '+vname+' read from '+filename
    if qg._verbose: print '... done'
#     # number of variables to read
#     Nv=len(vname)
#     # process rank
#     rank = qg.rank
# 
#     if rank == 0:
# 
#         ### create a netcdf file to store QG pv for inversion
#         rootgrp = Dataset(filename, 'r')
# 
#         # 3D variables
#         nc_V=[]
#         for name in vname:
#             nc_V.append(rootgrp.createVariable(name,dtype,('t','z','y','x',)))
#     
#     # loop around variables now and store them
#     Vn = qg.da.createNaturalVec()
#     for i in xrange(Nv):
#         qg.da.globalToNatural(V[i], Vn)
#         #qg.da.naturalToGlobal(Vn,V[i])
#         scatter, Vn0 = PETSc.Scatter.toAll(Vn)
#         if rank == 0:
#             Vf = nc_V[i][-1,...]
#         scatter.scatter(Vn0, Vn, False, PETSc.Scatter.Mode.FORWARD)
#         qg.comm.barrier()
#       
#     if rank == 0:
#         # close the netcdf file
#         rootgrp.close()
        
        
        
def read_nc(vnames, filename):
    """ Read variables from a netcdf file
    Parameters:
        vnames list of variable names
        filename
    """

    # open netdc file
    rootgrp = Dataset(filename, 'r')
    
    # loop around variables to load
    V=[]
    for name in vnames:
        V.append(rootgrp.variables[name])
    
    # close the netcdf file
    rootgrp.close()
    
    return V
        