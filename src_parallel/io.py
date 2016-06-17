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
        

