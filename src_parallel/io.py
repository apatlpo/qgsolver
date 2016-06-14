#!/usr/bin/python
# -*- encoding: utf8 -*-

from petsc4py import PETSc

import numpy as np
from netCDF4 import Dataset


def write_nc(V, vname, filename, qg):    
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
    rank = qg.comm.getRank()

    if rank == 0:

        ### create a netcdf file to store QG pv for inversion
        rootgrp = Dataset(filename, 'w',
                          format='NETCDF4_CLASSIC', clobber=True)

        # create dimensions
        rootgrp.createDimension('x', qg.grid.Nx)
        rootgrp.createDimension('y', qg.grid.Ny)
        rootgrp.createDimension('z', qg.grid.Nz)
        
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
            nc_V.append(rootgrp.createVariable(name,dtype,('x','y','z',)))

    # loop around variables now and store them
    for i in xrange(Nv):
        Vn = qg.da.createNaturalVec()
        qg.da.globalToNatural(V[i], Vn)
        scatter, Vn0 = PETSc.Scatter.toZero(Vn)
        scatter.scatter(Vn, Vn0, False, PETSc.Scatter.Mode.FORWARD)
        if rank == 0:
            #nx, ny, nz = qg.da.sizes
            #print qg.da.sizes
            Vf = Vn0[...].reshape(qg.da.sizes, order='f')
            nc_V[i][:]=Vf[:,:,:]
        else:
            Vf=np.zeros(qg.da.sizes)
        qg.comm.barrier()
      
    if rank == 0:
        # close the netcdf file
        rootgrp.close()
        

