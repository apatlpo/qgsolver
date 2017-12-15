#!/usr/bin/python
# -*- encoding: utf8 -*-

import sys, os
from petsc4py import PETSc

import numpy as np
from netCDF4 import Dataset
import netCDF4
import time

#
#==================== Pure IO ============================================
#

def write_nc(V, vname, filename, qg, create=True):    
    """ Write a variable to a netcdf file

    Parameters:
        V list of petsc vectors (may contain None)
        vname list of corresponding names
        filename
        qg object
    """

    # number of variables to be stored
    Nv=len(vname)
    # process rank
    rank = qg.rank
    if qg.grid.mask:
        # get global mask for rank 0 (None for other proc)
        global_D = get_global(qg.grid.D, qg)

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
        if qg.grid.mask:
            # 2D mask
            nc_mask = rootgrp.createVariable('mask',dtype,('y','x'))
            nc_mask[:]= global_D[qg.grid._k_mask,:,:]
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
    for i in range(Nv):
        #  get global variable for rank 0 (None for other proc) 
        vglobal = get_global(V[i], qg)
        if rank == 0:
            if create:
                nc_V[i][:] = vglobal[np.newaxis,...]
            else:
                if i==0: it=nc_V[i].shape[0]
                nc_V[i][it,...] = vglobal[:]
         
    if rank == 0:
        # close the netcdf file
        rootgrp.close()

def read_nc_petsc(V, vname, filename, qg, fillmask=None):    
    """
    Read a variable from a netcdf file and stores it in a petsc Vector

    Parameters:
        V one(!) petsc vector
        vname corresponding name in netcdf file
        filename
        qg object
        fillmask : value to replace default netCDF fill value transformed in Nan
    """
    (xs, xe), (ys, ye), (zs, ze) = qg.da.getRanges()
    istart = xs + qg.grid.i0
    iend = xe + qg.grid.i0
    jstart = ys + qg.grid.j0
    jend = ye + qg.grid.j0
    kdown = zs + qg.grid.k0
    kup = ze + qg.grid.k0

    if os.path.isfile(filename): 
        v = qg.da.getVecArray(V)
        rootgrp = Dataset(filename, 'r')
        ndim=len(rootgrp.variables[vname].shape)
        if ndim>3:
        #v[i, j, k] = rootgrp.variables['q'][-1,k,j,i]
        # line above does not work for early versions of netcdf4 python library
        # print netCDF4.__version__  1.1.1 has a bug and one cannot call -1 for last index:
        # https://github.com/Unidata/netcdf4-python/issues/306
            vread = rootgrp.variables[vname][rootgrp.variables[vname].shape[0]-1,kdown:kup,jstart:jend,istart:iend]
        else:
            vread = rootgrp.variables[vname][kdown:kup,jstart:jend,istart:iend]        
        #vread = rootgrp.variables[vname][kdown:kup,jstart:jend,istart:iend]
        # replace the default fill value of netCDF by the input fillmask value
        if fillmask is not None:
            mx = np.ma.masked_values (vread, netCDF4.default_fillvals['f8'])
            vread[:]=mx.filled(fill_value=fillmask)
        for k in range(zs, ze):
            for j in range(ys, ye):
                for i in range(xs, xe):
                    v[i, j, k] = vread[k-zs,j-ys,i-xs]                  
    
        rootgrp.close()
    else:
        print('Error: read '+vname+': '+filename+' does not exist. Program will stop.')
        sys.exit()    
    qg.comm.barrier()

def read_nc_petsc_2D(V, vname, filename, level, qg):    
    """
    Read a 2D variable from a netcdf file and stores it in a petsc 3D Vector at k=level

    Parameters:
        V one(!) petsc vector
        vname corresponding name in netcdf file
        filename
        level vertical index to store 2D variable in 3D petsc vector
        qg object
    """
    v = qg.da.getVecArray(V)
    (xs, xe), (ys, ye), (zs, ze) = qg.da.getRanges()
    istart = xs + qg.grid.i0
    iend = xe + qg.grid.i0
    jstart = ys + qg.grid.j0
    jend = ye + qg.grid.j0

    if os.path.isfile(filename): 
        rootgrp = Dataset(filename, 'r')
        ndim=len(rootgrp.variables[vname].shape)
        if ndim==2:
            vread = rootgrp.variables[vname][jstart:jend,istart:iend]
        else:
            print('Error in read_nc_petsc_2D')
            sys.exit()
        for j in range(ys, ye):
            for i in range(xs, xe):
                v[i, j, level] = vread[j-ys,i-xs]                  
        rootgrp.close()
    else:
        print('Error: read '+vname+': '+filename+' does not exist. Program will stop.')
        sys.exit()         
               
def read_nc(vnames, filename,qg):
    """ Read variables from a netcdf file

    Parameters:
        vnames list of variable names
        filename
    """
   
    if os.path.isfile(filename):
        # open netdc file
        rootgrp = Dataset(filename, 'r')
        
        # loop around variables to load
        kstart = qg.grid.k0
        kend = qg.grid.k0 + qg.grid.Nz
    
        if isinstance(vnames, list):
            V=[]
            for name in vnames:
                if name == 'N2':
                    V.append(rootgrp.variables[name][kstart:kend])
                elif name == 'f0':
                    V.append(rootgrp.variables[name][:])
                # elif name == 'zt' or name == 'zw':
                #     V.append(rootgrp.variables[name][kstart:kend])
                else:
                    print('Error in read_nc: unknown variable '+name)
                    sys.exit()
        else:
            if vnames == 'N2':
                V = rootgrp.variables[vnames][kstart:kend]
            elif vnames == 'f0':
                V = rootgrp.variables[vnames][:]
            # elif vnames == 'zt' or vnames == 'zw':
            #     V = rootgrp.variables[vnames][kstart:kend]
            else:
                print('error in read_nc: unknown variable '+vnames)
                sys.exit()
    
        # close the netcdf file
        rootgrp.close()
        
        return V
    else:
        print('Error: read '+vnames+': '+filename+' does not exist. Program will stop.')
        sys.exit()


def read_hgrid_dimensions(hgrid_file):
    """ Reads grid dimension from netcdf file
    Could put dimension names as optional inputs ...
    """
    # open netcdf file
    rootgrp = Dataset(hgrid_file, 'r')
    Nx = len(rootgrp.dimensions['x'])
    Ny = len(rootgrp.dimensions['y'])    
    return Nx, Ny



def get_global(V, qg):
    """
    Returns a copy of the global V array on process 0, otherwise returns None
    """    
    rank = qg.rank
    Vn = qg.da.createNaturalVec()
    qg.da.globalToNatural(V, Vn)
    scatter, Vn0 = PETSc.Scatter.toZero(Vn)
    scatter.scatter(Vn, Vn0, False, PETSc.Scatter.Mode.FORWARD)
    if rank == 0:
        Vf = Vn0[...].reshape(qg.da.sizes[::-1], order='c')
    else:
        Vf = None

    scatter.destroy()
    Vn0.destroy()
    return Vf

#
#==================== Data input ============================================
#

class input(object):
    '''
    Hold data that will be used 
    Interpolate data in time
    '''

    def __init__(self,variable,files, da):
        ''' init data input
        should test if variables are 2D
        '''
        # browses files in order to figure out available times
        self.data = da.createGlobalVec()
        pass
    
    
    def update(self,time):
        ''' interpolate input data at a given time
        '''
        pass








