#!/usr/bin/python
# -*- encoding: utf8 -*-

""" Read result from a PV inversion
"""

import sys
import numpy as np
from netCDF4 import Dataset

def create_nc(filename, x, y, z):
    
    ### create a netcdf file
    rootgrp = Dataset(filename, 'w',
                      format='NETCDF4_CLASSIC', clobber=True)

    # create dimensions
    rootgrp.createDimension('x', x.shape[0])
    rootgrp.createDimension('y', y.shape[0])
    rootgrp.createDimension('z', z.shape[0])
    
    # create variables
    dtype='f8'
    nc_x = rootgrp.createVariable('x',dtype,('x'))
    nc_y = rootgrp.createVariable('y',dtype,('y'))
    nc_z = rootgrp.createVariable('z',dtype,('z'))
    
    nc_x[:] = x
    nc_y[:] = y
    nc_z[:] = z
        
    return rootgrp


if __name__ == "__main__":

    d2r = np.pi/180.
    dtype='f8'
     

    if __name__ == "__main__":
        
        # data path
        # datadir='/home/caparmor-work/aponte/qg_nemo/dev/data/'
        # datadir='/home/caparmor-work/slgentil/nemotest/dev/data/'
        datadir='/home/mulroy/slgentil/models/natl60/qgsolver/dev/data/'


        # NEMO output file
        output_file=datadir+'output.nc'
        nco = Dataset(output_file, 'r')
        x = nco.variables['x']
        y = nco.variables['y']
        z = nco.variables['z']
        q_out = nco.variables['q']   
        psi_out = nco.variables['psi']   

        # NEMO input psi and q
        input_file=datadir+'input.nc'
        nci = Dataset(input_file, 'r')
        q_in = nci.variables['q']
        psi_in = nci.variables['psi']
        
        diffnc = create_nc('data/diff.nc', x, y, z)
        diff_q = diffnc.createVariable('diff_q',dtype,('z','y','x'))
        diff_psi = diffnc.createVariable('diff_psi',dtype,('z','y','x'))
        diff_q[:] = q_in[:] - q_out[:]
        diff_psi[:] = psi_in[:] - psi_out[:]

        diffnc.close()
        nco.close()
        nci.close()