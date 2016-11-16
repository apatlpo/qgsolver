import sys
# netcdf
from netCDF4 import Dataset

def  outnc(varname,vardata):
    """
    :param varname: name of variable to save and of output file siffixed by .nc
    :param vardata: array to save
    :return:
    """
    ### create a netcdf file to store QG pv for inversion
    rootgrp = Dataset(varname+".nc", 'w', format='NETCDF4_CLASSIC', clobber=True)

    # create dimensions
    rootgrp.createDimension('x', vardata.shape[0])
    rootgrp.createDimension('y', vardata.shape[1])


    # create variables
    dtype='f8'
    nc_var = rootgrp.createVariable(varname,dtype,('x','y'))


    # fills in coordinate variables, keep truely periodic data
    nc_var[:]=vardata[:]

    # close the netcdf file
    rootgrp.close()

