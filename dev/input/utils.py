import sys
# netcdf
from netCDF4 import Dataset



def create_nc(filename, lon, lat, zt, zw):
    ''' Create an input netcdf file
    
    Parameters
    ----------
    filename: str
        name of the file created
    lon: ndarray
        2D array representing longitude ( lat x lon )
    lat: ndarray
        2D array representing latitude ( lat x lon )
    zt: 
    
    '''
    
    ### create a netcdf file
    rootgrp = Dataset(filename, 'w',
                      format='NETCDF4_CLASSIC', clobber=True)

    # create dimensions
    rootgrp.createDimension('x', lon.shape[1])
    rootgrp.createDimension('y', lat.shape[0])
    rootgrp.createDimension('zt', zt.size)
    rootgrp.createDimension('zw', zw.size)
    
    # create variables
    dtype='f8'
    nc_lon = rootgrp.createVariable('lon',dtype,('y','x'))
    nc_lat = rootgrp.createVariable('lat',dtype,('y','x'))
    nc_zt = rootgrp.createVariable('zt',dtype,('zt'))
    nc_zw = rootgrp.createVariable('zw',dtype,('zw'))
    
    nc_lon[:] = lon
    nc_lat[:] = lat
    nc_zt[:] = zt
    nc_zw[:] = zw
        
    return rootgrp


def query_yes_no(question, default="yes"):
    """Ask a yes/no question via raw_input() and return their answer.

    "question" is a string that is presented to the user.
    "default" is the presumed answer if the user just hits <Enter>.
        It must be "yes" (the default), "no" or None (meaning
        an answer is required of the user).

    The "answer" return value is True for "yes" or False for "no".
    """
    valid = {"yes": True, "y": True, "ye": True,
             "no": False, "n": False}
    if default is None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)

    while True:
        sys.stdout.write(question + prompt)
        choice = input().lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' "
                             "(or 'y' or 'n').\n")