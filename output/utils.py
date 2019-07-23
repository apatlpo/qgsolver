import sys
from glob import glob
import xarray as xr
from netCDF4 import Dataset


# ----------------------------------- xarray ------------------------------------

def open_input(path, substract_bg=True):
    D = []
    for f in glob(path+'input/*.nc'):
        if '_bg' in f:
            D.append(xr.open_dataset(f).rename({'psi': 'psi_bg', 'q': 'q_bg', 'rho': 'rho_bg'}))
        else:
            D.append(xr.open_dataset(f))

    ds = xr.merge(D)
    ds = ds.rename({'zt': 'z'})
    ds = ds.assign_coords(x_km=ds.x/1e3, y_km=ds.y/1e3)    
    if substract_bg:
        # substract bg from solution in order to compare it with outputs
        ds['psi'] = ds['psi'] - ds['psi_bg']
        ds['q'] = ds['q'] - ds['q_bg']
    # derive useful diagnostics
    ds['zeta'] = vorticity(ds.psi)
    ds['u'], ds['v'] = velocity(ds.psi)    
    return ds

def open_output(file_out, grd):
    ds = xr.open_dataset(file_out)
    ds = ds.drop(['x','y'])
    ds = ds.assign_coords(x=grd.x, y=grd.y)
    ds['zeta'] = vorticity(ds.psi)
    ds['u'], ds['v'] = velocity(ds.psi)
    ds = ds.assign_coords(x_km=ds.x/1e3, y_km=ds.y/1e3)
    return ds

def velocity(psi):
    x, y = psi.x, psi.y
    # assumes a uniform grid
    dx, dy = (x.shift(x=1)-x), (y.shift(y=1)-y)
    u = (psi.shift(x=1)-psi.shift(x=-1))/dx/2.
    v = (psi.shift(y=1)-psi.shift(y=-1))/dy/2.
    return u.rename('u'), v.rename('v')

def vorticity(psi):
    x, y = psi.x, psi.y
    # assumes a uniform grid
    dx, dy = (x.shift(x=1)-x), (y.shift(y=1)-y)
    xi = (psi.shift(x=1)-2.*psi+psi.shift(x=-1))/dx/dx \
         + (psi.shift(y=1)-2.*psi+psi.shift(y=-1))/dy/dy
    return xi.rename('zeta')

def xrsel(da, coord, x, threshold=1.):
    ''' select xarray object around coordinate that is not a dimension

    Parameters
    ----------
    da: xarray DataArray/Dataset
    coord: str
        coordinate used for selection
    x: float/int
        value selected
    threshold: float/int
        threshold value for selection, optional, default to 1
    '''
    return da.where( abs(da[coord]-x)<threshold, drop=True).isel(**{da[coord].dims[0]: 0})


# ----------------------------------- netcdf4 ------------------------------------

def outnc(varname,vardata):
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

