import sys
import xarray as xr
from netCDF4 import Dataset
import numpy as np

# ------------------------------------------------------------------------------------------------

# useful parameters
g=9.81


def fix_coords(ds):
    ''' fix coordinates in ROMS output files
    '''
    # indices are wrong in netcdf files
    ds = (ds.reset_index(['xi_rho','eta_rho','xi_u','eta_v'])
                 .drop(['xi_rho_','eta_rho_','xi_u_','eta_v_']))

    xrho = ds.x_rho.isel(eta_rho=0).drop('y_rho').rename('xi_rho')
    yr = ds.y_rho.isel(xi_rho=0).drop('x_rho').rename('eta_rho')

    xu = ( ((xrho.shift(xi_rho=1)+xrho)*0.5)
          .isel(xi_rho=slice(1,None))
          .rename({'xi_rho':'xi_u'}) )

    yv = ( ((yr.shift(eta_rho=1)+yr)*0.5)
          .isel(eta_rho=slice(1,None))
          .rename({'eta_rho':'eta_v'}) )

    # cannot add xi_psi and eta_psi, they are not dimensions
    #xp = xu.rename({'xi_u': 'xi_psi'})
    #yp = yv.rename({'eta_v': 'eta_psi'})

    return ds.assign_coords(xi_rho=xrho,eta_rho=yr, xi_u=xu, eta_v=yv)


@xr.register_dataset_accessor('grd')
class datasetgrd(object):
    ''' add grd accessor in order to perform common manipulations on ROMS outputs
    '''
    def __init__(self, xarray_obj):
        self._obj = xarray_obj
        # does not work from here
        #self.fix_coords()
        
    def rho2u(self, vin):
        ds = self._obj
        if isinstance(vin, str):
            v = self._obj[vin]
        else:
            v = vin
        vout = (v.shift(xi_rho=1)+v)*0.5
        vout = vout.isel(xi_rho=slice(1,None))
        vout = vout.rename({'xi_rho':'xi_u'}).assign_coords(xi_u=ds.xi_u)
        return vout

    def rho2v(self, vin):
        ds = self._obj
        if isinstance(vin, str):
            v = self._obj[vin]
        else:
            v = vin
        vout = (v.shift(eta_rho=1)+v)*0.5
        vout = vout.isel(eta_rho=slice(1,None))
        vout = vout.rename({'eta_rho':'eta_v'}).assign_coords(eta_v=ds.eta_v)
        return vout

    def rho2psi(self, vin):
        ds = self._obj
        if isinstance(vin, str):
            v = self._obj[vin]
        else:
            v = vin
        vpsi = (v.shift(xi_rho=1)+v+v.shift(xi_rho=1,eta_rho=1)+v.shift(eta_rho=1))*0.25
        vpsi = vpsi.isel(xi_rho=slice(1,None), eta_rho=slice(1,None))
        vpsi = vpsi.rename({'xi_rho':'xi_psi', 'eta_rho':'eta_psi'})
        vpsi = vpsi.assign_coords(xi_psi=ds.xi_psi, eta_psi=ds.eta_psi)
        return vpsi

    def psi2rho(self, vin):
        ds = self._obj
        if isinstance(vin, str):
            v = self._obj[vin]
        else:
            v = vin
        vout = (v.shift(xi_psi=1)+v+v.shift(xi_psi=1,eta_psi=1)+v.shift(eta_psi=1))*0.25
        vout = vout.isel(xi_psi=slice(1,None), eta_psi=slice(1,None))
        vout = vout.rename({'xi_psi':'xi_rho', 'eta_psi':'eta_rho'})
        vout = vout.assign_coords(xi_rho=ds.xi_rho[1:-1], eta_rho=ds.eta_rho[1:-1])
        vout, _ = xr.align(vout, ds.h, join='outer')
        return vout

    def get_z(self, zeta=None, h=None, vgrid='r', hgrid='rho'):
        ''' compute vertical coordinates
            zeta should have the size of the final output
            vertical coordinate should be first
        '''
        ds = self._obj
        if zeta is not None:
            _zeta = zeta
        else:
            _zeta = ds.zeta
        if h is not None:
            _h = h
        else:
            _h = ds.h
        #
        if hgrid is 'u':
            _zeta = self.rho2u(_zeta)
            _h = self.rho2u(_h)
        elif hgrid is 'v':
            _zeta = self.rho2v(_zeta)
            _h = self.rho2v(_h)
        #
        sc=ds['sc_'+vgrid]
        cs=ds['Cs_'+vgrid]
        #
        z0 = (ds.hc * sc + _h * cs) / (ds.hc + _h)
        z = _zeta + (_zeta + _h) * z0
        # manually swap dims, could also perform align with T,S
        if z.ndim == 3:
            z = z.transpose(sc.dims[0], _zeta.dims[0], _zeta.dims[1])
        elif z.ndim == 2:
            z = z.transpose(sc.dims[0], _zeta.dims[0])            
        return z.rename('z_'+vgrid)

    def get_uv_from_psi(self, psi):
        # note that u, v are computed at rho points
        ds = self._obj
        x, y = ds.xi_rho, ds.eta_rho
        #
        u = - 0.5*(psi.shift(eta_rho=1)-psi)/(y.shift(eta_rho=1)-y) \
            - 0.5*(psi-psi.shift(eta_rho=-1))/(y-y.shift(eta_rho=-1))
        #
        v =   0.5*(psi.shift(xi_rho=1)-psi)/(x.shift(xi_rho=1)-x) \
            + 0.5*(psi-psi.shift(xi_rho=-1))/(x-x.shift(xi_rho=-1))
        return u, v
    
def get_p(rho, zeta, rho0, rho_a=None):
    #
    if rho_a is None:
        _rho_a = 0.
    else:
        _rho_a = rho_a
    #
    _rho = g*(rho.shift(z_r=1)+rho)*.5*(rho.z_r.shift(z_r=1)-rho.z_r)
    _rho = _rho.sortby(_rho.z_r, ascending=False).shift(z_r=1).fillna(0.)
    p0 = g*(rho0+_rho_a.sel(z_r=0,method='nearest')+rho.sel(z_r=0,method='nearest'))*zeta
    p = _rho.cumsum('z_r').sortby(_rho.z_r, ascending=True) + p0
    p = p.rename('p')
    return p

def get_pv(u, v, rho, rho_a, f, f0, zr, zw, ds):

    # relative vorticity
    xi_v =  ( (v.diff('xi_rho')/ds.xi_rho.diff('xi_rho'))
             .rename({'xi_rho':'xi_psi', 'eta_v':'eta_psi'})
             .assign_coords(xi_psi=ds.xi_u.rename({'xi_u':'xi_psi'}), 
                            eta_psi=ds.eta_v.rename({'eta_v':'eta_psi'})))

    xi_u = -( (u.diff('eta_rho')/ds.eta_rho.diff('eta_rho'))
             .rename({'eta_rho':'eta_psi', 'xi_u':'xi_psi'})
             .assign_coords(xi_psi=ds.xi_u.rename({'xi_u':'xi_psi'}), 
                            eta_psi=ds.eta_v.rename({'eta_v':'eta_psi'})))
    xi = ds.grd.psi2rho(xi_u+xi_v)

    # stretching term
    drho = (rho.shift(z_r=1)+rho)*0.5 * zr.diff('z_r')/rho_a.diff('z_r')
    drho = drho.rename({'z_r':'z_w'}).assign_coords(z_w=zw[:-1])
    Sint = drho.diff('z_w')/zw[:-1].diff('z_w')
    Sint = Sint.rename({'z_w':'z_r'}).assign_coords(z_r=zr[1:-1])

    # note that top and bottom values are not used in the solver
    S = f0 * xr.concat([0.*rho.isel(z_r=0), Sint, 0.*rho.isel(z_r=-1)], dim='z_r') #.transpose('z_r','eta_rho','xi_rho')

    # assemble pb
    q = (xi + S + f - f0 ).rename('q') # order of variable conditions dimension order

    return q

def interp2z(z0, z, v, extrap):
    ''' Interpolate on a horizontally uniform grid
    '''
    import fast_interp3D as fi  # OpenMP accelerated C based interpolator
    #
    if v.ndim == 1 or z.ndim == 1 :
        lz = z.squeeze()[:,None,None]
        lv = v.squeeze()[:,None,None]
    elif v.ndim == 2 :
        lz = z[...,None]
        lv = v[...,None]
    else:
        lz = z[...]
        lv = v[...]
    #
    if extrap:
        zmin = np.min(z0)-1.
        lv = np.vstack((lv[[0],...], lv))
        lz = np.vstack((zmin+0.*lz[[0],...], lz))
    #
    vi = fi.interp(z0.astype('float64'), lz.astype('float64'), lv.astype('float64'))
    return vi

def interp2z_xr(z0, z, v,  hgrid='rho', extrap=True):
    _xmap = {'rho': 'rho', 'u': 'u', 'v': 'rho'}
    _ymap = {'rho': 'rho', 'u': 'rho', 'v': 'v'}
    _xgrid, _ygrid = _xmap[hgrid], _ymap[hgrid]
    return (xr.DataArray(interp2z(z0.values, z.values, v.values, extrap), 
                         dims=['z_r','eta_'+_ygrid,'xi_'+_xgrid],
                         coords={'z_r': z0.values, 'xi_'+_xgrid: v['xi_'+_xgrid], 
                                 'eta_'+_ygrid: v['eta_'+_ygrid]}) )


# ------------------------------------------------------------------------------------------------

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