#!/usr/bin/python
# -*- encoding: utf8 -*-

""" Create metrics, psi, pv fields for a curvlinear grid
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import netCDF4
from netCDF4 import Dataset

# maybe temporary
# import matplotlib as mpl
# from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.fftpack._fftpack import zfft

d2r = np.pi/180.
fillvalue = netCDF4.default_fillvals['f8']

def create_nc(filename, lon, lat, zt, zw):
    
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
        
    # rootgrp.createVariable(name,dtype,('zt','y','x',)))
    return rootgrp



if __name__ == "__main__":
    
    # doesn't keep levels under the depth mask_depth
    mask_depth = -3000.
    # mask_depth = -7000.


    ### NEMO grid file
    datadir='/home7/pharos/othr/NATL60/'
    # datadir='data/'
    griddir=datadir+'NATL60-I/BOXES/'
    
    
    ### horizontal grid
    print "read horizontal grid"
    hgrid_file=griddir+'NATL60LMX_coordinates_v4.nc'
    
    hgrid = Dataset(hgrid_file, 'r')
    lon = hgrid.variables['nav_lon'][:]
    lat = hgrid.variables['nav_lat'][:]
    dxt = hgrid.variables['e1t'][:]
    dyt = hgrid.variables['e2t'][:]
    dxu = hgrid.variables['e1u'][:]
    dyu = hgrid.variables['e2u'][:]
    dxv = hgrid.variables['e1v'][:]
    dyv = hgrid.variables['e2v'][:]
    
    L = lon.shape[1]
    M = lat.shape[0]
    
    ### vertical grid, caution: level 0 in nemo correspond to the surface, level N correspond to positive depth
    print "read vertical grid"
    vgrid_file=griddir+'NATL60LMX_v4.1_cdf_mesh_zgr.nc'
    
    vgrid = Dataset(vgrid_file, 'r')
    zt = -vgrid.variables['gdept_0'][0,::-1]
    zw = -vgrid.variables['gdepw_0'][0,::-1]    
    N = zt.shape[0]
    dzt = vgrid.variables['e3t_0'][0,::-1]
    dzw = vgrid.variables['e3w_0'][0,::-1]

    # find nearest index in zt for mask_depth
    index_mask_depth = min(range(len(zt)), key=lambda i: abs(zt[i]-mask_depth))
    print "mask reference at level:",index_mask_depth

        
    # store metric terms
    #zt = np.hstack((zt,zt[[-1]]))
    print "create metrics grid"
    # metricsout = create_nc('data/nemo_metrics.nc', lon, lat, zt[index_mask_depth:], zw[index_mask_depth:])
    metricsout = create_nc('data/window_metrics.nc', lon, lat, zt, zw)
    #
    dtype='f8'
    # 
    nc_dxt = metricsout.createVariable('dxt',dtype,('y','x'))
    nc_dxt[:] = dxt
    nc_dyt = metricsout.createVariable('dyt',dtype,('y','x'))
    nc_dyt[:] = dyt
    nc_dxu = metricsout.createVariable('dxu',dtype,('y','x'))
    nc_dxu[:] = dxu
    nc_dyu = metricsout.createVariable('dyu',dtype,('y','x'))
    nc_dyu[:] = dyu
    nc_dxv = metricsout.createVariable('dxv',dtype,('y','x'))
    nc_dxv[:] = dxv
    nc_dyv = metricsout.createVariable('dyv',dtype,('y','x'))
    nc_dyv[:] = dyv
    nc_dzt = metricsout.createVariable('dzt',dtype,('zt'))
    nc_dzt[:] = dzt
    nc_dzw = metricsout.createVariable('dzw',dtype,('zw'))
    nc_dzw[:] = dzw

                
    ### load and store psi    
    print "read psi"
    psi_file = datadir+'DIAG_DIMUP/psi0/LMX/LMX_y2007m01d01_psi0_split.nc'
    psiin = Dataset(psi_file, 'r')
    psi = psiin.variables['psi0_hydrostatic']   
  
    #create 2D mask at reference level index_mask_depth (land=1, water=0)
    print "store mask"
    nc_mask = metricsout.createVariable('mask',dtype,('y','x'), fill_value=psi._FillValue)
    nc_mask[:] = psi[N-index_mask_depth-1,:,:]
    nc_mask[:] = np.where(nc_mask == nc_mask._FillValue, nc_mask, 1.) 
    nc_mask[:] = np.where(nc_mask != nc_mask._FillValue, nc_mask, 0.) 

    # enlarge the mask: if the i,j point has an adjacent land point then it becames land
    dummy = nc_mask[1:-1,1:-1]+nc_mask[:-2,1:-1]+nc_mask[2:,1:-1]+nc_mask[1:-1,:-2]+nc_mask[1:-1,2:]
    nc_mask[1:-1,1:-1] = np.where(dummy == 5, nc_mask[1:-1,1:-1], 0.)
    metricsout.close()
    psiin.close()



