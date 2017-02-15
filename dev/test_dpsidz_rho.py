#!/usr/bin/python
# -*- encoding: utf8 -*-

""" Create metrics, psi, pv fields for a curvlinear grid
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from netCDF4 import Dataset

# maybe temporary
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.fftpack._fftpack import zfft

d2r = np.pi/180.
reflev=250

def create_nc(filename, lon, lat, zc, zf):
    
    ### create a netcdf file
    rootgrp = Dataset(filename, 'w',
                      format='NETCDF4_CLASSIC', clobber=True)

    # create dimensions
    rootgrp.createDimension('x', lon.shape[1])
    rootgrp.createDimension('y', lat.shape[0])
    rootgrp.createDimension('zc', zc.size)
    rootgrp.createDimension('zf', zf.size)
    
    # create variables
    dtype='f8'
    nc_lon = rootgrp.createVariable('lon',dtype,('y','x'))
    nc_lat = rootgrp.createVariable('lat',dtype,('y','x'))
    nc_zc = rootgrp.createVariable('zc',dtype,('zc'))
    nc_zf = rootgrp.createVariable('zf',dtype,('zf'))
    
    nc_lon[:] = lon
    nc_lat[:] = lat
    nc_zc[:] = zc
    nc_zf[:] = zf
        
    #rootgrp.createVariable(name,dtype,('z','y','x',)))
    return rootgrp



if __name__ == "__main__":
    
    
    ### NEMO grid file
    # datadir='/home7/pharos/othr/NATL60/'
    datadir='data/'
    griddir=datadir+'NATL60-I/BOXES/'
    
    
    ### horizontal grid
    hgrid_file=griddir+'NATL60LMX_coordinates_v4.nc'
    
    hgrid = Dataset(hgrid_file, 'r')
    lon = hgrid.variables['nav_lon'][:]
    lat = hgrid.variables['nav_lat'][:]
    e1 = hgrid.variables['e1t'][:]
    e2 = hgrid.variables['e2t'][:]
    
    
    ### vertical grid
    vgrid_file=griddir+'NATL60LMX_v4.1_cdf_mesh_zgr.nc'
    
    vgrid = Dataset(vgrid_file, 'r')
    zc = -vgrid.variables['gdept_0'][0,::-1]
    zf = -vgrid.variables['gdepw_0'][0,::-1]   
    myzc = -vgrid.variables['gdept_0'][0,:]
    myzf = -vgrid.variables['gdepw_0'][0,:]  
    e3 = np.diff(zf)
    dzc = np.diff(zc)
    dzf = np.diff(zf)
    print 'lon:',lon.shape
    print 'lat:',lat.shape
    print 'zc:',zc.shape, zc[0], zc[-1]
    print 'zf:',zf.shape

    
    
    # compute the Coriolis frequency and a reference value
    # from oocgcm/oocgcm/parameters/physicalparameters.py
    grav = 9.81                  # acceleration due to gravity (m.s-2)
    omega = 7.292115083046061e-5 # earth rotation rate (s-1)
    earthrad = 6371229            # mean earth radius (m)
    f = 2. * omega * np.sin(lat * d2r)
    f0 = np.mean(f)    
    rho0 = 1000.
    dtype='f8'


    ### load and store rho
    
    rho_file = datadir+'DIAG_DIMUP/density/LMX/LMX_y2007m01d01_density.nc'
    rhoin = Dataset(rho_file, 'r')
    rho = rhoin.variables['density'][:]
 

    # load background density
    rhobg_file = datadir+'DIAG_DIMUP/2007_2008/LMX/bg/LMX_2007_2008_density_bg_mindepth10.nc'
    bgin = Dataset(rhobg_file, 'r')
    rhobg = bgin.variables['density_bg'][:]
    print 'rhobg:',rhobg.shape


    rhoout = create_nc('data/nemo_rho.nc', lon, lat, zc, zf)
    #
    nc_rho = rhoout.createVariable('rho',dtype,('zc','y','x'))
    for k in xrange(zc.size):
        nc_rho[k,...] = rho[rho.shape[0]-1-k,...] - rhobg[rho.shape[0]-1-k,None,None]

    # create netcdf file with rho at level reflev        
    rhomy = create_nc('data/rho.nc', lon, lat, zc[reflev], zf[reflev])
    nc_rholev = rhomy.createVariable('rholev',dtype,('zc','y','x'))
    nc_rholev[0,:,:] = - grav*0.5*(nc_rho[reflev,:,:] + nc_rho[reflev+1,:,:])/(rho0*f0)
    # nc_rholev[0,:,:] = - grav*0.5*(nc_rho[reflev+1,:,:] + nc_rho[reflev+2,:,:])/(rho0*f0)
    # nc_rholev[0,:,:] = - grav*0.5*(nc_rho[reflev-1,:,:] + nc_rho[reflev,:,:])/(rho0*f0)
    # nc_rholev[0,:,:] = - grav*0.5*(nc_rho[reflev-2,:,:] + nc_rho[reflev-1,:,:])/(rho0*f0)



    ### load and store psi
    
    psi_file = datadir+'DIAG_DIMUP/psi0/LMX/LMX_y2007m01d01_psi0.nc'
    psiin = Dataset(psi_file, 'r')
    psi = psiin.variables['psi0'][:] 
    print 'psi:',psi.shape

    # compute my psi from rho
    mypsi = np.empty_like(rho) 
    mypsi[0,:,:] = -0.5*grav*(myzf[0]-myzc[0])*(rho[0,:,:]- rhobg[0])/rho0/f0
    for k in range(1,zc.size):
        mypsi[k,:,:]=mypsi[k-1,:,:] - 0.5*grav*(myzc[k-1]-myzc[k])*(rho[k-1,:,:]-rhobg[k-1]+ \
            rho[k,:,:]- rhobg[k])/rho0/f0

    # store psi
    psiout = create_nc('data/nemo_psi.nc', lon, lat, zc, zf)
    nc_psi = psiout.createVariable('psi',dtype,('zc','y','x'))
    for k in xrange(zc.size):
        nc_psi[k,...] = psi[psi.shape[0]-1-k,...]
        # nc_psi[k,...] = mypsi[psi.shape[0]-1-k,...]
        

    # create netcdf file with dpsi/dz at level reflev
    dpsidzmy = create_nc('data/dpsidz.nc', lon, lat, zc[reflev], zf[reflev])
    nc_dpsidz = dpsidzmy.createVariable('dpsidz',dtype,('zc','y','x'))
    nc_dpsidz[0,:,:] = (nc_psi[reflev+1,:,:]-nc_psi[reflev,:,:])/dzc[reflev]
  

    # create netcdf file with dpsi/dz - rho at level reflev
    diffmy = create_nc('data/diff.nc', lon, lat, zc[reflev], zf[reflev])
    diff = diffmy.createVariable('diff',dtype,('zc','y','x'))
    diff[:] = nc_dpsidz[:]- nc_rholev[:]
 
    # close file

    rhoin.close()  
    bgin.close()   
    psiin.close() 

    rhoout.close() 
    psiout.close() 

    rhomy.close() 
    diffmy.close()    
    dpsidzmy.close() 

