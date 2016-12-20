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
    datadir='/home7/pharos/othr/NATL60/'
    griddir=datadir+'NATL60-I/BOXES/'
    
    
    ### horizontal grid
    hgrid_file=griddir+'NATL60LMX_coordinates_v4.nc'
    
    hgrid = Dataset(hgrid_file, 'r')
    lon = hgrid.variables['nav_lon'][:]
    lat = hgrid.variables['nav_lat'][:]
    e1 = hgrid.variables['e1t'][:]
    e2 = hgrid.variables['e2t'][:]
    
    
    # plot hgrid
    lims=[-80,-55, 30, 45]
    lon_tcks = range(-80,-55, 5)
    lat_tcks = range(30,45,5)
    #
    plt.figure(figsize=(8,3))
    ax=plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent(lims,ccrs.PlateCarree())
    ax.plot(lon[::10,::10],lat[::10,::10],'k-',transform=ccrs.PlateCarree())
    ax.plot(lon.transpose()[::10,::10],lat.transpose()[::10,::10],'k-',transform=ccrs.PlateCarree())
    plt.title('Horizontal grid (every 10 points)', size=10) # to modify the title
    ax.set_xticks(lon_tcks, crs=ccrs.PlateCarree())
    ax.set_yticks(lat_tcks, crs=ccrs.PlateCarree())
    ax.coastlines(resolution='50m') # Currently can be one of “110m”, “50m”, and “10m”
    ax.gridlines()
    #figname='figs/snapshot_'+vkey.replace (" ", "_")+'_magnitude.jpg'
    plt.savefig('figs/nemo_input_hgrid.jpg', dpi=300)
    #print figname+' printed'
    
    
    ### vertical grid
    vgrid_file=griddir+'NATL60LMX_v4.1_cdf_mesh_zgr.nc'
    
    vgrid = Dataset(vgrid_file, 'r')
    zc = -vgrid.variables['gdept_0'][0,::-1]
    zf = -vgrid.variables['gdepw_0'][0,::-1]    
    e3 = np.diff(zf)
    

    # plot vertical grid
    plt.figure()
    plt.plot(zf,'k+')
    plt.plot(zf,'k.')
    plt.grid()
    plt.savefig('figs/nemo_input_vgrid.jpg', dpi=300)
        
    # plot horizontal metrics
    plt.figure(figsize=(8,3))
    ax=plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent(lims,ccrs.PlateCarree())
    im = ax.pcolormesh(lon,lat,e1,transform=ccrs.PlateCarree())
    cbar = plt.colorbar(im, format="%.1f")
    plt.title('e1 [m]', size=10) # to modify the title
    ax.set_xticks(lon_tcks, crs=ccrs.PlateCarree())
    ax.set_yticks(lat_tcks, crs=ccrs.PlateCarree())
    ax.coastlines(resolution='50m') # Currently can be one of “110m”, “50m”, and “10m”
    ax.gridlines()    
    plt.savefig('figs/nemo_input_e1.jpg', dpi=300)

    # plot horizontal metrics
    plt.figure(figsize=(8,3))
    ax=plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent(lims,ccrs.PlateCarree())
    im = ax.pcolormesh(lon,lat,e2,transform=ccrs.PlateCarree())
    cbar = plt.colorbar(im, format="%.1f")
    plt.title('e2 [m]', size=10) # to modify the title
    ax.set_xticks(lon_tcks, crs=ccrs.PlateCarree())
    ax.set_yticks(lat_tcks, crs=ccrs.PlateCarree())
    ax.coastlines(resolution='50m') # Currently can be one of “110m”, “50m”, and “10m”
    ax.gridlines()    
    plt.savefig('figs/nemo_input_e2.jpg', dpi=300)

        
    # store metric terms
    #zc = np.hstack((zc,zc[[-1]]))
    rootgrp = create_nc('data/nemo_metrics.nc', lon, lat, zc, zf)
    #
    dtype='f8'
    #nc_zf = rootgrp.createVariable('zf',dtype,('zf'))
    #nc_zf[:] = zf
    # 
    nc_e1 = rootgrp.createVariable('e1',dtype,('y','x'))
    nc_e1[:] = e1
    nc_e2 = rootgrp.createVariable('e2',dtype,('y','x'))
    nc_e2[:] = e2
    #nc_e3 = rootgrp.createVariable('e3',dtype,('zc'))
    ##nc_e3[:] = e3.append(e3[-1])
    #nc_e3[:] = np.hstack((e3,e3[[-1]]))
    
    rootgrp.close()
        
    
    
    
    ###
    
    # compute the Coriolis frequency and a reference value
    # from oocgcm/oocgcm/parameters/physicalparameters.py
    grav = 9.81                  # acceleration due to gravity (m.s-2)
    omega = 7.292115083046061e-5 # earth rotation rate (s-1)
    earthrad = 6371229            # mean earth radius (m)
    f = 2. * omega * np.sin(lat * d2r)
    f0 = np.mean(f)
    
    # load stratification profile
    N2_file = datadir+'DIAG_DIMUP/2007_2008/LMX/bg/LMX_2007_2008_bruntvaisala_bg_mindepth10.nc'
    nc = Dataset(N2_file, 'r')
    N2 = nc.variables['bvfreq_bg'][:]
    nc.close()
    #N2[1:] = N2[:-1]
    np.hstack((N2[0],N2))
    #np.hstack((N2,N2[-1]))
    print '!!! Shape of N2 is %d' %N2.shape

    # load PV
    pv_file = datadir+'DIAG_DIMUP/qgpv/LMX/test_good/LMX_y2007m01d01_qgpv_v1.nc'
    nc = Dataset(pv_file, 'r')
    q = nc.variables['qgpv_v1']


    # store variables
    rootgrp = create_nc('data/nemo_pv.nc', lon, lat, zc, zf)
    #
    nc_f = rootgrp.createVariable('f',dtype,('y','x'))
    nc_f[:] = f    
    nc_f0 = rootgrp.createVariable('f0',dtype)
    nc_f0[:] = f0
    #
    nc_N2 = rootgrp.createVariable('N2',dtype,('zf'))
    nc_N2[:] = N2[::-1]
    #
    nc_q = rootgrp.createVariable('q',dtype,('zc','y','x'))
    for k in xrange(zc.size):
        nc_q[k,...] = q[q.shape[0]-1-k,...]
    #
    rootgrp.close()


    # plot pv
    plt.figure(figsize=(8,3))
    ax=plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent(lims,ccrs.PlateCarree())
    im = ax.pcolormesh(lon,lat,q[5,:,:]/f0,transform=ccrs.PlateCarree())
    cbar = plt.colorbar(im, format="%.2f")
    plt.title('q/f0(z=%0.0f) [1]' %zc[-5-1], size=10) # to modify the title
    ax.set_xticks(lon_tcks, crs=ccrs.PlateCarree())
    ax.set_yticks(lat_tcks, crs=ccrs.PlateCarree())
    ax.coastlines(resolution='50m') # Currently can be one of “110m”, “50m”, and “10m”
    ax.gridlines()
    plt.savefig('figs/nemo_input_q.jpg', dpi=300)

    # close file
    nc.close()
    

  
    ### load and store psi
    
    psi_file = datadir+'DIAG_DIMUP/psi0/LMX/LMX_y2007m01d01_psi0.nc'
    nc = Dataset(psi_file, 'r')
    psi = nc.variables['psi0'][:]
    
    rootgrp = create_nc('data/nemo_psi.nc', lon, lat, zc, zf)
    #
    nc_psi = rootgrp.createVariable('psi',dtype,('zc','y','x'))
    for k in xrange(zc.size):
        nc_psi[k,...] = psi[psi.shape[0]-1-k,...]    
    
  
    # plot psi
    plt.figure(figsize=(8,3))
    ax=plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent(lims,ccrs.PlateCarree())
    im = ax.pcolormesh(lon,lat,psi[5,:,:],transform=ccrs.PlateCarree())
    cbar = plt.colorbar(im, format="%.2f")
    plt.title('psi(z=%0.0f) [m^2/s]' %zc[-5-1], size=10) # to modify the title
    ax.set_xticks(lon_tcks, crs=ccrs.PlateCarree())
    ax.set_yticks(lat_tcks, crs=ccrs.PlateCarree())
    ax.coastlines(resolution='50m') # Currently can be one of “110m”, “50m”, and “10m”
    ax.gridlines()
    plt.savefig('figs/nemo_input_psi.jpg', dpi=300)
    
    # close file
    nc.close()    



    ### load and store rho
    
    rho_file = datadir+'DIAG_DIMUP/density/LMX/LMX_y2007m01d01_density.nc'
    nc = Dataset(rho_file, 'r')
    rho = nc.variables['density'][:]
    
    # load background density
    rhobg_file = datadir+'DIAG_DIMUP/2007_2008/LMX/bg/LMX_2007_2008_density_bg_mindepth10.nc'
    ncbg = Dataset(rhobg_file, 'r')
    rhobg = ncbg.variables['density_bg'][:]
    
    
    
    rootgrp = create_nc('data/nemo_rho.nc', lon, lat, zc, zf)
    #
    nc_rho = rootgrp.createVariable('rho',dtype,('zc','y','x'))
    for k in xrange(zc.size):
        nc_rho[k,...] = rho[rho.shape[0]-1-k,...] - rhobg[rho.shape[0]-1-k,None,None]
    
    # plot horizontal metrics
    plt.figure(figsize=(8,3))
    ax=plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent(lims,ccrs.PlateCarree())
    im = ax.pcolormesh(lon,lat,rho[5,:,:]-rhobg[5],transform=ccrs.PlateCarree())
    cbar = plt.colorbar(im, format="%.2f")
    plt.title('rho(z=%0.0f) [kg/m^3]' %zc[-5-1], size=10) # to modify the title
    ax.set_xticks(lon_tcks, crs=ccrs.PlateCarree())
    ax.set_yticks(lat_tcks, crs=ccrs.PlateCarree())
    ax.coastlines(resolution='50m') # Currently can be one of “110m”, “50m”, and “10m”
    ax.gridlines()
    plt.savefig('figs/nemo_input_rho.jpg', dpi=300)
    
    # close file
    nc.close()    


    
    plt.ion()
    plt.show(); 