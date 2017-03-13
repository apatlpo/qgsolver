#!/usr/bin/python
# -*- encoding: utf8 -*-

""" Create metrics, psi, pv fields for a curvlinear grid
    Fields are created from original (U,V,T,S,SSH) Arakawa-C grid NEMO fields.
    Density is computed with JMD95.

    Processing take advantage of oocgcm and natl60lib libraries, using xarray
    (chunked data). Most of the variables will be xarray DataArrays. One should
    get its values attribute to visualize data.

    The routine allows the processing of successive dates if user asks for several dates
    in call functions. Date will be then an extra dimension of all arrays.

    8/3/16: output writing deals only with numpy arrays, since we must reverse
            vertical coordinates (reverse might be possible with xarray arrays.
            (testing needed).
"""

import sys, os
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from netCDF4 import Dataset

# xarray 0.8.2
import xarray as xr

# natl60lib
import natl60lib.initlib_natl60 as init
import natl60lib.iolib_natl60 as iolib
import natl60lib.phylib_natl60_new_clean as phy

# oocgcm
import oocgcm.parameters.physicalparameters as oopp

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

    stest='mode2_fcFalse_fvertTrue'

    # - Control parameters

    # region
    region='LMX'

    # date (day and mth may be set to None
    day='01'
    mth='01'
    yr='2007'

    # manual chunk (set fchunk to True or change fvert)
    fchunk=False
    chunksx=None
    chunksy=None
    chunksz=None
    fvert=True

    # filtering options for pressure calculation
    # filter is defined in natl60lib.phylib_natl60 (Lanczos, width=6pts)
    dfilt=False    #  density filtering (input of pressure calculation)
    pfilt=False    #  pressure filtering (output pressure)

    # - Main

    iolib.mkdirp('figs/'+stest)
    iolib.mkdirp('data/'+stest)

    # Container object, to carry variables of the NEMO processing.
    # Chunking is chosen to preserve full vertical columns, in order 
    # to facilitate vertical differencing.
    xp=init.container(region,fvert)
    
    # Manual chunks: if user wants to test other chunking options
    # (time chunk is hard-coded as 1, since QG solver does not need time differencing).
    if fchunk:
       xp.init_chunks(xp.region,chunksy=chunksy,chunksx=chunksx,chunksz=chunksz)

    ### NEMO dir
    datadir=xp.natl60_path

    ### horizontal grid 
    lon=xp.grd._arrays['longitude_at_t_location']
    lat=xp.grd._arrays['latitude_at_t_location']
    e1=xp.grd._arrays['cell_x_size_at_t_location'] 
    e2=xp.grd._arrays['cell_y_size_at_t_location'] 
   
    vlon=lon.to_masked_array()
    vlat=lat.to_masked_array()
    ve1=e1.to_masked_array()
    ve2=e2.to_masked_array()
     
    # plot hgrid
    lims=[-80,-55, 30, 45]
    lon_tcks = range(-80,-55, 5)
    lat_tcks = range(30,45,5)
    #
    plt.figure(figsize=(8,3))
    ax=plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent(lims,ccrs.PlateCarree())
    ax.plot(vlon[::10,::10],vlat[::10,::10],'k-',transform=ccrs.PlateCarree())
    ax.plot(vlon.transpose()[::10,::10],vlat.transpose()[::10,::10],'k-',transform=ccrs.PlateCarree())
    plt.title('Horizontal grid (every 10 points)', size=10) # to modify the title
    ax.set_xticks(lon_tcks, crs=ccrs.PlateCarree())
    ax.set_yticks(lat_tcks, crs=ccrs.PlateCarree())
    ax.coastlines(resolution='50m') # Currently can be one of “110m”, “50m”, and “10m”
    ax.gridlines()
    #figname='figs/'+stest+'/snapshot_'+vkey.replace (" ", "_")+'_magnitude.jpg'
    plt.savefig('figs/'+stest+'/nemo_input_hgrid.jpg', dpi=300)
    #print figname+' printed'
    
    ### vertical grid
    ### NATL60 grid is homogeneous with depth: metric coefficients are stored as 1-D. 
    ### (but could be different with other NEMO simulations)
    zc_xr = -xp.vgrd._arrays['depth_at_t_location'].isel(t=0)
    zc = zc_xr.to_masked_array()
    zc = zc[::-1]
    zf_xr = -xp.vgrd._arrays['depth_at_w_location'].isel(t=0)
    zf = zf_xr.to_masked_array()
    zf = zf[::-1]
 
    # e3 is defined positive in NEMO
    e3_xr = xp.vgrd._arrays['cell_z_size_at_t_location'].isel(x=0,y=0)
    e3 = e3_xr.to_masked_array()
    dzc_xr = xp.vgrd._arrays['cell_z_size_at_w_location'].isel(x=0,y=0)
    dzc = dzc_xr.to_masked_array()

    # plot vertical grid
    plt.figure()
    plt.plot(zf,'k+')
    plt.plot(zf,'k.')
    plt.grid()
    plt.savefig('figs/'+stest+'/nemo_input_vgrid.jpg', dpi=300)
        
    # plot horizontal metrics
    plt.figure(figsize=(8,3))
    ax=plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent(lims,ccrs.PlateCarree())
    im = ax.pcolormesh(vlon,vlat,ve1,transform=ccrs.PlateCarree())
    cbar = plt.colorbar(im, format="%.1f")
    plt.title('e1 [m]', size=10) # to modify the title
    ax.set_xticks(lon_tcks, crs=ccrs.PlateCarree())
    ax.set_yticks(lat_tcks, crs=ccrs.PlateCarree())
    ax.coastlines(resolution='50m') # Currently can be one of “110m”, “50m”, and “10m”
    ax.gridlines()    
    plt.savefig('figs/'+stest+'/nemo_input_e1.jpg', dpi=300)

    # plot horizontal metrics
    plt.figure(figsize=(8,3))
    ax=plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent(lims,ccrs.PlateCarree())
    im = ax.pcolormesh(vlon,vlat,ve2,transform=ccrs.PlateCarree())
    cbar = plt.colorbar(im, format="%.1f")
    plt.title('e2 [m]', size=10) # to modify the title
    ax.set_xticks(lon_tcks, crs=ccrs.PlateCarree())
    ax.set_yticks(lat_tcks, crs=ccrs.PlateCarree())
    ax.coastlines(resolution='50m') # Currently can be one of “110m”, “50m”, and “10m”
    ax.gridlines()    
    plt.savefig('figs/'+stest+'/nemo_input_e2.jpg', dpi=300)

        
    # store metric terms
    #zc = np.hstack((zc,zc[[-1]]))
    rootgrp = create_nc('data/'+stest+'/nemo_metrics.nc', vlon, vlat, zc, zf)
    #
    dtype='f8'
    #nc_zf = rootgrp.createVariable('zf',dtype,('zf'))
    #nc_zf[:] = zf
    # 
    nc_e1 = rootgrp.createVariable('e1',dtype,('y','x'))
    nc_e1[:] = ve1
    nc_e2 = rootgrp.createVariable('e2',dtype,('y','x'))
    nc_e2[:] = ve2
    #nc_e3 = rootgrp.createVariable('e3',dtype,('zc'))
    ##nc_e3[:] = e3.append(e3[-1])
    #nc_e3[:] = np.hstack((e3,e3[[-1]]))
    
    rootgrp.close()
        
    ###
    
    # compute the Coriolis frequency and a reference value
    # from oocgcm/oocgcm/parameters/physicalparameters.py
    grav = oopp.grav              # acceleration due to gravity (m.s-2)
    omega = oopp.omega            # earth rotation rate (s-1)
    earthrad = oopp.earthrad      # mean earth radius (m)
    f = xp.grd._arrays['coriolis_parameter_at_t_location'].to_masked_array()
    f0 = phy.average_coriolis_parameter2(xp.region)
    
    # load stratification profile
    N2_file = datadir+'DIAG_DIMUP/2007_2008/'+xp.region+'/bg/'+xp.region\
       +'_2007_2008_bruntvaisala_bg_mindepth10_new.nc'
    rhobg_file = datadir+'DIAG_DIMUP/2007_2008/'+xp.region+'/bg/'+xp.region\
       +'_2007_2008_density_bg_mindepth10.nc'
    rhobg_xr,N2_xr=phy.call_background(xp,xp.region,xp.grd,xp.vgrd,z=None,mindepth=10,
           N2file=N2_file,rhobgfile=rhobg_file)
    #N2[1:] = N2[:-1]
    print N2_file
    print N2_xr
    print rhobg_file
    print rhobg_xr
    N2=N2_xr.to_masked_array()
    rhobg = rhobg_xr.to_masked_array()

    np.hstack((N2[0],N2))
    #np.hstack((N2,N2[-1]))
    print '!!! Shape of N2 is %d' %N2.shape

    # compute density, psi, PV
    rho_xr,psi_xr,q_xr= phy.call_qgpv(xp,day,mth,yr,z=None,fcdens=False,fccurl=None,fcstretch=False,
    #q_xr= phy.call_qgpv(xp,day,mth,yr,z=None,fcdens=False,fccurl=False,fcstretch=False,
        dfilt=dfilt,pfilt=pfilt,mode=2,fout=True,N2file=N2_file,rhobgfile=rhobg_file) 
    print rho_xr,psi_xr,q_xr
    rho = rho_xr.to_masked_array()
    psi = psi_xr.to_masked_array()
    q = q_xr.to_masked_array()


    # store variables
    rootgrp = create_nc('data/'+stest+'/nemo_pv.nc', vlon, vlat, zc, zf)
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
    im = ax.pcolormesh(vlon,vlat,q[5,:,:]/f0,transform=ccrs.PlateCarree())
    cbar = plt.colorbar(im, format="%.2f")
    plt.title('q/f0(z=%0.0f) [1]' %zc[-5-1], size=10) # to modify the title
    ax.set_xticks(lon_tcks, crs=ccrs.PlateCarree())
    ax.set_yticks(lat_tcks, crs=ccrs.PlateCarree())
    ax.coastlines(resolution='50m') # Currently can be one of “110m”, “50m”, and “10m”
    ax.gridlines()
    plt.savefig('figs/'+stest+'/nemo_input_q.jpg', dpi=300)

    rootgrp = create_nc('data/'+stest+'/nemo_psi.nc', vlon, vlat, zc, zf)
    #
    nc_psi = rootgrp.createVariable('psi',dtype,('zc','y','x'))
    for k in xrange(zc.size):
        nc_psi[k,...] = psi[psi.shape[0]-1-k,...]  

    rootgrp.close()
  
    # plot psi
    plt.figure(figsize=(8,3))
    ax=plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent(lims,ccrs.PlateCarree())
    im = ax.pcolormesh(vlon,vlat,psi[5,:,:],transform=ccrs.PlateCarree())
    cbar = plt.colorbar(im, format="%.2f")
    plt.title('psi(z=%0.0f) [m^2/s]' %zc[-5-1], size=10) # to modify the title
    ax.set_xticks(lon_tcks, crs=ccrs.PlateCarree())
    ax.set_yticks(lat_tcks, crs=ccrs.PlateCarree())
    ax.coastlines(resolution='50m') # Currently can be one of “110m”, “50m”, and “10m”
    ax.gridlines()
    plt.savefig('figs/'+stest+'/nemo_input_psi.jpg', dpi=300)
    

    ### load and store rho
    rho_xr = arho_xr*oopp.rau0+rhobg_xr
    rho=rho_xr.to_masked_array()
    
    # load background density
    
    rootgrp = create_nc('data/'+stest+'/nemo_rho.nc', vlon, vlat, zc, zf)
    #
    nc_rho = rootgrp.createVariable('rho',dtype,('zc','y','x'))
    for k in xrange(zc.size):
        nc_rho[k,...] = rho[rho.shape[0]-1-k,...] - rhobg[rho.shape[0]-1-k,None,None]
    
    rootgrp.close()

    # plot rho                      
    plt.figure(figsize=(8,3))
    ax=plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent(lims,ccrs.PlateCarree())
    im = ax.pcolormesh(vlon,vlat,vrho[5,:,:]-vrhobg[5],transform=ccrs.PlateCarree())
    cbar = plt.colorbar(im, format="%.2f")
    plt.title('rho(z=%0.0f) [kg/m^3]' %zc[-5-1], size=10) # to modify the title
    ax.set_xticks(lon_tcks, crs=ccrs.PlateCarree())
    ax.set_yticks(lat_tcks, crs=ccrs.PlateCarree())
    ax.coastlines(resolution='50m') # Currently can be one of “110m”, “50m”, and “10m”
    ax.gridlines()
    plt.savefig('figs/'+stest+'/nemo_input_rho.jpg', dpi=300)
    
    plt.ion()
    plt.show(); 
