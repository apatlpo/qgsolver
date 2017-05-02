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


import numpy as np
from netCDF4 import Dataset

# natl60lib
import natl60lib.initlib_natl60 as init
import natl60lib.iolib_natl60 as iolib
import natl60lib.phylib_natl60_new_clean as phy

# oocgcm
import oocgcm.parameters.physicalparameters as oopp

# maybe temporary
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import xarray as xr
import sys, os
#import matplotlib as mpl
#from mpl_toolkits.axes_grid1 import make_axes_locatable
#from scipy.fftpack._fftpack import zwft

d2r = np.pi/180.

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
        
    #rootgrp.createVariable(name,dtype,('z','y','x',)))
    return rootgrp



if __name__ == "__main__":

    stest='mode2_fcTrue_fvertTrue'

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

    # doesn't keep levels under the depth mask_depth
    mask_depth = -3000.
    # mask_depth = -7000.

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
    print "read horizontal grid"
    lon=xp.grd._arrays['longitude_at_t_location']
    lat=xp.grd._arrays['latitude_at_t_location']
    dxt=xp.grd._arrays['cell_x_size_at_t_location'].to_masked_array() 
    dyt=xp.grd._arrays['cell_y_size_at_t_location'].to_masked_array() 
    dxt = xp.grd._arrays['cell_x_size_at_t_location'].to_masked_array()
    dyt = xp.grd._arrays['cell_y_size_at_t_location'].to_masked_array()
    dxu = xp.grd._arrays['cell_x_size_at_u_location'].to_masked_array()
    dyu = xp.grd._arrays['cell_y_size_at_u_location'].to_masked_array()
    dxv = xp.grd._arrays['cell_x_size_at_v_location'].to_masked_array()
    dyv = xp.grd._arrays['cell_y_size_at_v_location'].to_masked_array()

    vlon=lon.to_masked_array()
    vlat=lat.to_masked_array()

    ## plot hgrid
    #lims=[-80,-55, 30, 45]
    #lon_tcks = range(-80,-55, 5)
    #lat_tcks = range(30,45,5)
    ##
    #plt.figure(figsize=(8,3))
    #ax=plt.axes(projection=ccrs.PlateCarree())
    #ax.set_extent(lims,ccrs.PlateCarree())
    #ax.plot(vlon[::10,::10],vlat[::10,::10],'k-',transform=ccrs.PlateCarree())
    #ax.plot(vlon.transpose()[::10,::10],vlat.transpose()[::10,::10],'k-',transform=ccrs.PlateCarree())
    #plt.title('Horizontal grid (every 10 points)', size=10) # to modify the title
    #ax.set_xticks(lon_tcks, crs=ccrs.PlateCarree())
    #ax.set_yticks(lat_tcks, crs=ccrs.PlateCarree())
    #ax.coastlines(resolution='50m') # Currently can be one of “110m”, “50m”, and “10m”
    #ax.gridlines()
    ##figname='figs/'+stest+'/snapshot_'+vkey.replace (" ", "_")+'_magnitude.jpg'
    #plt.savefig('figs/'+stest+'/nemo_input_hgrid.jpg', dpi=300)
    ##print figname+' printed'
    
    ### vertical grid, caution: level 0 in nemo correspond to the surface, level N correspond to positive depth
    print "read vertical grid"
    ### NATL60 grid is homogeneous with depth: metric coefficients are stored as 1-D. 
    ### (but could be different with other NEMO simulations)
    zt_xr = -xp.vgrd._arrays['depth_at_t_location'].isel(t=0)
    zti = zt_xr.to_masked_array()
    zt = np.flipud(zti)
    N = zt.shape[0]
    zw_xr = -xp.vgrd._arrays['depth_at_w_location'].isel(t=0)
    zwi = zw_xr.to_masked_array()
    zw = np.flipud(zwi)
 
    # dzt/dzw is defined positive in NEMO
    dzt_xr = xp.vgrd._arrays['cell_z_size_at_t_location'].isel(t=0,x=0,y=0)
    dzti = dzt_xr.to_masked_array()
    dzt = np.flipud(dzti)
    dzw_xr = xp.vgrd._arrays['cell_z_size_at_w_location'].isel(t=0,x=0,y=0)
    dzwi = dzw_xr.to_masked_array()
    dzw = np.flipud(dzwi)
    
    print zt
    print zw
    print dzt
    print dzw

    # find nearest index in zt for mask_depth
    index_mask_depth = min(range(len(zt)), key=lambda i: abs(zt[i]-mask_depth))
    print "mask reference at level:",index_mask_depth

    # plot vertical grid
    #plt.figure()
    #plt.plot(zw,'k+')
    #plt.plot(zw,'k.')
    #plt.grid()
    #plt.savefig('figs/'+stest+'/nemo_input_vgrid.jpg', dpi=300)
        
    ## plot horizontal metrics
    #plt.figure(figsize=(8,3))
    #ax=plt.axes(projection=ccrs.PlateCarree())
    #ax.set_extent(lims,ccrs.PlateCarree())
    #im = ax.pcolormesh(vlon,vlat,dxt,transform=ccrs.PlateCarree())
    #cbar = plt.colorbar(im, format="%.1f")
    #plt.title('dxt [m]', size=10) # to modify the title
    #ax.set_xticks(lon_tcks, crs=ccrs.PlateCarree())
    #ax.set_yticks(lat_tcks, crs=ccrs.PlateCarree())
    #ax.coastlines(resolution='50m') # Currently can be one of “110m”, “50m”, and “10m”
    #ax.gridlines()    
    #plt.savefig('figs/'+stest+'/nemo_input_dxt.jpg', dpi=300)

    ## plot horizontal metrics
    #plt.figure(figsize=(8,3))
    #ax=plt.axes(projection=ccrs.PlateCarree())
    #ax.set_extent(lims,ccrs.PlateCarree())
    #im = ax.pcolormesh(vlon,vlat,dyt,transform=ccrs.PlateCarree())
    #cbar = plt.colorbar(im, format="%.1f")
    #plt.title('dyt [m]', size=10) # to modify the title
    #ax.set_xticks(lon_tcks, crs=ccrs.PlateCarree())
    #ax.set_yticks(lat_tcks, crs=ccrs.PlateCarree())
    #ax.coastlines(resolution='50m') # Currently can be one of “110m”, “50m”, and “10m”
    #ax.gridlines()    
    #plt.savefig('figs/'+stest+'/nemo_input_dyt.jpg', dpi=300)

        
    # store metric terms
    #zt = np.hstack((zt,zt[[-1]]))
    print "create metrics grid"
    # metricsout = create_nc('data/nemo_metrics.nc', lon, lat, zt[index_mask_depth:], zw[index_mask_depth:])
    metricsout = create_nc('data/'+stest+'/nemo_metrics.nc', vlon, vlat, zt, zw)    #

    dtype='f8'
    #nc_zw = rootgrp.createVariable('zw',dtype,('zw'))
    #nc_zw[:] = zw
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
    # metricsout.close()
        
    ###
    
    # compute the Coriolis frequency and a reference value
    # from oocgcm/oocgcm/parameters/physicalparameters.py
    grav = oopp.grav              # acceleration due to gravity (m.s-2)
    omega = oopp.omega            # earth rotation rate (s-1)
    earthrad = oopp.earthrad      # mean earth radius (m)
    f = xp.grd._arrays['coriolis_parameter_at_t_location'].to_masked_array()
    f0 = phy.average_coriolis_parameter2(xp.region)
    
    # load stratification profile
    print "read stratification"
    N2_file = datadir+'DIAG_DIMUP/2007_2008/'+xp.region+'/bg/'+xp.region\
       +'_2007_2008_bruntvaisala_bg_mindepth10_new.nc'
    rhobg_file = datadir+'DIAG_DIMUP/2007_2008/'+xp.region+'/bg/'+xp.region\
       +'_2007_2008_density_bg_mindepth10.nc'
    rhobg_xr,N2_xr=phy.call_background(xp,xp.region,xp.grd,xp.vgrd,z=None,mindepth=10,
           N2file=N2_file,rhobgfile=rhobg_file,fN2gridT=False)
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

    # load PV
    print "compute PV"
    print "Execution might be longer than 30 minutes, please wait..."
    # compute density, psi, PV
    # fout: activate rho and psi outputting
    # fc* : activate online processing (no intermediate reading)
    arho_xr,psi_xr,q_xr= phy.call_qgpv(xp,day,mth,yr,z=None,fcdens=True,fccurl=None,fcstretch=True,
    #q_xr= phy.call_qgpv(xp,day,mth,yr,z=None,fcdens=False,fccurl=False,fcstretch=False,
        dfilt=dfilt,pfilt=pfilt,mode=2,fout=True,N2file=N2_file,rhobgfile=rhobg_file)
    print arho_xr,psi_xr,q_xr
    rho_xr = arho_xr*1000.
    #rho_xr = arho_xr*oopp.rau0
    q = q_xr.to_masked_array()

    # force FillValue
    q._FillValue=-999.
  
    # store variables

    # store PV
    print "store pv"
    pvout = create_nc('data/'+stest+'/nemo_pv.nc', vlon, vlat, zt, zw)
    #
    nc_f = pvout.createVariable('f',dtype,('y','x'))
    nc_f[:] = f    
    nc_f0 = pvout.createVariable('f0',dtype)
    nc_f0[:] = f0
    #
    nc_N2 = pvout.createVariable('N2',dtype,('zw'))
    nc_N2[:] = np.flipud(N2)
    #
    nc_q = pvout.createVariable('q',dtype,('zt','y','x'))
    nc_q[:] = np.flipud(q)

    #create 2D mask at reference level index_mask_depth (land=1, water=0)
    print "store mask"
    nc_mask = metricsout.createVariable('mask',dtype,('y','x'), fill_value=q._FillValue)
    nc_mask[:] = q[N-index_mask_depth-1,:,:]
    nc_mask[:] = np.where(nc_mask == nc_mask._FillValue, nc_mask, 1.) 
    nc_mask[:] = np.where(nc_mask != nc_mask._FillValue, nc_mask, 0.) 
    #
    # enlarge the mask: if the i,j point has an adjacent land point then it becames land
    dummy = nc_mask[1:-1,1:-1]+nc_mask[:-2,1:-1]+nc_mask[2:,1:-1]+nc_mask[1:-1,:-2]+nc_mask[1:-1,2:]
    nc_mask[1:-1,1:-1] = np.where(dummy == 5, nc_mask[1:-1,1:-1], 0.)
    metricsout.close()
    pvout.close()

    ### load and store psi    
    print "read psi"
    psi = psi_xr.to_masked_array()
    # store psi
    print "create psi file"
    # psiout = create_nc('data/nemo_psi.nc', lon, lat, zt[index_mask_depth:], zw[index_mask_depth:])
    psiout = create_nc('data/'+stest+'/nemo_psi.nc', vlon, vlat, zt, zw)
    nc_psi = psiout.createVariable('psi',dtype,('zt','y','x'))
    nc_psi[:] = np.flipud(psi) 

    ## plot psi
    #plt.figure(figsize=(8,3))
    #ax=plt.axes(projection=ccrs.PlateCarree())
    #ax.set_extent(lims,ccrs.PlateCarree())
    #im = ax.pcolormesh(vlon,vlat,psi[5,:,:],transform=ccrs.PlateCarree())
    #cbar = plt.colorbar(im, format="%.2f")
    #plt.title('psi(z=%0.0f) [m^2/s]' %zt[-5-1], size=10) # to modify the title
    #ax.set_xticks(lon_tcks, crs=ccrs.PlateCarree())
    #ax.set_yticks(lat_tcks, crs=ccrs.PlateCarree())
    #ax.coastlines(resolution='50m') # Currently can be one of “110m”, “50m”, and “10m”
    #ax.gridlines()
    #plt.savefig('figs/'+stest+'/nemo_input_psi.jpg', dpi=300)
  
    psiout.close()

    ### load rho and background  

    print "read rho"
    rho = rho_xr.to_masked_array()   
    # store rho - background
    print "create rho file"
    rhoout = create_nc('data/'+stest+'/nemo_rho.nc', vlon, vlat, zt, zw)
    #
    nc_rho = rhoout.createVariable('rho',dtype,('zt','y','x'))
    nc_rho[:] = np.flipud(rho[:,:,:])

    rhoout.close()

    # plot rho                      
    #plt.figure(figsize=(8,3))
    #ax=plt.axes(projection=ccrs.PlateCarree())
    #ax.set_extent(lims,ccrs.PlateCarree())
    #im = ax.pcolormesh(vlon,vlat,rho[5,:,:]-rhobg[5],transform=ccrs.PlateCarree())
    #cbar = plt.colorbar(im, format="%.2f")
    #plt.title('rho(z=%0.0f) [kg/m^3]' %zt[-5-1], size=10) # to modify the title
    #ax.set_xticks(lon_tcks, crs=ccrs.PlateCarree())
    #ax.set_yticks(lat_tcks, crs=ccrs.PlateCarree())
    #ax.coastlines(resolution='50m') # Currently can be one of “110m”, “50m”, and “10m”
    #ax.gridlines()
    #plt.savefig('figs/'+stest+'/nemo_input_rho.jpg', dpi=300)
    
   # plt.ion()
   # plt.show(); 
