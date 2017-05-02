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

    
    # # plot hgrid
    # lims=[-80,-55, 30, 45]
    # lon_tcks = range(-80,-55, 5)
    # lat_tcks = range(30,45,5)
    # #
    # plt.figure(figsize=(8,3))
    # ax=plt.axes(projection=ccrs.PlateCarree())
    # ax.set_extent(lims,ccrs.PlateCarree())
    # ax.plot(lon[::10,::10],lat[::10,::10],'k-',transform=ccrs.PlateCarree())
    # ax.plot(lon.transpose()[::10,::10],lat.transpose()[::10,::10],'k-',transform=ccrs.PlateCarree())
    # plt.title('Horizontal grid (every 10 points)', size=10) # to modify the title
    # ax.set_xticks(lon_tcks, crs=ccrs.PlateCarree())
    # ax.set_yticks(lat_tcks, crs=ccrs.PlateCarree())
    # ax.coastlines(resolution='50m') # Currently can be one of “110m”, “50m”, and “10m”
    # ax.gridlines()
    # #figname='figs/snapshot_'+vkey.replace (" ", "_")+'_magnitude.jpg'
    # plt.savefig('figs/nemo_input_hgrid.jpg', dpi=300)
    # #print figname+' printed'
    
    
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

    # # plot vertical grid
    # plt.figure()
    # plt.plot(zw,'k+')
    # plt.plot(zw,'k.')
    # plt.grid()
    # plt.savefig('figs/nemo_input_vgrid.jpg', dpi=300)
        
    # # plot horizontal metrics
    # plt.figure(figsize=(8,3))
    # ax=plt.axes(projection=ccrs.PlateCarree())
    # ax.set_extent(lims,ccrs.PlateCarree())
    # im = ax.pcolormesh(lon,lat,e1,transform=ccrs.PlateCarree())
    # cbar = plt.colorbar(im, format="%.1f")
    # plt.title('e1 [m]', size=10) # to modify the title
    # ax.set_xticks(lon_tcks, crs=ccrs.PlateCarree())
    # ax.set_yticks(lat_tcks, crs=ccrs.PlateCarree())
    # ax.coastlines(resolution='50m') # Currently can be one of “110m”, “50m”, and “10m”
    # ax.gridlines()    
    # plt.savefig('figs/nemo_input_e1.jpg', dpi=300)

    # # plot horizontal metrics
    # plt.figure(figsize=(8,3))
    # ax=plt.axes(projection=ccrs.PlateCarree())
    # ax.set_extent(lims,ccrs.PlateCarree())
    # im = ax.pcolormesh(lon,lat,e2,transform=ccrs.PlateCarree())
    # cbar = plt.colorbar(im, format="%.1f")
    # plt.title('e2 [m]', size=10) # to modify the title
    # ax.set_xticks(lon_tcks, crs=ccrs.PlateCarree())
    # ax.set_yticks(lat_tcks, crs=ccrs.PlateCarree())
    # ax.coastlines(resolution='50m') # Currently can be one of “110m”, “50m”, and “10m”
    # ax.gridlines()    
    # plt.savefig('figs/nemo_input_e2.jpg', dpi=300)

        
    # store metric terms
    #zt = np.hstack((zt,zt[[-1]]))
    print "create metrics grid"
    # metricsout = create_nc('data/nemo_metrics.nc', lon, lat, zt[index_mask_depth:], zw[index_mask_depth:])
    metricsout = create_nc('data/nemo_metrics.nc', lon, lat, zt, zw)
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

    
    # metricsout.close()
        
    
    ###
    
    # compute the Coriolis frequency and a reference value
    # from oocgcm/oocgcm/parameters/physicalparameters.py
    grav = 9.81                  # acceleration due to gravity (m.s-2)
    omega = 7.292115083046061e-5 # earth rotation rate (s-1)
    earthrad = 6371229            # mean earth radius (m)
    f = 2. * omega * np.sin(lat * d2r)
    f0 = 8.5158e-5    
    
    # load stratification profile
    print "read stratification"
    # N2_file = datadir+'DIAG_DIMUP/2007_2008/LMX/bg/LMX_2007_2008_bruntvaisala_bg_mindepth10.nc'
    N2_file = datadir+'DIAG_DIMUP/2007_2008/LMX/bg/LMX_2007_2008_bruntvaisala_bg_mindepth10_new.nc'
    nc = Dataset(N2_file, 'r')
    N2 = nc.variables['bvfreq_bg'][:]
    nc.close()

    # load PV
    print "read PV"
    # pv_file = datadir+'DIAG_DIMUP/qgpv/LMX/test_good_new/LMX_y2007m01d01_qgpv_v1.nc'
    # pv_file = datadir+'DIAG_DIMUP/qgpv/LMX/test/LMX_y2007m01d01_qgpv_v2_test.nc'    
    pv_file = datadir+'DIAG_DIMUP/qgpv/LMX/test/LMX_y2007m01d01_qgpv_v2_test_accurate.nc'
    pvin = Dataset(pv_file, 'r')
    # q = nc.variables['qgpv_v1']
    q = pvin.variables['qgpv_v2']
    # pvin.close()


    # store PV
    print "store pv"
    pvout = create_nc('data/nemo_pv.nc', lon, lat, zt, zw)
    #
    nc_f = pvout.createVariable('f',dtype,('y','x'))
    nc_f[:] = f  
    #  
    nc_f0 = pvout.createVariable('f0',dtype)
    nc_f0[:] = f0
    #
    nc_N2 = pvout.createVariable('N2',dtype,('zw'))
    nc_N2[:] = np.flipud(N2)
    #
    # nc_q = pvout.createVariable('q',dtype,('zt','y','x'), fill_value=q._FillValue)
    nc_q = pvout.createVariable('q',dtype,('zt','y','x'))
    nc_q[:] = np.flipud(q)
  

    #create 2D mask at reference level index_mask_depth (land=1, water=0)
    print "store mask"
    nc_mask = metricsout.createVariable('mask',dtype,('y','x'), fill_value=q._FillValue)
    nc_mask[:] = q[N-index_mask_depth-1,:,:]
    nc_mask[:] = np.where(nc_mask == nc_mask._FillValue, nc_mask, 1.) 
    nc_mask[:] = np.where(nc_mask != nc_mask._FillValue, nc_mask, 0.) 

    # enlarge the mask: if the i,j point has an adjacent land point then it becames land
    dummy = nc_mask[1:-1,1:-1]+nc_mask[:-2,1:-1]+nc_mask[2:,1:-1]+nc_mask[1:-1,:-2]+nc_mask[1:-1,2:]
    nc_mask[1:-1,1:-1] = np.where(dummy == 5, nc_mask[1:-1,1:-1], 0.)
    metricsout.close()
    pvin.close()
    pvout.close()


    # # plot pv
    # plt.figure(figsize=(8,3))
    # ax=plt.axes(projection=ccrs.PlateCarree())
    # ax.set_extent(lims,ccrs.PlateCarree())
    # im = ax.pcolormesh(lon,lat,q[5,:,:]/f0,transform=ccrs.PlateCarree())
    # cbar = plt.colorbar(im, format="%.2f")
    # plt.title('q/f0(z=%0.0f) [1]' %zt[-5-1], size=10) # to modify the title
    # ax.set_xticks(lon_tcks, crs=ccrs.PlateCarree())
    # ax.set_yticks(lat_tcks, crs=ccrs.PlateCarree())
    # ax.coastlines(resolution='50m') # Currently can be one of “110m”, “50m”, and “10m”
    # ax.gridlines()
    # plt.savefig('figs/nemo_input_q.jpg', dpi=300)
    
  
    ### load and store psi    
    print "read psi"
    old=False
    if old :
        psi_file = datadir+'DIAG_DIMUP/psi0/LMX/LMX_y2007m01d01_psi0.nc'
        psiin = Dataset(psi_file, 'r')
        psi = psiin.variables['psi0'][:] 
    else:
        psi_file = datadir+'DIAG_DIMUP/psi0/LMX/LMX_y2007m01d01_psi0_split.nc'
        psiin = Dataset(psi_file, 'r')
        psi = psiin.variables['psi0_hydrostatic'][:]    
        psi_surf = psiin.variables['psi0_surface_pressure'][:]

    # store psi
    print "create psi file"
    if old:
        # psiout = create_nc('data/nemo_psi.nc', lon, lat, zt[index_mask_depth:], zw[index_mask_depth:])
        psiout = create_nc('data/nemo_psi.nc', lon, lat, zt, zw)
        nc_psi = psiout.createVariable('psi',dtype,('zt','y','x'))
        nc_psi[:] = np.flipud(psi) 
    else:
        psiout = create_nc('data/nemo_psi.nc', lon, lat, zt, zw)
        nc_psi = psiout.createVariable('psi',dtype,('zt','y','x'))
        nc_psi[:] = np.flipud(psi[:,:,:]+psi_surf[None,:,:])
  
    # # plot psi
    # plt.figure(figsize=(8,3))
    # ax=plt.axes(projection=ccrs.PlateCarree())
    # ax.set_extent(lims,ccrs.PlateCarree())
    # im = ax.pcolormesh(lon,lat,psi[5,:,:],transform=ccrs.PlateCarree())
    # cbar = plt.colorbar(im, format="%.2f")
    # plt.title('psi(z=%0.0f) [m^2/s]' %zt[-5-1], size=10) # to modify the title
    # ax.set_xticks(lon_tcks, crs=ccrs.PlateCarree())
    # ax.set_yticks(lat_tcks, crs=ccrs.PlateCarree())
    # ax.coastlines(resolution='50m') # Currently can be one of “110m”, “50m”, and “10m”
    # ax.gridlines()
    # plt.savefig('figs/nemo_input_psi.jpg', dpi=300)
    
    # close file
    psiin.close()
    psiout.close()    



    ### load rho and background
    
    print "read rho"
    rho_file = datadir+'DIAG_DIMUP/density/LMX/LMX_y2007m01d01_density.nc'
    rhoin = Dataset(rho_file, 'r')
    rho = rhoin.variables['density'][:]
    
    # load background density
    print "read rho background"
    rhobg_file = datadir+'DIAG_DIMUP/2007_2008/LMX/bg/LMX_2007_2008_density_bg_mindepth10.nc'
    rhobgin = Dataset(rhobg_file, 'r')
    rhobg = rhobgin.variables['density_bg'][:]
    
    # store rho - background
    print "create rho file"
    rhoout = create_nc('data/nemo_rho.nc', lon, lat, zt, zw)
    nc_rho = rhoout.createVariable('rho',dtype,('zt','y','x'))
    nc_rho[:] = np.flipud(rho[:,:,:] - rhobg[:,None,None])
    

    # # plot horizontal metrics
    # plt.figure(figsize=(8,3))
    # ax=plt.axes(projection=ccrs.PlateCarree())
    # ax.set_extent(lims,ccrs.PlateCarree())
    # im = ax.pcolormesh(lon,lat,rho[5,:,:]-rhobg[5],transform=ccrs.PlateCarree())
    # cbar = plt.colorbar(im, format="%.2f")
    # plt.title('rho(z=%0.0f) [kg/m^3]' %zt[-5-1], size=10) # to modify the title
    # ax.set_xticks(lon_tcks, crs=ccrs.PlateCarree())
    # ax.set_yticks(lat_tcks, crs=ccrs.PlateCarree())
    # ax.coastlines(resolution='50m') # Currently can be one of “110m”, “50m”, and “10m”
    # ax.gridlines()
    # plt.savefig('figs/nemo_input_rho.jpg', dpi=300)
    
    # close file
    rhoin.close()
    # rhobgin.close()
    rhoout.close()    


    
    # plt.ion()
    # plt.show(); 