#!/usr/bin/python
# -*- encoding: utf8 -*-

""" Read result from a PV inversion
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


if __name__ == "__main__":
    
    # data path
    # datadir='/home/caparmor-work/aponte/qg_nemo/dev/data/'
    datadir='/home/caparmor-work/slgentil/nemotest/dev/data/'
    # datadir='/home/mulroy/slgentil/models/natl60/qgsolver/dev/data/'


    # NEMO output file
    output_file=datadir+'output.nc'
    nco = Dataset(output_file, 'r')
    q_out = nco.variables['q']   
    psi_out = nco.variables['psi']   

    # NEMO input psi and q
    input_file=datadir+'input.nc'
    nci = Dataset(input_file, 'r')
    q_in = nci.variables['q']
    psi_in = nci.variables['psi']
    
    # grid info
    grid_file=datadir+'nemo_metrics.nc'
    ncg = Dataset(grid_file, 'r')
    lon = ncg.variables['lon'][:,:]
    lat = ncg.variables['lat'][:,:]
    hdom = {'istart': 0, 'iend': 448-1, 'i0': 230,'jstart': 0, 'jend': 256-1,  'j0': 200}
    # hdom = {'istart': 0, 'iend': 100-1, 'i0': 100,'jstart': 0, 'jend': 100-1,  'j0': 400}
    
    lon = lon[hdom['j0']:hdom['j0']+hdom['jend']+1,hdom['i0']:hdom['i0']+hdom['iend']+1]
    lat = lat[hdom['j0']:hdom['j0']+hdom['jend']+1,hdom['i0']:hdom['i0']+hdom['iend']+1]
    #lon=nci.variables['x']
    #lat=nci.variables['y']

    z=nci.variables['z']
  
    # plt parameters

    xmin = np.floor(lon[0, 0])
    xmax = np.floor(lon[-1, -1]) + 1
    ymin = np.floor(lat[0, 0])
    ymax = np.floor(lat[-1, -1]) + 1
    lims=[xmin,xmax,ymin,ymax]
    lon_tcks = range(np.int(xmin),np.int(xmax),1)
    lat_tcks = range(np.int(ymin),np.int(ymax),1)

    # depth selection
    zlvl = [0, -200, -1000]
    print z[:]
    klvl = [np.argmin(z[:]<lz) for lz in zlvl]
    print klvl
  
    for k in klvl:
    #k=klvl[0]
    
        # plot psi input
        plt.figure(figsize=(8,3))
        ax=plt.axes(projection=ccrs.PlateCarree())
        ax.set_extent(lims,ccrs.PlateCarree())
        im = ax.pcolormesh(lon,lat,psi_in[0,k,:,:],transform=ccrs.PlateCarree())
        cbar = plt.colorbar(im, format="%.2f")
        plt.title('psi(z=%0.0f) in [m^2/s]' %z[k], size=10) # to modify the title
        ax.set_xticks(lon_tcks, crs=ccrs.PlateCarree())
        ax.set_yticks(lat_tcks, crs=ccrs.PlateCarree())
        ax.coastlines(resolution='50m') # Currently can be one of “110m”, “50m”, and “10m”
        ax.gridlines()
        plt.savefig('figs/nemo_input_v2_psi_in_k'+str(k)+'.jpg', dpi=300)
                
    
        # plot psi out
        plt.figure(figsize=(8,3))
        ax=plt.axes(projection=ccrs.PlateCarree())
        ax.set_extent(lims,ccrs.PlateCarree())
        im = ax.pcolormesh(lon,lat,psi_out[0,k,:,:],transform=ccrs.PlateCarree())
        cbar = plt.colorbar(im, format="%.2f")
        plt.title('psi(z=%0.0f) out [m^2/s]' %z[k], size=10) # to modify the title
        ax.set_xticks(lon_tcks, crs=ccrs.PlateCarree())
        ax.set_yticks(lat_tcks, crs=ccrs.PlateCarree())
        ax.coastlines(resolution='50m') # Currently can be one of “110m”, “50m”, and “10m”
        ax.gridlines()
        plt.savefig('figs/nemo_input_v2_psi_out_k'+str(k)+'.jpg', dpi=300)
    
    
    
        # plot psi and q input
        plt.figure(figsize=(8,8))
        #f, ax = plt.subplots(2,1, figsize=(8,8))

        #lax = ax[0]
        #lax=plt.axes(projection=ccrs.PlateCarree())
        lax = plt.subplot(2, 1, 1, projection=ccrs.PlateCarree())
        lax.set_extent(lims,ccrs.PlateCarree())
        im = lax.pcolormesh(lon,lat,q_in[0,k,:,:],transform=ccrs.PlateCarree())
        im_ctr = lax.contour(lon,lat,psi_in[0,k,:,:],20,colors='k',linewidths=0.5, transform=ccrs.PlateCarree())
        cbar = plt.colorbar(im, format="%.2e")
        lax.set_title('q psi(z=%0.0f) in [m^2/s]' %z[k], size=10) # to modify the title
        lax.set_xticks(lon_tcks, crs=ccrs.PlateCarree())
        lax.set_yticks(lat_tcks, crs=ccrs.PlateCarree())
        lax.coastlines(resolution='50m') # Currently can be one of “110m”, “50m”, and “10m”
        lax.gridlines()
        #plt.savefig('figs/nemo_input_v2_qpsi_in_k'+str(k)+'.jpg', dpi=300)

        #lax = ax[1]
        #lax=plt.axes(projection=ccrs.PlateCarree())
        lax = plt.subplot(2, 1, 2, projection=ccrs.PlateCarree())
        lax.set_extent(lims,ccrs.PlateCarree())
        im = lax.pcolormesh(lon,lat,q_out[0,k,:,:],transform=ccrs.PlateCarree())
        im_ctr = lax.contour(lon,lat,psi_out[0,k,:,:],20,colors='k',linewidths=0.5, transform=ccrs.PlateCarree())
        cbar = plt.colorbar(im, format="%.2e")
        lax.set_title('q psi(z=%0.0f) out [m^2/s]' %z[k], size=10) # to modify the title
        lax.set_xticks(lon_tcks, crs=ccrs.PlateCarree())
        lax.set_yticks(lat_tcks, crs=ccrs.PlateCarree())
        lax.coastlines(resolution='50m') # Currently can be one of “110m”, “50m”, and “10m”
        lax.gridlines()
        
        plt.savefig('figs/nemo_input_v2_qpsi_inout_k'+str(k)+'.jpg', dpi=300)

    sys.exit()
  
  
  
  
  
  
    # NEMO output file
    datadir='data/home/caparmor-work/aponte/qg_nemo/dev/data/output.nc'
    griddir=datadir+'NATL60-I/BOXES/'
    
    
    ### horizontal grid
    hgrid_file=griddir+'NATL60LMX_coordinates_v4.nc'
    
    hgrid = Dataset(hgrid_file, 'r')
    lon = hgrid.variables['nav_lon'][:]
    lat = hgrid.variables['nav_lat'][:]
    e1 = hgrid.variables['e1t'][:]
    e2 = hgrid.variables['e2t'][:]
    
    
    # plot hgrid
    # lims=[-80,-55, 30, 45]
    # lon_tcks = range(-80,-55, 5)
    # lat_tcks = range(30,45,5)
    lims=[-71,-65, 31, 38]
    lon_tcks = range(-71,-65, 5)
    lat_tcks = range(31,38,5)
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
    plt.savefig('figs/nemo_input_v2_hgrid.jpg', dpi=300)
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
    plt.savefig('figs/nemo_input_v2_vgrid.jpg', dpi=300)
        
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
    plt.savefig('figs/nemo_input_v2_e1.jpg', dpi=300)

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
    plt.savefig('figs/nemo_input_v2_e2.jpg', dpi=300)

        
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
    plt.savefig('figs/nemo_input_v2_q.jpg', dpi=300)

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
    plt.savefig('figs/nemo_input_v2_psi.jpg', dpi=300)
    
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
    plt.savefig('figs/nemo_input_v2_rho.jpg', dpi=300)
    
    # close file
    nc.close()    


    
    plt.ion()
    plt.show(); 