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

from scipy.fftpack._fftpack import zfft

d2r = np.pi/180.
fillvalue = netCDF4.default_fillvals['f8']
xs=230
xe=1000
ys=160
ye=650
level = 20

def create_nc(filename, lon, lat):
    
    ### create a netcdf file
    rootgrp = Dataset(filename, 'w',
                      format='NETCDF4_CLASSIC', clobber=True)

    # create dimensions
    rootgrp.createDimension('x', lon.shape[1])
    rootgrp.createDimension('y', lat.shape[0])
    
    # create variables
    dtype='f8'
    nc_lon = rootgrp.createVariable('lon',dtype,('y','x'))
    nc_lat = rootgrp.createVariable('lat',dtype,('y','x'))
    
    nc_lon[:] = lon
    nc_lat[:] = lat
        
    # rootgrp.createVariable(name,dtype,('zc','y','x',)))
    return rootgrp



if __name__ == "__main__":
    

    ### NEMO grid file
    # datadir='/home7/pharos/othr/NATL60/'
    datadir='data/'
    griddir=datadir+'NATL60-I/BOXES/'
    
    
    ### horizontal grid
    print "read horizontal grid"
    hgrid_file=griddir+'NATL60LMX_coordinates_v4.nc'
    
    hgrid = Dataset(hgrid_file, 'r')
    lon = hgrid.variables['nav_lon'][ys:ye,xs:xe]
    lat = hgrid.variables['nav_lat'][ys:ye,xs:xe]
    e1t = hgrid.variables['e1t'][ys:ye,xs:xe]
    e2t = hgrid.variables['e2t'][ys:ye,xs:xe]
    e2u = hgrid.variables['e2u'][ys:ye,xs:xe]
    e1u = hgrid.variables['e1u'][ys:ye,xs:xe]
    e1v = hgrid.variables['e1v'][ys:ye,xs:xe]
    e2v = hgrid.variables['e2v'][ys:ye,xs:xe]
    L = lon.shape[1]
    M = lat.shape[0]

    ### vertical grid, caution: level 0 in nemo correspond to the surface, level N correspond to positive depth
    print "read vertical grid"
    vgrid_file=griddir+'NATL60LMX_v4.1_cdf_mesh_zgr.nc'
    
    vgrid = Dataset(vgrid_file, 'r')
    zc = vgrid.variables['gdept_0'][0,level-2:level+3]
    zf = vgrid.variables['gdepw_0'][0,level-2:level+3] 
    N = zc.shape[0]
    e3f = vgrid.variables['e3w_0'][0,level-2:level+3]
    e3t = vgrid.variables['e3t_0'][0,level-2:level+3]
    print zc
    print zf
        
    
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
    N2_file = datadir+'DIAG_DIMUP/2007_2008/LMX/bg/LMX_2007_2008_bruntvaisala_bg_mindepth10_new.nc'
    nc = Dataset(N2_file, 'r')
    N2 = nc.variables['bvfreq_bg'][level-2:level+3]
    nc.close()

    # load PV
    print "read PV"
    # pv_file = datadir+'DIAG_DIMUP/qgpv/LMX/test/LMX_y2007m01d01_qgpv_v2_test.nc'
    pv_file = "/home7/pharos/othr/NATL60/DIAG_DIMUP/qgpv/LMX/test/LMX_y2007m01d01_qgpv_v2_test_accurate.nc"
    pvin = Dataset(pv_file, 'r')
    q = pvin.variables['qgpv_v2'][level-2:level+3,ys:ye,xs:xe]

    # pvin.close()


    ### load and store psi
    
    print "read psi"
    old=False
    if old :
        psi_file = datadir+'DIAG_DIMUP/psi0/LMX/LMX_y2007m01d01_psi0.nc'
        psiin = Dataset(psi_file, 'r')
        psi = psiin.variables['psi0'][level-2:level+3,ys:ye,xs:xe] 
    else:
        psi_file = datadir+'DIAG_DIMUP/psi0/LMX/LMX_y2007m01d01_psi0_split.nc'
        psiin = Dataset(psi_file, 'r')
        psi = psiin.variables['psi0_hydrostatic'][level-2:level+3,ys:ye,xs:xe]    
        psi_surf = psiin.variables['psi0_surface_pressure'][ys:ye,xs:xe]
        psi[:,:,:] = psi[:,:,:] + psi_surf[None,:,:]

    print "create psi file"

    # close file
    psiin.close()

 
    # Calculate relative vorticity LAP(psi)
    # e1t is dx at t point, e2v at dy on v point
    k=1
    q_relative=np.zeros_like(psi[k,:,:])
    for j in range(1,M-1):
        for i in range(1,L-1):
            q_relative[j,i] = \
                ( ( (psi[k,j,i+1]-psi[k,j,i])*e2u[j,i]/e1u[j,i] -\
                    (psi[k,j,i]-psi[k,j,i-1])*e2u[j,i-1]/e1u[j,i-1])  +\
                  ( (psi[k,j+1,i]-psi[k,j,i])*e1v[j,i]/e2v[j,i] -\
                    (psi[k,j,i]-psi[k,j-1,i])*e1v[j-1,i]/e2v[j-1,i])\
                ) / (e1t[j,i]*e2t[j,i])


    # Calculate stretching term, caution: level 0 is the surface
    q_stretch=np.zeros_like(psi[k,:,:])
    print f.shape, lat.shape
    for j in range(1,M-1):
        for i in range(1,L-1):
            q_stretch[j,i] = ( f0**2/N2[k+1] *  (psi[k+1,j,i]-psi[k,j,i])/e3f[k+1] - \
                               f0**2/N2[k] *    (psi[k,j,i]-psi[k-1,j,i])/e3f[k]) / e3t[k]

    # store Vorticity
    fcor=np.zeros_like(psi[k,:,:])
    fcor[1:,1:] = 0.25*(f[1:,1:]+f[1:,:-1]+f[:-1,:-1]+f[:-1,1:])
    q_sum = q_relative + q_stretch + f - f0

    # store original PV 
    dtype='f8'
    pv1 = create_nc('data/pv_orig.nc', lon, lat)
    nc_q1 = pv1.createVariable('q_orig',dtype,('y','x'))
    nc_q1[:] = q[k,:,:]

    # store calculated PV 
    pv2 = create_nc('data/pv_calc.nc', lon, lat)
    nc_q2 = pv2.createVariable('q_calc',dtype,('y','x'))
    nc_q2[:] = q_sum

    # store difference PV 
    pv3 = create_nc('data/pv_diff.nc', lon, lat)
    nc_q3 = pv3.createVariable('q_diff',dtype,('y','x'))
    nc_q3[:] = nc_q1[:] - nc_q2[:]
