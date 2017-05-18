#!/usr/bin/python
# -*- encoding: utf8 -*-

""" Create metrics, psi, pv fields for a curvlinear grid
"""

import os,sys
import shutil
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import netCDF4
from netCDF4 import Dataset

# lpolib
sys.path.append('/home/slyne/aponte/natl60/lporoms')
from lpolib.lporun import LPORun
from lpolib.utils import *

# maybe temporary
# import matplotlib as mpl
# from mpl_toolkits.axes_grid1 import make_axes_locatable
#from scipy.fftpack._fftpack import zfft

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

    # check number of arguments
    if  len(sys.argv) < 2:
        print '[syntaxe] : python create_input_roms.py rundir runscript'
        print 'rundir = directory created (relative to local dir)'
        print 'runscript = script that will be executed'
        quit()
   
    # get useful dirs
    startdir=os.getcwd()+'/'
    # get args
    RPATH = startdir+sys.argv[1]+'/'

    if os.path.exists(RPATH) :
        os.system('rm -Rf '+RPATH)
    os.mkdir(RPATH)
    os.chdir(RPATH)
    shutil.copytree(startdir+'../../src_parallel/',RPATH+'qgsolver')
    os.mkdir(RPATH+'input')
    os.mkdir(RPATH+'output')
    runscript = sys.argv[2]
    shutil.copy(startdir+runscript,RPATH)


    #
    # Start now to build input files
    #

    # doesn't keep levels under the depth mask_depth
    #mask_depth = -3000.
    #mask_depth = -7000.

    ### ROMS
    rpath="/home2/pharos/othr/aponte/roms_ird/"
    # rpath="/home/mulroy/slgentil/models/roms/Roms_Agrif/"
    # dirname=rpath+"curie1/jet_cfg1_wp5_2km_k1e7_TSUP5_2000a3000j"
    dirname=rpath+"caparmor/jet_cfg1_wp5_4km_k3.2e8_0a1500j"
    
    # suffix of variables used in the output
    suff='_t_dirac'
    #suff='_a'

    # time index is also added to the suffix
    #it=5
    it=11
    #it=15
    
    # Create a lporun class
    d=LPORun(dirname,verbose=True, open_nc=[], tdir_max=10)
    
    
    ### horizontal grid
    print "read horizontal grid"
    (xr,yr)=d.getXY('rho')
    (xp,yp)=d.getXY('psi')
    #print xr.shape
    
    x = xr[:-2,1:-1]
    y = yr[:-2,1:-1]
    dxt = np.diff(xr[1:-1,:-1],axis=1)
    dyt = np.diff(yr[1:,:-2],axis=0)
    #print x.shape, y.shape, dxt.shape, dyt.shape
    dxu = dxt
    dyu = dyt
    dxv = dxt
    dyv = dyt
    
    L = x.shape[1]
    M = x.shape[0]

    
    ### vertical grid, caution: level 0 in nemo correspond to the surface, level N correspond to positive depth
    print "read vertical grid"
    ssh=d.his.variables['ssh'+suff][it,:,:]
    (xr,yr)=d.getXY('rho')
    zr0=d.getZ('rho')        # Nominal z grid at t points, 1D array
    zw0=d.getZ('w')        # Nominal z grid at t points, 1D array
    zr =d.getZ('rho',ssh)    # Actual z grid at t points, 3D array
    #zw =d.getZ('w',ssh)    # Actual z grid at w points, 3D array
            
    zt = zr0
    zw = zw0
    N = zt.shape[0]
    dzt = np.diff(zw)
    dzw = np.hstack(([zt[1]-zt[0]],np.diff(zt),zt[-1]-zt[-2]))


    # find nearest index in zt for mask_depth
    #index_mask_depth = min(range(len(zt)), key=lambda i: abs(zt[i]-mask_depth))
    #print "mask reference at level:",index_mask_depth
       
    # store metric terms
    #zt = np.hstack((zt,zt[[-1]]))
    print "create metrics grid"
    # metricsout = create_nc('data/nemo_metrics.nc', lon, lat, zt[index_mask_depth:], zw[index_mask_depth:])
    metricsout = create_nc(RPATH+'input/roms_metrics.nc', x, y, zt, zw)
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
    #grav = 9.81                  # acceleration due to gravity (m.s-2)
    g=9.81
    #omega = 7.292115083046061e-5 # earth rotation rate (s-1)
    #earthrad = 6371229            # mean earth radius (m)
    #f = 2. * omega * np.sin(lat * d2r)
    #f0 = 8.5158e-5
    f = d.hgrid.f[1:-1,:-2]
    f0 = d.hgrid.f0
    
    # load density, and velocities 
    rho=d.his.variables['T'+suff][it,...]
    u=d.his.variables['u'+suff][it,...]
    v=d.his.variables['v'+suff][it,...]    
    # interpolate on the flat grid
    u=tvs_to_s(u,ssh,d)
    v=tvs_to_s(v,ssh,d)
    rho=tvs_to_s(rho,ssh,d)
    
    # compute density reference profile and take it away from the 3D density
    rho_a = rho.mean(axis=2).mean(axis=1)
    rho += -rho_a.reshape((d.N,1,1))
    
    # store stratification profile
    print "compute stratification"
    N2 = -g*np.diff(rho_a)/d.rho0/np.diff(zt)
    N2 = np.hstack((N2[0],N2,N2[-1]))
    print N2.shape
    
    # compute p from rho manually (does not agree with ROMS)
    p = np.zeros_like(rho)
    p[-1,...] = g*(d.rho0+rho_a[-1,None,None]+rho[-1,...])*ssh[None,:,:]
    for k in range(d.N-2,-1,-1):
        p[k,...] = p[k+1,...] + g*(rho[k,...]+rho[k+1,...])*0.5*(zr0[k+1,None,None]-zr0[k,None,None])


    ### store psi    
    print "create psi file"
    psiout = create_nc(RPATH+'input/roms_psi.nc', x, y, zt, zw)
    nc_psi = psiout.createVariable('psi',dtype,('zt','y','x'))
    nc_psi[:] = p[:,1:-1,:-2]/d.hgrid.f0/d.rho0
    # close file
    #psiout.close()
   
    # compute PV
    print "compute and store PV"

    pvout = create_nc(RPATH+'input/roms_pv.nc', x, y, zt, zw)
    #
    nc_f = pvout.createVariable('f',dtype,('y','x'))
    nc_f[:] = f  
    #  
    nc_f0 = pvout.createVariable('f0',dtype)
    nc_f0[:] = f0
    #
    nc_N2 = pvout.createVariable('N2',dtype,('zw'))
    nc_N2[:] = N2
    #
    nc_q = pvout.createVariable('q',dtype,('zt','y','x'))

    print zr0
    print zw0
    print N2
    #print f[:,0]
    print f0
    
    flag_stretching = False
    
    for k in np.arange(d.N):
    
        if flag_stretching :
            
            # compute relative vorticity
            lu = u[k, :, :]
            lv = v[k, :, :]
            xi = psi2rho(vorticity(lu, lv, d.hgrid))
            
            # compute vortex stretching
            if ( k==0 ):
                # use bottom density
                S = ( (rho[k+1,...]+rho[k,...])*0.5 \
                        *(zr0[k+1]-zr0[k])/(rho_a[k+1]-rho_a[k]) \
                     - rho[k,...] \
                        *(zr0[k+1]-zr0[k])/(rho_a[k+1]-rho_a[k]) \
                     ) /(zw0[k+1]-zw0[k])
            elif ( k==d.N-1 ):
                # use top density
                S = ( rho[k,...] \
                        *(zr0[k]-zr0[k-1])/(rho_a[k]-rho_a[k-1]) \
                     -(rho[k,...]+rho[k-1,...])*0.5 \
                        *(zr0[k]-zr0[k-1])/(rho_a[k]-rho_a[k-1]) \
                     ) /(zw0[k+1]-zw0[k])
            else:
                S = ( (rho[k+1,...]+rho[k,...])*0.5 \
                        *(zr0[k+1]-zr0[k])/(rho_a[k+1]-rho_a[k]) \
                     -(rho[k,...]+rho[k-1,...])*0.5 \
                        *(zr0[k]-zr0[k-1])/(rho_a[k]-rho_a[k-1]) \
                     ) /(zw0[k+1]-zw0[k])
            #S = S * d.hgrid.f
            S = S * d.hgrid.f0
    
            # assemble q
            pv= f-f0 + xi[1:-1,:-2] + S[1:-1,:-2]
            # store q
            nc_q[k,:,:] = pv
            
        else:
            
            # compute relative vorticity
            dx2 = d.hgrid.dx[:,:]*d.hgrid.dx[:,:]
            dy2 = d.hgrid.dy[:,:]*d.hgrid.dy[:,:]
    
            psi = p[k,:,:]/d.hgrid.f0/d.rho0
    
            psixx = np.zeros(dx2.shape)
            psixx[:,1:-1] = np.diff(psi[:,:], n=2, axis=1)
            psixx[:,0] = psixx[:,-2]
            psixx[:,-1]= psixx[:,1]
            psixx = np.divide(psixx,dx2)
    
            psiyy = np.zeros(dy2.shape)
            psiyy[1:-1,:] = np.diff(psi[:,:], n=2, axis=0)
            psiyy[0,:] = psiyy[1,:]
            psiyy[-1,:] = psiyy[-2,:]
            psiyy = np.divide(psiyy,dy2)
    
            xi = psixx + psiyy

            if (k == 0):
                # bottom bdy condition not used in the solver
                S = np.zeros_like(xi)[1:-1,:-2]
            elif (k == d.N - 1):
                # top bdy condition not used in the solver
                S = np.zeros_like(xi)[1:-1,:-2]
            else:
                # interior pv
                S =  ( (d.hgrid.f0**2/nc_N2[k+1])*(nc_psi[k+1,:,:]-nc_psi[k,:,:])/(zr0[k+1]-zr0[k]) - \
                       (d.hgrid.f0**2/nc_N2[k])*(nc_psi[k,:,:]-nc_psi[k-1,:,:])/(zr0[k ]-zr0[k-1]) \
                     )/(zw0[k+1]-zw0[k])
    
            # assemble q
            pv= f-f0 + xi[1:-1,:-2] + S[:]
            # store q
            nc_q[k,:,:] = pv
    
    print nc_q[:,0,0]
    psiout.close()
    
    #datadir='/home7/pharos/othr/NATL60/'
    #pv_file = datadir+'DIAG_DIMUP/qgpv/LMX/test/LMX_y2007m01d01_qgpv_v2_test_accurate.nc'
    #pvin = Dataset(pv_file, 'r')
    #q = pvin.variables['qgpv_v2']
    #print q._FillValue

    #create 2D mask at reference level index_mask_depth (land=1, water=0)
    print "store mask"
    print 
    nc_mask = metricsout.createVariable('mask',dtype,('y','x'), fill_value=-999.0)
    nc_mask[:] = 1.
    #nc_mask[:] = np.where(nc_mask == nc_mask._FillValue, nc_mask, 1.) 
    #nc_mask[:] = np.where(nc_mask != nc_mask._FillValue, nc_mask, 0.) 
    #print nc_mask[:].shape
    nc_mask[:5,:]=0.
    nc_mask[-5:,:]=0.

    # enlarge the mask: if the i,j point has an adjacent land point then it becames land
    #dummy = nc_mask[1:-1,1:-1]+nc_mask[:-2,1:-1]+nc_mask[2:,1:-1]+nc_mask[1:-1,:-2]+nc_mask[1:-1,2:]
    #nc_mask[1:-1,1:-1] = np.where(dummy == 5, nc_mask[1:-1,1:-1], 0.)
    metricsout.close()
    #pvin.close()
    pvout.close()

  
    ### store psi    
    #print "create psi file"
    #psiout = create_nc(RPATH+'input/roms_psi.nc', x, y, zt, zw)
    #nc_psi = psiout.createVariable('psi',dtype,('zt','y','x'))
    #nc_psi[:] = p[:,1:-1,:-2]/d.hgrid.f0/d.rho0
    ## close file
    #psiout.close()    


    ### store rho - rho_background   
    print "create rho file"
    rhoout = create_nc(RPATH+'input/roms_rho.nc', x, y, zt, zw)
    nc_rho = rhoout.createVariable('rho',dtype,('zt','y','x'))
    nc_rho[:] = rho[:,1:-1,:-2]    
    rhoout.close()


    # commands to execute code
    print 'You need to execute the following commands: '
    print '  bash'
    print '  source activate petsc_env'
    print '  cd '+RPATH
    print '  mpirun -np 8  python  '+runscript
    print '  (the number of mpi processes may need to adjusted here or in the run script)'

    
    # plt.ion()
    # plt.show(); 
