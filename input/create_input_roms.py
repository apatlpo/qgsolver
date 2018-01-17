#!/usr/bin/python
# -*- encoding: utf8 -*-

""" Create metrics, psi, pv fields from a ROMS output in an idealized setup 
"""

import os,sys
import argparse
import shutil
import numpy as np
#import netCDF4
from netCDF4 import Dataset, default_fillvals

from utils import *

# lpolib
sys.path.append('/home/slyne/aponte/lporoms')
from lpolib.lporun import LPORun
from lpolib.utils import *


d2r = np.pi/180.
dtype='f8'
fillvalue = default_fillvals[dtype]


def metrics():
    ''' Compute and store metric terms
    '''
    
    ### horizontal grid
    print('read horizontal grid')
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

    ### vertical grid
    print('read vertical grid')
    ssh=d.his.variables['ssh'+suff][it,:,:]
    (xr,yr)=d.getXY('rho')
    zr0=d.getZ('rho')        # Nominal z grid at t points, 1D array
    zw0=d.getZ('w')        # Nominal z grid at t points, 1D array
    zr =d.getZ('rho',ssh)    # Actual z grid at t points, 3D array
    #zw =d.getZ('w',ssh)    # Actual z grid at w points, 3D array
            
    zt = zr0
    # skip the first w point (the bottom) to keep the same convention as in nemo. T[0] is below w[0]
    zw = zw0[1:]
    N = zt.shape[0]
    # dzt = np.diff(zw)
    dzt = np.hstack((zw[0]-zt[0],np.diff(zw)))
    # dzw = np.hstack(([zt[1]-zt[0]],np.diff(zt),zt[-1]-zt[-2]))
    dzw = np.hstack((np.diff(zt),zw[-1]-zt[-1]))
       
    # store metric terms
    #zt = np.hstack((zt,zt[[-1]]))
    print('create metrics grid')
    metricsout = create_nc(outdir+'/roms_metrics.nc', x, y, zt, zw)
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

    return metricsout, ssh, x, y, zt, zw


def store_psi(background=False):
    if not background:
        print('create psi file')
        psiout = create_nc(outdir+'roms_psi.nc', x, y, zt, zw)
        nc_psi = psiout.createVariable('psi',dtype,('zt','y','x'))
        nc_psi[:] = p[:,1:-1,:-2]/d.hgrid.f0/d.rho0
        return psiout    
    else:
        print('create a file for background fields')
        bgout = create_nc(outdir+'roms_bg.nc', x, y, zt, zw)
        nc_psi = bgout.createVariable('psi',dtype,('zt','y','x'))
        nc_psi[:] = np.tile(p[:,1:-1,:-2].mean(axis=2,keepdims=True),(1,1,nc_psi[:].shape[2])) \
                    /d.hgrid.f0/d.rho0
        return bgout
        
def store_rho(background=False):
    if not background:
        print('create rho file')
        rhoout = create_nc(outdir+'roms_rho.nc', x, y, zt, zw)
        nc_rho = rhoout.createVariable('rho',dtype,('zt','y','x'))
        nc_rho[:] = rho[:,1:-1,:-2]
        return rhoout
    else:
        print('store background density')
        nc_rho = bgout.createVariable('rho',dtype,('zt','y','x'))
        nc_rho[:] = np.tile(rho[:,1:-1,:-2].mean(axis=2,keepdims=True),(1,1,nc_rho[:].shape[2]))

def store_pv(background=False):
    if not background:
        print('compute and store PV')
        
        # compute stratification profile
        N2 = -g*np.diff(rho_a)/d.rho0/np.diff(zt)
        # N2 = np.hstack((N2[0],N2,N2[-1]))
        N2 = np.hstack((N2,N2[-1]))
        
        pvout = create_nc(outdir+'roms_pv.nc', x, y, zt, zw)
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
        lnc_q = pvout.createVariable('q',dtype,('zt','y','x'))
    
        print('zt from ',zt[0],' to ',zt[-1])
        print('zw from ',zw[0],' to ',zw[-1])
        
        flag_pv = 0    
        for k in np.arange(d.N):
        
            if flag_pv == 0 :
                
                # compute relative vorticity
                lu = u[k, :, :]
                lv = v[k, :, :]
                xi = psi2rho(vorticity(lu, lv, d.hgrid))
                
                # compute vortex stretching
                if ( k==0 ):
                    # use bottom density                
                    # bottom bdy condition not used in the solver
                    S = np.zeros_like(rho[0,:,:])
                elif ( k==d.N-1 ):
                    # use top density              
                    # bottom bdy condition not used in the solver
                    S = np.zeros_like(rho[0,:,:])
                else:
                    S = ( (rho[k+1,...]+rho[k,...])*0.5 \
                            *(zt[k+1]-zt[k])/(rho_a[k+1]-rho_a[k]) \
                         -(rho[k,...]+rho[k-1,...])*0.5 \
                            *(zt[k]-zt[k-1])/(rho_a[k]-rho_a[k-1]) \
                         ) /(zw[k]-zw[k-1])
                #S = S * d.hgrid.f
                S = S * d.hgrid.f0
        
                # assemble q
                pv= f-f0 + xi[1:-1,:-2] + S[1:-1,:-2]
                # store q
                lnc_q[k,:,:] = pv
                
            elif flag_pv == 1:
                
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
                    S =  ( (d.hgrid.f0**2/nc_N2[k])*(nc_psi[k+1,:,:]-nc_psi[k,:,:])/(zt[k+1]-zt[k]) - \
                           (d.hgrid.f0**2/nc_N2[k-1])*(nc_psi[k,:,:]-nc_psi[k-1,:,:])/(zt[k ]-zt[k-1]) \
                         )/(zw[k]-zw[k-1])
        
                # assemble q
                pv= f-f0 + xi[1:-1,:-2] + S[:]
    
                # store q
                lnc_q[k,:,:] = pv
                
            elif flag_pv ==2:
                
                # compute relative vorticity
                lu = u[k, :, :]
                lv = v[k, :, :]
                xi = psi2rho(vorticity(lu, lv, d.hgrid))
                
                # compute vortex stretching
                if ( k==0 ):
                    # use bottom density
                    S = np.zeros_like(xi)
                    # S = 0.
                elif ( k==d.N-1 ):
                    # use top density
                    S = np.zeros_like(xi)
                    # S = 0.
                else:
                    S = d.hgrid.f0 * ( -g*(rho[k+1,...]+rho[k,...])*0.5 /nc_N2[k] \
                                       +g*(rho[k,...]+rho[k-1,...])*0.5 /nc_N2[k-1] \
                                       ) /(zw[k]-zw[k-1])
        
                # assemble q
                pv= f-f0 + xi[1:-1,:-2] + S[1:-1,:-2]
                # store q
                lnc_q[k,:,:] = pv
        return pvout, lnc_q
    else:
        print('compute and store background PV')
        #
        nc_qbg = bgout.createVariable('q',dtype,('zt','y','x'))
        nc_qbg[:] = np.tile(nc_q[:].mean(axis=2,keepdims=True),(1,1,nc_q[:].shape[2]))
        
   


if __name__ == "__main__":
    ''' Main code to read a ROMS output and create a qgsolver input
    
    Usage
    -----
    python create_input_roms.py 
    python create_input_roms.py input_dir/  # place files in input_dir/
     
    '''

    # default directories
    romsdir = '/home2/pharos/othr/aponte/roms_ird/caparmor/jet_cfg1_wp5_4km_k3.2e8_0a1500j/'
    outdir = './input/'    

    # parse input arguments
    parser = argparse.ArgumentParser(description='Create input files from ROMS output')
    parser.add_argument('-r', dest='romsdir', help='directory where ROMS outputs are found',
                        default=romsdir)
    parser.add_argument('-o', dest='outdir', help='destination directory',
                        default=outdir)    
    args = parser.parse_args()
    #
    romsdir = args.romsdir    
    #
    outdir = args.outdir   
    if os.path.exists(outdir):
        if query_yes_no(outdir+' exists, erase and replace? ', default="yes"):
            os.system('rm -Rf '+outdir)
            os.mkdir(outdir)
        else:
            print('Choose a different name then')
            sys.exit()
    else:
        os.mkdir(outdir)
    #
    print('Input directory: '+outdir)
    print('ROMS output directory: \n\t'+romsdir)
            

    #
    # Start now to build input files
    #
    
    # suffix of variables used in the output
    suff='_t_dirac'
    #suff='_a'

    # time index is also added to the suffix
    #it=5
    it=11
    #it=15
    
    # Create a lporun class
    d=LPORun(romsdir,verbose=True, open_nc=[], tdir_max=10)
    
    # metric terms
    metricsout, ssh, x, y, zt, zw = metrics()
            
    # compute the Coriolis frequency and a reference value
    g=9.81
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
    
    # compute p from rho manually (does not agree with ROMS)
    p = np.zeros_like(rho)
    p[-1,...] = g*(d.rho0+rho_a[-1,None,None]+rho[-1,...])*ssh[None,:,:]
    for k in range(d.N-2,-1,-1):
        p[k,...] = p[k+1,...] + g*(rho[k,...]+rho[k+1,...])*0.5*(zt[k+1,None,None]-zt[k,None,None])

    # store psi
    psiout = store_psi()
    bgout = store_psi(background=True)

    # compute PV
    pvout, nc_q = store_pv()
    store_pv(background=True)
    
    #create 2D mask at reference level index_mask_depth (land=1, water=0)
    print('store mask')
    nc_mask = metricsout.createVariable('mask',dtype,('y','x'), fill_value=-999.0)
    nc_mask[:] = 1.
    nc_mask[:5,:]=0.
    nc_mask[-5:,:]=0.
    
    # close files
    metricsout.close()
    psiout.close()    
    pvout.close()
    
    # store rho - rho_background  
    rhoout = store_rho()
    store_rho(background=True)
    rhoout.close()

    # commands to execute code
    print('You need to execute the following commands: ')
    print('  bash')
    print('  source activate petsc')
    print('  cd '+outdir)
    print('  mpirun -np 8  python run.py')
    print('  (the number of mpi processes may need to adjusted here or in the run script)')

