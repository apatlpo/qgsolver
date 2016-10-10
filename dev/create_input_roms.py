#
# Compute quasigeostrophic PV and store it for latter inversion
# and time stepping
# ROMS variable are interpolated on a flat grid.
#
# requires lporoms library: Vertical grid are (flat) zlevels:
#   https://apatlpo@bitbucket.org/mdunphy/lporoms.git

#
import sys;
from lpolib.lporun import LPORun
#from lpolib.vmodes import VModes
from lpolib.utils import *
#from IPython import embed
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
# for interpolation
#from scipy import interpolate
# for colorbar positioning
from mpl_toolkits.axes_grid1 import make_axes_locatable
# netcdf
from netCDF4 import Dataset
# solves the sparse linear pb
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve


#change default fontsize
mpl.rcParams['font.size']=10
mpl.use('Agg')

# flag for printing
prtfig=True

# main path to data
rpath="/home2/pharos/othr/aponte/roms_ird/"
dirname=rpath+"curie1/jet_cfg1_wp5_2km_k1e7_TSUP5_2000a3000j"

# prefix use for figure outputs
pref=dirname.replace(rpath,"")
pref=pref.replace("caparmor/","")
pref=pref.replace("fermi/","")
pref=pref.replace("curie1/","")
pref=pref+"_zlvl_"

# suffix of variables used in the output
suff='_t_dirac'
#suff='_a'

# time index is also added to the suffix
#it=5
it=11
#it=15

# output file
pvfile='data/'+pref+'pv.nc'
#pvfile='data/'+pref+'pv'+suff+'_t'+str(it)+'.nc'
print "Output file is: "+pvfile

# Create a lporun class
d=LPORun(dirname,verbose=True, open_nc=[], tdir_max=10)

# time variable
#it=10
ti=d.his.variables['time_instant'][it]/86400
print 'Input time selected = %.0f' %(ti)

# get ssh and compute the vertical grid
ssh=d.his.variables['ssh'+suff][it,:,:]
(xr,yr)=d.getXY('rho')
zr0=d.getZ('rho')        # Nominal z grid at t points, 1D array
zw0=d.getZ('w')        # Nominal z grid at t points, 1D array
zr =d.getZ('rho',ssh)    # Actual z grid at t points, 3D array
zw =d.getZ('w',ssh)    # Actual z grid at w points, 3D array

# target z levels
zr1=zr0.reshape((d.N,1,1))+np.zeros_like(zr)
zw1=zw0.reshape((d.N+1,1,1))+np.zeros_like(zw)

# print nominal grid info
print "Target z levels zr1 at the surface = " + '%f' %zr1[-1,0,0]
print "Target z levels zw1 at the surface = " + '%f' %zw1[-1,0,0]

# plot ssh from zw variable
plt.figure(1)
plt.ion()
plt.show()
plt.contourf(xr/1e3,yr/1e3,zw[-1,:,:],20)
plt.axes().set_aspect('equal', 'box','C')
plt.colorbar();
plt.title('sea level [m], ' + str(ti) +"d")


# load density, and velocities 
rho=d.his.variables['T'+suff][it,...]
u=d.his.variables['u'+suff][it,...]
v=d.his.variables['v'+suff][it,...]

# interpolate on the flat grid
#rho=tvs_to_s(rho,ssh,d,ztarget=zr0)
#u=tvs_to_s(u,ssh,d,ztarget=zr0)
#v=tvs_to_s(v,ssh,d,ztarget=zr0)
#
rho=tvs_to_s_fast(rho,ssh,d)
u=tvs_to_s_fast(u,ssh,d)
v=tvs_to_s_fast(v,ssh,d)


# compute density reference profile and take it away from the 3D density
# average horizontally
rho_a = rho.mean(axis=2).mean(axis=1)
rho += -rho_a.reshape((d.N,1,1))


# plot the background stratification profile
plt.figure(2)
plt.ion()
plt.show()
plt.plot(rho_a,zr0,'k', lw=2, label='Mean profile')
plt.grid()
if prtfig:
    plt.savefig("figs/"+pref+"rho_bg.pdf",)


# plot a transect of density anomaly with respect to the reference
i=d.Lm/2
yrz=np.tile(yr[:,i],[d.N,1])
lvls=20

# plot total density
f, ax = plt.subplots(2, sharex=True, figsize=(6,5))
im0=ax[0].contourf(yrz/1e3,zr1[:,:,i],rho[:,:,i]+rho_a.reshape((d.N,1)) \
                    ,lvls,cmap=plt.cm.RdYlBu)
divider0 = make_axes_locatable(ax[0])
cax0 = divider0.append_axes("right", size="5%", pad=0.05)
cbar0 = plt.colorbar(im0, cax=cax0, format="%.2f")
cbar0.set_label('[kg/m^3]')
ax[0].set_title('Total density')
ax[0].grid(True)
ax[0].set_ylabel('z [m]')

# plot anomalies
im1=ax[1].contourf(yrz/1e3,zr[:,:,i],rho[:,:,i],lvls,cmap=plt.cm.RdYlBu)
divider1 = make_axes_locatable(ax[1])
cax1 = divider1.append_axes("right", size="5%", pad=0.05)
cbar1 = plt.colorbar(im1, cax=cax1, format="%.2f")
cbar1.set_label('[kg/m^3]')
ax[1].set_title('Density anomaly')
ax[1].grid(True)
ax[1].set_xlabel('y [km]')
ax[1].set_ylabel('z [m]')

if prtfig:
    plt.savefig("figs/"+pref+"rho_ano_yz"+suff+".pdf")



### create a netcdf file to store QG pv for inversion
rootgrp = Dataset(pvfile, 'w', format='NETCDF4_CLASSIC', clobber=True)

# create dimensions
rootgrp.createDimension('zc', d.N)
rootgrp.createDimension('zf', d.N+1)
rootgrp.createDimension('x', d.L-1)
rootgrp.createDimension('y', d.M+1)

# create variables
dtype='f8'
nc_zc = rootgrp.createVariable('zc',dtype,('zc',))
nc_zf = rootgrp.createVariable('zf',dtype,('zf',))
nc_x = rootgrp.createVariable('x',dtype,('y','x'))
nc_y = rootgrp.createVariable('y',dtype,('y','x'))
nc_f0 = rootgrp.createVariable('f0',dtype)
#nc_rho0 = rootgrp.createVariable('rho0',dtype)
#nc_it = rootgrp.createVariable('it',dtype)
# pv inv variables
nc_pv = rootgrp.createVariable('pv',dtype,('zc','y','x',))
nc_rho = rootgrp.createVariable('rho',dtype,('zc','y','x',))
nc_N2 = rootgrp.createVariable('N2',dtype,('zf',))
nc_f = rootgrp.createVariable('f',dtype,('y','x',))

# streamfunction, i.e. the solution presumably
nc_psi = rootgrp.createVariable('psi',dtype,('zc','y','x',))


# fills in coordinate variables, keep truely periodic data
nc_zc[:]=zr0[:]
nc_zf[:]=zw1[:]
nc_x[:]=d.hgrid.x_rho[:,:-2]
nc_y[:]=d.hgrid.y_rho[:,:-2]
nc_f[:]=d.hgrid.f[:]
#nc_f[:]=d.hgrid.f0
nc_f0[:]=d.hgrid.f0
#nc_rho0[:]=d.rho0
#nc_it[:]=it


# useful parameters
g=9.81

# start computation of pv and N2
for k in np.arange(d.N):

    # compute relative vorticity
    lu=u[k,:,:]
    lv=v[k,:,:]
    xi=psi2rho(vorticity(lu,lv,d.hgrid))
    # compute vortex stretching
    if ( k==0 ):
        # use bottom density
        S = ( (rho[k+1,...]+rho[k,...])*0.5
                *(zr0[k+1]-zr0[k])/(rho_a[k+1]-rho_a[k])
             - rho[k,...]
                *(zr0[k+1]-zr0[k])/(rho_a[k+1]-rho_a[k])
             ) /(zw0[k+1]-zw0[k]) ;
    elif ( k==d.N-1 ):
        # use top density
        S = ( rho[k,...]
                *(zr0[k]-zr0[k-1])/(rho_a[k]-rho_a[k-1])
             -(rho[k,...]+rho[k-1,...])*0.5
                *(zr0[k]-zr0[k-1])/(rho_a[k]-rho_a[k-1])
             ) /(zw0[k+1]-zw0[k]) ;
    else:
        S = ( (rho[k+1,...]+rho[k,...])*0.5
                *(zr0[k+1]-zr0[k])/(rho_a[k+1]-rho_a[k])
             -(rho[k,...]+rho[k-1,...])*0.5
                *(zr0[k]-zr0[k-1])/(rho_a[k]-rho_a[k-1])
             ) /(zw0[k+1]-zw0[k]) ;
    #S = S * d.hgrid.f
    S = S * d.hgrid.f0
    # assemble pv
    pv=xi+S
    # store pv 
    nc_pv[k,:,:]=pv[:,:-2]
    # store N2
    if k>0:
        nc_N2[k] = -g*(rho_a[k]-rho_a[k-1])/(zr0[k]-zr0[k-1])
    # sync data to netcdf file
    rootgrp.sync()
    print k, "depth level done"


# extrapolate N2 top and bottom values
nc_N2[0]=nc_N2[1]
nc_N2[-1]=nc_N2[-2]

### fill's in density
nc_rho[:] = rho[:,:,:-2]

# add a good estimate of the streamfunction from pressure
p=get_p(rho, zr, zw, d)
p=tvs_to_s(p,ssh,d,ztarget=zr0)
nc_psi[:]=p[:,:,:-2]/d.hgrid.f0/d.rho0
#nc_psi[:]=p[:,:,:-2]/d.hgrid.f/d.rho0

# close the netcdf file
rootgrp.close()



