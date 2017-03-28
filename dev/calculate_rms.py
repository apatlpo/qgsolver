#!/usr/bin/python
# -*- encoding: utf8 -*-

import sys
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from util import outnc
from matplotlib.backends.backend_pdf import PdfPages


# ------------------------------
# Read mask
# ------------------------------
rootgrp = Dataset('data/output.nc', 'r') 
mask2d = rootgrp.variables['mask'] [...]  
rootgrp.close()

# ------------------------------
# Calculate RMS ( L*PSI - RHS)
# ------------------------------

print "------------"
print "L*PSI vs RHS"
print "------------"


# read files
rootgrp = Dataset("data/lpsiout.nc", 'r')
dummy = rootgrp.variables["Lpsi"][0,...]
Nx = dummy.shape[2]
Ny = dummy.shape[1]
Nz = dummy.shape[0]
print "Nx, Ny,Nz =",Nx,Ny,Nz
rootgrp.close()
# expand 2d mask in Z axis
mask3d = np.tile(mask2d, (Nz,1, 1))
# mask Lpsi with mask3d
Lpsi = np.ma.masked_where(mask3d==1, dummy)

rootgrp = Dataset("data/rhs.nc", 'r')
dummy = rootgrp.variables["rhs"][0,...]
rootgrp.close()
rhs = np.ma.masked_where(mask3d==1, dummy)

# Surface
tab_lpsi = np.array(Lpsi[Nz-1,1:-1,1:-1])
tab_rhs = np.array(rhs[Nz-1,1:-1,1:-1])
max_lpsi = np.abs(tab_lpsi).max()
max_rhs = np.abs(tab_rhs).max()
rms = np.sqrt(np.mean(np.square(tab_lpsi)))
rms_diff = np.sqrt(np.mean(np.square(tab_lpsi-tab_rhs)))       
print ""
print "Surface"
print "   max Lpsi = ", max_lpsi
print "   max rhs = ", max_rhs
print "   rms Lpsi = ", rms
print "   rms Lpsi-rhs =", rms_diff
print "   diff ordre de grandeur = ","%e"%(rms_diff/max_lpsi)

# Bottom
tab_lpsi = np.array(Lpsi[0,1:-1,1:-1])
tab_rhs = np.array(rhs[0,1:-1,1:-1])
max_lpsi = np.abs(tab_lpsi).max()
max_rhs = np.abs(tab_rhs).max()
rms = np.sqrt(np.mean(np.square(tab_lpsi)))
rms_diff = np.sqrt(np.mean(np.square(tab_lpsi-tab_rhs)))         
print ""
print "Bottom"
print "   max Lpsi = ", max_lpsi
print "   max rhs = ", max_rhs
print "   rms Lpsi = ", rms
print "   rms Lpsi-rhs =", rms_diff
print "   diff ordre de grandeur = ","%e"%(rms_diff/max_lpsi)

# Lateral
tab_lpsi = np.array([Lpsi[0:Nz,:,0],Lpsi[0:Nz,:,Nx-1],Lpsi[0:Nz,0,:],Lpsi[0:Nz,Ny-1,:]])
tab_rhs = np.array([rhs[0:Nz,:,0],rhs[0:Nz,:,Nx-1],rhs[0:Nz,0,:],rhs[0:Nz,Ny-1,:]])
max_lpsi = np.abs(tab_lpsi).max()
max_rhs = np.abs(tab_rhs).max()
rms = np.sqrt(np.mean(np.square(tab_lpsi)))
rms_diff = np.sqrt(np.mean(np.square(tab_lpsi-tab_rhs)))          
print ""
print "Lateral"
print "   max Lpsi = ", max_lpsi
print "   max rhs = ", max_rhs
print "   rms Lpsi = ", rms
print "   rms Lpsi-rhs =", rms_diff
print "   diff ordre de grandeur = ","%e"%(rms_diff/max_lpsi)

# Interior
tab_lpsi = np.array(Lpsi[1:Nz-1,1:-1,1:-1])
tab_rhs = np.array(rhs[1:Nz-1,1:-1,1:-1])
max_lpsi = np.abs(tab_lpsi).max()
max_rhs = np.abs(tab_rhs).max()
rms = np.sqrt(np.mean(np.square(tab_lpsi)))
rms_diff = np.sqrt(np.mean(np.square(tab_lpsi-tab_rhs)))        
print ""
print "Interior"
print "   max Lpsi = ", max_lpsi
print "   max rhs = ", max_rhs
print "   rms Lpsi = ", rms
print "   rms Lpsi-rhs =", rms_diff
print "   diff ordre de grandeur = ","%e"%(rms_diff/max_lpsi)

# ------------------------------
# Calculate RMS ( PSIout - PSIin)
# ------------------------------

print ""
print "-------------------------------------"
print "Psi from NATL60 vs Psi from inversion"
print "-------------------------------------"

rootgrp = Dataset('data/output.nc', 'r') 
z = rootgrp.variables['z'] [:] 
dummy = rootgrp.variables['psi'] [0,...]  
rootgrp.close()
psi_out = np.ma.masked_where(mask3d==1, dummy)

rootgrp = Dataset('data/input.nc', 'r')
dummy = rootgrp.variables['psi'][0,...]
rootgrp.close()
psi_in = np.ma.masked_where(mask3d==1, dummy)

# Surface
tab_in = np.array(psi_in[Nz-1,1:-1,1:-1])
tab_out = np.array(psi_out[Nz-1,1:-1,1:-1])
max_psi_in = np.abs(tab_in).max()
max_psi_out = np.abs(tab_out).max()
rms = np.sqrt(np.mean(np.square(tab_in)))
rms_diff = np.sqrt(np.mean(np.square(tab_in-tab_out)))          
print ""
print "Surface"
print "-------"
print "   max psi_in = ", max_psi_in
print "   max psi_out = ", max_psi_out
print "   rms psi_in = ", rms
print "   rms psi_in-psi_out =", rms_diff
print "   diff ordre de grandeur = ","%e"%(rms_diff/max_psi_in)

# Bottom
tab_in = np.array(psi_in[0,1:-1,1:-1])
tab_out = np.array(psi_out[0,1:-1,1:-1])
max_psi_in = np.abs(tab_in).max()
max_psi_out = np.abs(tab_out).max()
rms = np.sqrt(np.mean(np.square(tab_in)))
rms_diff = np.sqrt(np.mean(np.square(tab_in-tab_out)))          
print ""
print "Bottom"
print "-------"
print "   max psi_in = ", max_psi_in
print "   max psi_out = ", max_psi_out
print "   rms psi_in = ", rms
print "   rms psi_in-psi_out =", rms_diff
print "   diff ordre de grandeur = ","%e"%(rms_diff/max_psi_in)

# Lateral
tab_in = np.array([psi_in[0:Nz,:,0],psi_in[0:Nz,:,Nx-1],psi_in[0:Nz,0,:],psi_in[0:Nz,Ny-1,:]])
tab_out = np.array([psi_out[0:Nz,:,0],psi_out[0:Nz,:,Nx-1],psi_out[0:Nz,0,:],psi_out[0:Nz,Ny-1,:]])
max_psi_in = np.abs(tab_in).max()
max_psi_out = np.abs(tab_out).max()
rms = np.sqrt(np.mean(np.square(tab_in)))
rms_diff = np.sqrt(np.mean(np.square(tab_in-tab_out)))         
print ""
print "Lateral"
print "-------"
print "   max psi_in = ", max_psi_in
print "   max psi_out = ", max_psi_out
print "   rms psi_in = ", rms
print "   rms psi_in-psi_out =", rms_diff
print "   diff ordre de grandeur = ","%e"%(rms_diff/max_psi_in)

# Interior
tab_in = np.array(psi_in[1:Nz-1,1:-1,1:-1])
tab_out = np.array(psi_out[1:Nz-1,1:-1,1:-1])
max_psi_in = np.abs(tab_in).max()
max_psi_out = np.abs(tab_out).max()
rms = np.sqrt(np.mean(np.square(tab_in)))
rms_diff = np.sqrt(np.mean(np.square(tab_in-tab_out)))         
print ""
print "Interior"
print "-------"
print "   max psi_in = ", max_psi_in
print "   max psi_out = ", max_psi_out
print "   rms psi_in = ", rms
print "   rms psi_in-psi_out =", rms_diff
print "   diff ordre de grandeur = ","%e"%(rms_diff/max_psi_in)

# courbes en z de RMS(psi_in) et RMS(psi-in - psi_out)
rms = np.zeros_like(psi_in[0:Nz,0,0])
rms_diff = np.zeros_like(psi_in[0:Nz,0,0])
for k in range(0,Nz):
	rms[k] = np.sqrt(np.mean(np.square(psi_in[k,:,:])))
	rms_diff[k] = np.sqrt(np.mean(np.square(psi_in[k,:,:] - psi_out[k,:,:])))

# plot 

# read psi out for rtol= 1e-10,1e-7,1e-4,1e-1
rootgrp = Dataset('data/output10.nc', 'r')
psi10 = rootgrp.variables['psi'][0,...]
rootgrp.close()
rootgrp = Dataset('data/output7.nc', 'r')
psi7 = rootgrp.variables['psi'][0,...]
rootgrp.close()
rootgrp = Dataset('data/output4.nc', 'r')
psi4 = rootgrp.variables['psi'][0,...]
rootgrp.close()
rootgrp = Dataset('data/output3.nc', 'r')
psi3 = rootgrp.variables['psi'][0,...]
rootgrp.close()
rootgrp = Dataset('data/output2.nc', 'r')
psi2 = rootgrp.variables['psi'][0,...]
rootgrp.close()
rootgrp = Dataset('data/output1.nc', 'r')
psi1 = rootgrp.variables['psi'][0,...]
rootgrp.close()
# calculate for each psi
rms10 = np.zeros_like(psi_in[0:Nz,0,0])
rms7  = np.zeros_like(psi_in[0:Nz,0,0])
rms4  = np.zeros_like(psi_in[0:Nz,0,0])
rms3  = np.zeros_like(psi_in[0:Nz,0,0])
rms2  = np.zeros_like(psi_in[0:Nz,0,0])
rms1  = np.zeros_like(psi_in[0:Nz,0,0])
for k in range(0,Nz):
	rms10[k] = np.sqrt(np.mean(np.square(psi10[k,:,:])))
	rms7[k] = np.sqrt(np.mean(np.square(psi7[k,:,:])))
	rms4[k] = np.sqrt(np.mean(np.square(psi4[k,:,:])))
	rms3[k] = np.sqrt(np.mean(np.square(psi4[k,:,:])))
	rms2[k] = np.sqrt(np.mean(np.square(psi4[k,:,:])))
	rms1[k] = np.sqrt(np.mean(np.square(psi1[k,:,:])))
# create a PdfPages object
pdf = PdfPages('figs/rms.pdf')
fig1=plt.figure()
plt.plot(rms,z[0:Nz], 'k')
plt.plot(rms10,z[0:Nz], 'r', label='rtol=1.e-10')
plt.plot(rms7,z[0:Nz], 'b', label='rtol=1.e-7')
plt.plot(rms4,z[0:Nz], 'm', label='rtol=1.e-4')
# plt.plot(rms3,z[0:Nz], 'y', label='rtol=1.e-3')
# plt.plot(rms2,z[0:Nz], 'm', label='rtol=1.e-2')
plt.plot(rms1,z[0:Nz], 'c', label='rtol=1.e-1')
plt.legend(loc='best')
plt.title('RMS PSI_IN vs PSI_INV for several rtol')
pdf.savefig(fig1)

plt.savefig('figs/rms.jpg', dpi=300)
# fig2=plt.figure()
# plt.plot(rms_diff,z[0:Nz], 'k.')
# plt.title('RMS (PSI - PSI_INV)')
# pdf.savefig(fig2)
plt.ion()
plt.show()
pdf.close()


### create a netcdf file with psi_in and psi_out masked
rootgrp = Dataset('psiinmask.nc', 'w',
                  format='NETCDF4_CLASSIC', clobber=True)
# create dimensions
rootgrp.createDimension('x', Nx)
rootgrp.createDimension('y', Ny)
rootgrp.createDimension('z', Nz)
# create variables
dtype='f8'
psiin = rootgrp.createVariable('psi',dtype,('z','y','x'))