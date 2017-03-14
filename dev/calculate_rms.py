#!/usr/bin/python
# -*- encoding: utf8 -*-

import sys
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from util import outnc
from matplotlib.backends.backend_pdf import PdfPages

# ------------------------------
# Calculate RMS ( L*PSI - RHS)
# ------------------------------

print "------------"
print "L*PSI vs RHS"
print "------------"


# read files
rootgrp = Dataset("data/lpsi.nc", 'r')
Lpsi = rootgrp.variables["Lpsi"][0,...]
Nx = Lpsi.shape[2]
Ny = Lpsi.shape[1]
Nz = Lpsi.shape[0]
print "Nx, Ny,Nz =",Nx,Ny,Nz
rootgrp.close()

rootgrp = Dataset("data/rhs.nc", 'r')
rhs = rootgrp.variables["rhs"][0,...]
rootgrp.close()

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
psi_out = rootgrp.variables['psi'] [0,...]  
rootgrp.close()

rootgrp = Dataset('data/input.nc', 'r')
psi_in = rootgrp.variables['psi'][0,...]
rootgrp.close()


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
# create a PdfPages object
pdf = PdfPages('figs/rms.pdf')
fig1=plt.figure()
plt.plot(rms,z[0:Nz], 'k--')
plt.title('RMS PSI from NATL60')
pdf.savefig(fig1)
# plt.savefig('figs/rms.jpg', dpi=300)
fig2=plt.figure()
plt.plot(rms_diff,z[0:Nz], 'k.')
plt.title('RMS (PSI - PSI_INV)')
pdf.savefig(fig2)
plt.ion()
plt.show()
pdf.close()