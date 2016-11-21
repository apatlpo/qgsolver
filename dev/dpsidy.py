#!/usr/bin/python
# -*- encoding: utf8 -*-

import sys
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from util import outnc

# choose variable
# variable = "psi"
variable = "u"
# variable = "q"
# variable = "bdy"


dy=4000.

if variable=="psi":
    # Level to plot
    level=-1
    # open netdc file
    # rootgrp1 = Dataset("data/output_run1.nc", 'r')
    # rootgrp2 = Dataset("data/output_run2.nc", 'r')
    # rootgrp1 = Dataset("data/output_run3.nc", 'r')
    # rootgrp2 = Dataset("data/output_run4.nc", 'r')
    # rootgrp1 = Dataset("data/output_run1bis.nc", 'r')
    # rootgrp2 = Dataset("data/output_run3bis.nc", 'r')
    rootgrp1 = Dataset("data/output_run2bis.nc", 'r')
    rootgrp2 = Dataset("data/output_run4bis.nc", 'r')
    # titles
    # title1 = "run1"
    # title2 = "run2"
    # title1 = "run3"
    # title2 = "run4"
    # title1 = "run1bis"
    # title2 = "run3bis"
    title1 = "run2bis"
    title2 = "run4bis"
    # read variable
    var1 = rootgrp1.variables["psi"][-1,level,:,:]
    var2 = rootgrp2.variables["psi"][-1,level,:,:]
elif variable=="u":
    # Level to plot
    level=-1
    # open netdc file
    # rootgrp1 = Dataset("data/output_run1.nc", 'r')
    # rootgrp1 = Dataset("data/output_run2.nc", 'r')
    # rootgrp1 = Dataset("data/output_run3.nc", 'r')
    # rootgrp2 = Dataset("data/output_run4.nc", 'r')
    rootgrp1 = Dataset("data/output_run5.nc", 'r')
    rootgrp2 = Dataset("data/output_run6.nc", 'r')
    # rootgrp1 = Dataset("data/output_run1bis.nc", 'r')
    # rootgrp2 = Dataset("data/output_run3bis.nc", 'r')
    # rootgrp1 = Dataset("data/output_run2bis.nc", 'r')
    # rootgrp2 = Dataset("data/output_run4bis.nc", 'r')
    # titles
    # title1 = "run1"
    # title1 = "run2"
    # title1 = "run3"
    # title2 = "run4"

    title1 = "run5"
    title2 = "run6"
    #  title1 = "run1bis"
    # title2 = "run3bis"
    # title1 = "run2bis"
    # title2 = "run4bis"
    # read variable
    psi1 = rootgrp1.variables["psi"][-1,level,:,:]
    psi2 = rootgrp2.variables["psi"][-1,level,:,:]
    var1 = np.diff(psi1,axis=0)/dy
    var2 = np.diff(psi2,axis=0)/dy
elif variable=="q":
    # Level to plot
    level=-2
    # open netdc file
    rootgrp1 = Dataset("data/jet_cfg1_wp5_4km_k3.2e8_0a1500j_zlvl_pv_from_psi.nc", 'r')
    rootgrp2 = Dataset("data/jet_cfg1_wp5_4km_k3.2e8_0a1500j_zlvl_pv_from_rho.nc", 'r')
    # titles
    title1 = "qpsi"
    title2 = "qrho"
    var1 = rootgrp1.variables["q"][level,:,:]
    var2 = rootgrp2.variables["q"][level,:,:]
elif variable=="bdy":
    # open netdc file
    rootgrp1 = Dataset("data/bdypsi.nc", 'r')
    rootgrp2 = Dataset("data/bdyrho.nc", 'r')
    # titles
    title1 = "bdypsi"
    title2 = "bdyrho"
    var1 = rootgrp1.variables["bdypsi"][:,:]
    var2 = rootgrp2.variables["bdyrho"][:,:]

# close the netcdf file
rootgrp1.close()
rootgrp2.close()

diffvar = var2 - var1

# plot psi
fig=plt.figure()
plt.ion()

vmax = var1.max()
vmin = -vmax

ax = fig.add_subplot(3,1,1)
# mesh = ax.pcolormesh(var1,vmin=vmin, vmax=vmax)
mesh = ax.pcolormesh(var1)
ax.set_aspect('equal', 'box','C')
plt.colorbar(mesh, ax=ax)
plt.title(title1)

ax = fig.add_subplot(3,1,2)
mesh = ax.pcolormesh(var2)
ax.set_aspect('equal', 'box','C')
plt.colorbar(mesh, ax=ax)

# mesh = ax.pcolormesh(var2,vmin=vmin, vmax=vmax)

plt.title(title2)

ax = fig.add_subplot(3,1,3)
mesh = ax.pcolormesh(diffvar)
ax.set_aspect('equal', 'box','C')
# mesh = ax.pcolormesh(var2,vmin=vmin, vmax=vmax)
plt.title(title2+' - '+ title1)

fig.colorbar(mesh, ax=ax)
fig.suptitle(variable)

plt.show()
plt.savefig("figs/"+variable+'_'+title1+'_'+title2+".png")


# plt.savefig("figs/"+variable+'_'+title1+'_'+title2+".pdf")



