#!/usr/bin/python
# -*- encoding: utf8 -*-

import sys
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt

# open netdc file
rootgrp = Dataset("data/output.nc", 'r')

# loop around variables to load
psi = rootgrp.variables["psi"][:]
print psi.shape
dpsidy = np.diff(psi,axis=2)/4000.
print dpsidy.shape,dpsidy[0,-1,:,:].min(),dpsidy[0,-1,:,:].max()

# close the netcdf file
rootgrp.close()

plt.figure()
plt.ion()
plt.show()
plt.pcolormesh(dpsidy[0,-1,:,:])
plt.colorbar()
plt.show()
plt.savefig("figs/dpsidy.pdf")