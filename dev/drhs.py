#!/usr/bin/python
# -*- encoding: utf8 -*-

import sys
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from util import outnc

rootgrp = Dataset("data/lpsi.nc", 'r')
Lpsi = rootgrp.variables["Lpsi"][0,...]
print Lpsi.shape
rootgrp.close()

rootgrp = Dataset("data/rhsinv.nc", 'r')
rhs = rootgrp.variables["rhsinv"][0,...]
print rhs.shape
rootgrp.close()

#outnc('drhs',Lpsi-rhs)
vardata=Lpsi-rhs

rootgrp = Dataset("drhs.nc", 'w', format='NETCDF4_CLASSIC', clobber=True)
# create dimensions
rootgrp.createDimension('z', vardata.shape[0])
rootgrp.createDimension('y', vardata.shape[1])
rootgrp.createDimension('x', vardata.shape[2])


 # create variables
dtype='f8'
nc_var = rootgrp.createVariable('drhs',dtype,('z','y','x'))


# fills in coordinate variables, keep truely periodic data
nc_var[:]=vardata[:]

# close the netcdf file
rootgrp.close()

