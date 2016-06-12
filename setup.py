#!/usr/bin/python
# -*- encoding: utf8 -*-

# python2.7 setup.py

import os
import shutil

# ------
# tests whether peptsc4py is installed and copy the proper source
# directory as qgsolver/

_qgdir='./qgsolver'
# delete directory if existing
if os.path.exists(_qgdir):
    print 'Deletes existing '+_qgdir
    shutil.rmtree(_qgdir)

# test existence of petsc4
try:
    import petsc4py
    #from petsc4py import version
    print 'petsc4py is available'
    shutil.copytree('./src_parallel/',_qgdir)

except:
    print 'petsc4py is not available, install serial code'
    shutil.copytree('./src_serial/',_qgdir)
    


