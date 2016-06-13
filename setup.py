#!/usr/bin/python
# -*- encoding: utf8 -*-

# python2.7 setup.py
# python2.7 setup.py build_ext --inplace

import os
import shutil

# for cython code
# should test the existence of cython !!!
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

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
    
# ext_modules = [Extension('qgsolver.set_L_fast', 
#                          sources=['qgsolver/set_L_fast.pyx'])]
# 
# setup(cmdclass = {'build_ext': build_ext}, 
#       ext_modules = ext_modules)



