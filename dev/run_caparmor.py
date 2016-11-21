#!python
# -*- encoding: utf8 -*-

# python2.7 setup.py
# python2.7 setup.py build_ext --inplace

import os
import shutil
import sys

# check number of arguments

if  len(sys.argv) < 2:
    print '[syntaxe] : run_caparmor workdir'
    quit()

workdir=sys.argv[1]

# Création du répertoire workdir dans /work

startdir=os.getcwd()
HOMEDIR = "/home/mulroy/slgentil/models/natl60/qgsolver"
WORKDIR = os.getenv("workdir")
RPATH = WORKDIR+'/'+workdir

if os.path.exists(RPATH) :
    os.system('rm -Rf '+RPATH)
os.mkdir(RPATH)
os.chdir(RPATH)


# test existence of petsc4py
try:
    import petsc4py
    #from petsc4py import version
    print 'petsc4py is available'
    shutil.copytree(HOMEDIR+'/src_parallel/','./qgsolver')
    os.mkdir(RPATH+'/dev')
    os.mkdir(RPATH+'/dev/data')
    shutil.copy(HOMEDIR+'/dev/test_basic.py','./dev')
    shutil.copy(HOMEDIR+'/dev/job_caparmor','./dev')
    shutil.copy(HOMEDIR+'/dev/data/jet_cfg1_wp5_4km_k3.2e8_0a1500j_zlvl_pv.nc','./dev/data')

    print 'cd '+RPATH+'/dev'
    print 'qsub job_caparmor'

except:
    print 'petsc4py is not available, install serial code'
    shutil.copytree(HOMEDIR+'/src_serial/','./qgsolver')

