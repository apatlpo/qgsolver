#!python
# -*- encoding: utf8 -*-

# python2.7 setup.py
# python2.7 setup.py build_ext --inplace

import os
import shutil
import sys

# check number of arguments
if  len(sys.argv) < 2:
    print '[syntaxe] : run_caparmor workdir case'
    print 'workdir = directory created in /work/username'
    print 'case = roms or nemo'
    quit()

workdir = sys.argv[1]
casename = sys.argv[2]
valid_cases = ['uniform','roms','nemo']
if not any(casename in case for case in valid_cases):
    print 'unknown case (uniform or roms or nemo)'
    sys.exit()
goodcase=False

# Search the number of cores in test_basic.py
pyfile = open( 'test_basic.py', 'r' )
for line in pyfile:
    if 'casename' in line and  casename in line:
        goodcase=True
    if 'ncores_x' in line and goodcase==True :
        ncores_x = int(line[line.index('=')+1:])
    if 'ncores_y' in line and goodcase==True :
        ncores_y = int(line[line.index('=')+1:])
        break

#select batch queue
nb_cores = ncores_x*ncores_y
nb_nodes = ((nb_cores-1)/8)+1
if nb_cores<=8:
    queue='parallel8'
elif nb_cores<=32:
    queue='parallel32'
elif nb_cores<=64:
    queue='parallel64'
elif nb_cores<=256:
    queue='parallel256'

# Création du répertoire workdir dans /work
startdir=os.getcwd()
HOMEDIR = startdir+"/.."
WORKDIR = os.getenv("WORKDIR")
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

    # copy python files in workdir
    shutil.copytree(HOMEDIR+'/src_parallel/','./qgsolver')
    os.mkdir(RPATH+'/dev')
    os.mkdir(RPATH+'/dev/data')
    shutil.copy(HOMEDIR+'/dev/test_basic.py','./dev')


    # make job.caparmor
    os.chdir(RPATH+'/dev')
    fo = open('job_caparmor','w')
    fo.write('#!/bin/csh\n')
    fo.write('#PBS -N qgsolver\n')
    fo.write('#PBS -q '+queue+'\n')
    fo.write('#PBS -l select='+str(nb_nodes)+':ncpus=8:mpiprocs=8\n')
    fo.write('\n')
    fo.write('# cd to the directory you submitted your job\n')
    fo.write('cd $PBS_O_WORKDIR\n')
    fo.write('\n')
    fo.write('# get the path for python\n')
    fo.write('setenv PATH ${HOME}/.miniconda2/envs/petsc/bin:${PATH}\n')
    fo.write('setenv PYTHONPATH $PBS_O_WORKDIR/..\n')
    fo.write('\n')
    # fo.write('time mpirun -np '+str(nb_cores)+' python test_basic.py  >& output.mpi\n')
    fo.write('time mpirun -np '+str(nb_cores)+' python test_basic.py -ksp_view -ksp_monitor -ksp_converged_reason >& output.mpi\n')
    fo.close()

    # copy data file in workdir
    if casename=='roms':
        shutil.copy(HOMEDIR+'/dev/data/jet_cfg1_wp5_4km_k3.2e8_0a1500j_zlvl_pv.nc',RPATH+'/dev/data')
    elif casename=='nemo':
        shutil.copy(HOMEDIR+'/dev/data/nemo_metrics.nc',RPATH+'/dev/data')
        shutil.copy(HOMEDIR+'/dev/data/nemo_psi.nc',RPATH+'/dev/data')
        shutil.copy(HOMEDIR+'/dev/data/nemo_pv.nc',RPATH+'/dev/data')
        shutil.copy(HOMEDIR+'/dev/data/nemo_rho.nc',RPATH+'/dev/data')

    os.system('cd '+RPATH+'/dev')
    os.system('qsub job_caparmor')
    print 'Log file:  output.mpi in '+RPATH+'/dev'

except:
    print 'petsc4py is not available, install serial code'
    shutil.copytree(HOMEDIR+'/src_serial/','./qgsolver')

