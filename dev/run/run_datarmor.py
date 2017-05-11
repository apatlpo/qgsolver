#!python
# -*- encoding: utf8 -*-

# python2.7 setup.py
# python2.7 setup.py build_ext --inplace

import os
import shutil
import sys
import importlib

        
def read_pyscript(pyscript='test_basic.py'):
        
    # append qgsolver directory to import package
    sys.path.append(RPATH)
    # Search the number of cores in test_basic.py
    mod = importlib.import_module(pyscript.replace('.py',''))
    ncores_x, ncores_y = mod.main(ping_mpi_cfg=True)
    
    #select batch queue
    nb_cores = ncores_x*ncores_y
    nb_nodes = ((nb_cores-1)/28)+1
    memory = 120
    print ncores_x, ncores_y, nb_cores, nb_nodes
    return ncores_x, ncores_y, nb_cores, nb_nodes, memory
    


def copy_scripts():
    
    if os.path.exists(RPATH) :
        os.system('rm -Rf '+RPATH)
    os.mkdir(RPATH)
    os.chdir(RPATH)
    
    # copy python files in workdir
    shutil.copytree(HOMEDIR+'/src_parallel/','./qgsolver')
    os.mkdir(RPATH+'/input')
    os.mkdir(RPATH+'/output')
    # shutil.copy(HOMEDIR+'/dev/run/test_basic.py','./dev')
    shutil.copy(HOMEDIR+'/dev/run/test_omega.py','./qgsolver')

def write_batchfile():
    # make job.datarmor
    os.chdir(RPATH+'/qgsolver')
    fo = open('job_datarmor','w')
    fo.write('#!/bin/csh\n')
    fo.write('#PBS -q mpi\n')
    fo.write('#PBS -l select='+str(nb_nodes)+':ncpus=28:mpiprocs=28:mem='+str(memory)+'G\n')
    fo.write('#PBS \n')
    fo.write('#PBS -N qgsolver\n')
    fo.write('#PBS -l walltime=01:00:00\n')
    fo.write('# cd to the directory you submitted your job\n')
    fo.write('cd $PBS_O_WORKDIR\n')
    fo.write('\n')
    fo.write('# get the path for python\n')
    fo.write('setenv PATH ${HOME}/.miniconda2/envs/petsc/bin:${PATH}\n')
    fo.write('setenv PYTHONPATH $PBS_O_WORKDIR/..\n')
    fo.write('\n')
    # fo.write('time mpirun -np '+str(nb_cores)+' python test_basic.py  >& output.mpi\n')
    # fo.write('time mpirun -np '+str(nb_cores)+' python test_basic.py -ksp_view -ksp_monitor -ksp_converged_reason >& output.mpi\n')
    fo.write('time mpirun -np '+str(nb_cores)+' python test_omega.py -ksp_view -ksp_monitor -ksp_converged_reason >& output.mpi\n')
    fo.close()
    fo.close()
    
    return
    
    
def move_input_files():
    
    # copy data file in workdir
    if casename=='roms':
        shutil.copy(HOMEDIR+'/dev/data/jet_cfg1_wp5_4km_k3.2e8_0a1500j_zlvl_pv.nc',RPATH+'/input')
    elif casename=='nemo':
        #shutil.copy(HOMEDIR+'/dev/data/nemo_metrics.nc',RPATH+'/dev/data')
        #shutil.copy(HOMEDIR+'/dev/data/nemo_psi.nc',RPATH+'/dev/data')
        #shutil.copy(HOMEDIR+'/dev/data/nemo_pv.nc',RPATH+'/dev/data')
        #shutil.copy(HOMEDIR+'/dev/data/nemo_rho.nc',RPATH+'/dev/data')
        os.system('ln -s /home1/datawork/slgentil/nemo_mask_data/*.nc '+RPATH+'/input')
    return



if __name__ == "__main__":


    # check number of arguments
    if  len(sys.argv) < 2:
        print '[syntaxe] : run_caparmor workdir case'
        print 'rundir = directory created in /work/username'
        print 'case = uniform or roms or nemo'
        quit()
    
    # get useful dirs
    startdir=os.getcwd()
    HOMEDIR = startdir+"/../.."
    WORKDIR = os.getenv("DATAWORK")
    # WORKDIR = os.getenv("SCRATCH")
    # get args
    rundir = sys.argv[1]
    casename = sys.argv[2]
    valid_cases = ['uniform','roms','nemo']
    if not any(casename in case for case in valid_cases):
        print 'unknown case (uniform or roms or nemo)'
        sys.exit()
    goodcase=False
    RPATH = WORKDIR+'/'+rundir

    # copy python scripts in workdir
    copy_scripts()

    # read py script
    ncores_x, ncores_y, nb_cores, nb_nodes, memory = read_pyscript(pyscript='test_omega.py')
    
    try:
        import petsc4py
        #from petsc4py import version
        print 'petsc4py is available'    
    
        # write batch file
        write_batchfile()
        
        # move input files
        move_input_files()
        
        # submit job
        os.system('cd '+RPATH+'/qgsolver')
        # os.system('qsub job_datarmor')
        print 'cd '+RPATH+'/qgsolver'
        print 'qsub job_datarmor'
        print 'Log file:  output.mpi'    

    except:
        print 'petsc4py is not available, install serial code'
        


 
