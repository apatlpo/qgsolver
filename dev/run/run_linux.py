#!python
# -*- encoding: utf8 -*-

# python2.7 setup.py
# python2.7 setup.py build_ext --inplace

import os
import shutil
import sys
import importlib

   
def create_dirs():  

    # get useful dirs
    startdir=os.getcwd()
    HOMEDIR = startdir+"/../.."
    USER = os.getenv("USERNAME")
    WORKDIR = HOMEDIR+"/users/"+USER

    # get args
    rundir = sys.argv[1]
    casename = sys.argv[2]
    valid_cases = ['uniform','roms','nemo']
    if not any(casename in case for case in valid_cases):
        print 'unknown case (uniform or roms or nemo)'
        sys.exit()
    goodcase=False

    RPATH = WORKDIR+'/'+rundir
    if os.path.exists(RPATH) :
        os.system('rm -Rf '+RPATH)
    os.makedirs(RPATH)
    os.chdir(RPATH)
    RPATH = os.getcwd()
    os.mkdir("input")
    os.mkdir("output")
    # os.mkdir("qgsolver")
    print "Your workking directory is "+RPATH
    return HOMEDIR,RPATH,casename

def copy_scripts(pyscript='test_basic.py'):

    os.chdir(RPATH)

    # copy python files in working directory
    shutil.copytree(HOMEDIR+'/src_parallel/','./qgsolver')    
    shutil.move("./qgsolver/omega.py", "./qgsolver/pvinv.py")
    if casename == 'nemo':
        shutil.copy(HOMEDIR+'/dev/input/create_input_nemo.py','./input')
    elif casename == 'roms':
        shutil.copy(HOMEDIR+'/dev/input/create_input_roms.py','./input')
    # shutil.copy(HOMEDIR+'/dev/run/test_basic.py','./qgsolver')
    shutil.copy(HOMEDIR+'/dev/run/test_omega.py','./qgsolver')


def read_pyscript(pyscript='test_basic.py'):
        
    # append qgsolver directory to import package
    sys.path.append(HOMEDIR)
    # Search the number of cores in test_basic.py
    mod = importlib.import_module(pyscript.replace('.py',''))
    ncores_x, ncores_y = mod.main(ping_mpi_cfg=True)
    
    #select batch queue
    nb_cores = ncores_x*ncores_y
    nb_nodes = ((nb_cores-1)/28)+1
    memory = 120
    
    return ncores_x, ncores_y, nb_cores, nb_nodes, memory
    


def submit_mpi():
    
    # make job.datarmor
    os.chdir(RPATH+'/qgsolver')
    os.system('export PYTHONPATH='+RPATH)
    os.system('mpirun -np '+str(nb_cores)+' python test_basic.py -ksp_view -ksp_monitor -ksp_converged_reason >& output.mpi')
    print 'Log file:  output.mpi in '+RPATH+'/qgsolver'    

    return
    
    
def move_input_files():
    
    # copy data file in workdir
    if casename=='roms':
        shutil.copy(HOMEDIR+'/dev/data/jet_cfg1_wp5_4km_k3.2e8_0a1500j_zlvl_pv.nc',RPATH+'/dev/data')
    elif casename=='nemo':
        shutil.copy(HOMEDIR+'/users/slgentil/input_files/nemo_metrics.nc',RPATH+'/input')
        shutil.copy(HOMEDIR+'/users/slgentil/input_files/nemo_psi.nc',RPATH+'/input')
        shutil.copy(HOMEDIR+'/users/slgentil/input_files/nemo_pv.nc',RPATH+'/input')
        shutil.copy(HOMEDIR+'/users/slgentil/input_files/nemo_rho.nc',RPATH+'/input')
	print "Copy data files"

    return



if __name__ == "__main__":


    # check number of arguments
    if  len(sys.argv) < 2:
        print '[syntaxe] : run_linux workdir case'
        print 'rundir = directory created in /work/username'
        print 'case = uniform or roms or nemo'
        quit()

    # create workking directories
    HOMEDIR,RPATH,casename = create_dirs()

    # copy_scripts
    copy_scripts(HOMEDIR,RPATH)

    # read py script
    ncores_x, ncores_y, nb_cores, nb_nodes, memory = read_pyscript(HOMEDIR,RPATH)
    
    try:
        import petsc4py
        #from petsc4py import version
        print 'petsc4py is available'    
        
        # move input files
        move_input_files()   
        
        # submit job
        submit_mpi()

    except:
        print 'petsc4py is not available, install serial code'
        


 
