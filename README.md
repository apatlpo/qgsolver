

# Overview

This library implements a 3D z-level QG solver in regular grid or in lon/lat 
grid.
It relies on petsc4py for parallelized PV inversions but can also work in serial
(not at the moment)


# Install

## qgsolver

Download and install with:
```
git clone https://apatlpo@bitbucket.org/apatlpo/qgsolver.git
cd qgsolver
python setup.py
```

## libraries required

qgsolver requires petsc4py (and thus petsc) and netcdf4

### with conda on standard linux platform

We use conda for the install of python libraries required by qgsolver:
```
bash
source activate petsc_env
export PYTHONPATH=$PYTHONPATH:/home/slyne/aponte/natl60/python/oocgcm/
```

Proper conda install on Linux:
```
wget http://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh
(specify .miniconda2 and not miniconda2 as target dir for conda)
bash
conda update conda
conda create --name petsc_env python
source activate petsc_env
conda install -c juanlu001 petsc4py=3.6.0
conda install -y netcdf4
```

### with conda on caparmor (pb size < 512x252x100)

Proper conda install on Caparmor:
```
conda create --name petsc_env python
source activate petsc_env
conda install -c sed-pro-inria petsc4py=3.4
conda install -y netcdf4
```

Use of qgsolver on Caparmor
```
bash
source activate petsc_env
cd .../qgsolver/dev
python run_caparmor.py workdir
```
run\_caparmor.py creates "workdir" in directory /work/username with subdirectories dev and qgsolver.

The .bashrc file in the caparmor home directory could look like:
```
#alias
# User specific aliases and functions
alias qs="qstat|grep aponte"

# added by Miniconda2 4.0.5 installer
export PATH="/home1/caparmor/aponte/.miniconda2/bin:$PATH"

# add path to launch batch
export PATH="/usr/pbs/bin/:$PATH"

# for qgsolver
export WORKDIR="/work/aponte/"

source activate petsc
```

### with pip on caparmor (pb size >= 512x252x100)

When the problem is larger than 256x256x100 (approximately),
petsc needs to have been compiled with the option --with-64-bit-indices.
If not, petsc will crash systematically with errors (bus error or out of memory).

```csh
module load python/2.7.10_gnu-4.9.2
# install pip:
# https://pip.pypa.io/en/stable/installing/#id10
wget https://bootstrap.pypa.io/get-pip.py
python get-pip.py --user
# edit .cshrc:
set path = ($path $home/.local/bin)
setenv  LD_LIBRARY_PATH /home1/caparmor/aponte/.local/lib:${LD_LIBRARY_PATH}

# install target libraries
setenv MPICC mpiicc
# https://bitbucket.org/mpi4py/mpi4py/issues/53/building-with-intel-compiler-failed
pip install --user --upgrade --ignore-installed --no-cache-dir mpi4py
pip install --user --upgrade --ignore-installed --no-cache-dir numpy

setenv PETSC_CONFIGURE_OPTIONS '--with-64-bit-indices --with-fc=0 --download-f2cblaslapack'
pip install --user --upgrade --ignore-installed --no-cache-dir  petsc petsc4py
```

I should also be possible to compile petsc4py within a petsc compilation with
something like (not tested):
```
module load python/2.7.10_gnu-4.9.2
setenv MPICC mpiicc
setenv PETSC_DIR /home1/caparmor/aponte/petsc/petsc-3.7.4
setenv PETSC_ARCH linux-gnu-intel

wget http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.7.4.tar.gz
tar -zxvf petsc-3.7.4.tar.gz
cd petsc-3.7.4/

./configure PETSC_ARCH=linux-gnu-intel --with-cc=mpiicc --with-fc=mpiifort --with-blas-lapack-dir=/appli/intel/Compiler/11.1/073/mkl  --with-64-bit-indices   --download-petsc4py
make all test
```



