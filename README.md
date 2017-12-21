# Overview

This library implements a 3D z-level QG solver in regular grid or in lon/lat 
grid.
It relies on petsc4py for parallelized PV inversions but should also work in serial
in the future


# Install

## qgsolver

Download and install with:
```csh
git clone https://github.com/apatlpo/qgsolver.git
cd qgsolver
```

## libraries required

qgsolver requires [petsc4py](https://bitbucket.org/petsc/petsc4py) (and thus petsc) and netcdf4

### Install with conda 

Download Miniconda3 from the [conda website](https://conda.io/miniconda.html)
```csh
bash Miniconda3-latest-Linux-x86_64.sh
(specify .miniconda3 and not miniconda3 as target dir for conda)
bash
conda update conda
conda create --name petsc python
source activate petsc
conda install -c conda-forge petsc
conda install -c conda-forge petsc4py
conda install -c conda-forge netcdf4
conda install -c conda-forge matplotlib
conda install -c conda-forge jupyter 
```

Use of qgsolver on Datarmor
```csh
bash
source activate petsc
cd .../qgsolver/dev
python run_caparmor.py workdir
```
run\_caparmor.py creates "workdir" in directory /work/username with subdirectories dev and qgsolver.


