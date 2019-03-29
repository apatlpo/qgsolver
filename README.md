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

### Install dependencies with conda

The principal dependcy is [petsc4py](https://bitbucket.org/petsc/petsc4py)

Download Miniconda3 from the [conda website](https://conda.io/miniconda.html)
```csh
bash Miniconda3-latest-Linux-x86_64.sh
(specify .miniconda3 and not miniconda3 as target dir for conda)
bash
conda update conda
conda create -n petsc -c conda-forge python=3.6 petsc4py netcdf4 matplotlib snakeviz xarray
source activate petsc
```

# Run on datarmor

Generate input files if necessary on your workstation
```csh
cd qgsolver/input/
python create_input_roms.py -o roms_in/
```

On datarmor, copy input files if necessary:
```csh
cp -r /home/slyne/aponte/natl60/qgsolver/input/roms_in /home1/datawork/aponte/qgsolver/
```

Run qgsolver on Datarmor:
```csh
bash
source activate petsc
cd /home/slyne/aponte/natl60/qgsolver/run/
python run_datarmor.py qgsolver/roms_out roms.py /home1/datawork/aponte/qgsolver/roms_in/
```

# API

See [doc](http://qgsolver-doc.readthedocs.io/en/latest/)

