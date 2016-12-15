
This library implements a 3D z-level QG solver in regular grid or in lon/lat 
grid.
It relies on petsc4py for parallelized PV inversions but can also work in serial
(not at the moment)

Download and install with:
```
git clone https://apatlpo@bitbucket.org/apatlpo/qgsolver.git
cd qgsolver
python setup.py
```


We use conda for the python install:
```
bash
source activate natl60
export PYTHONPATH=$PYTHONPATH:/home/slyne/aponte/natl60/python/oocgcm/
```

Proper conda install on Linux:
```
bash /home/mulroy/slgentil/tarfiles/Miniconda2-latest-Linux-x86_64.sh
bash
conda update conda
conda create --name natl60 python
source activate natl60
conda install dask
conda install xarray
conda install -c juanlu001 petsc4py=3.6.0
conda install libgfortran=1.0
conda install -c scitools cartopy=0.13.1
conda install basemap
conda install -c asmeurer pango
```

Proper conda install on Caparmor:
```
conda create --name petsc python
source activate petsc
conda install -c sed-pro-inria petsc4py=3.4
conda install -y netcdf4
```

Use on Caparmor
```
bash
source activate petsc
cd .../qgsolver/dev
python run_caparmor.py workdir

The script creates "workdir" in directory /work/username with subdirectories dev and qgsolver.

```

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


