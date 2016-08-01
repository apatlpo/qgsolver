
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

Proper conda install:
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


