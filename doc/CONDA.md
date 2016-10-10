


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
conda install netcdf4=1.1.1
conda install -c asmeurer pango 
```



Miniconda in general:

# Overview

Miniconda installers contain the conda package manager and Python.
Once Miniconda is installed, you can use the conda command to install any other packages and create environments.
There are two variants of the installer: Miniconda is Python 2 based and Miniconda3 is Python 3 based.
The other difference is that the Python 3 version of Miniconda will default to Python 3 when creating new environments and building packages
Miniconda
Installation

After downloading Miniconda2-latest-Linux-x86_64.sh
bash Miniconda2-latest-Linux-x86_64.sh

Miniconda must be used with bash. If you want to use it with csh, add in your .cshrc (pas terrible!!!)
```
#
#----------------------------------------------------------------
# alias Miniconda
#----------------------------------------------------------------
#
setenv PATH ${PATH}: /home/mulroy/slgentil/miniconda2/bin
alias source_activate 'setenv OLDPATH ${PATH};setenv PATH /home/mulroy/slgentil/miniconda2/envs/\!*/bin:${PATH}'
alias source_deactivate 'setenv PATH $OLDPATH'
```

#Main commands:
```
conda --version
conda update conda
conda create --name natl60 python Create new environment natl60
source activate natl60 Switch to another environment (activate/deactivate) (or source_activate in csh)
source deactivate To change your path from the current environment back to the root (or source_deactivate in csh)
conda info --envs List all environments
conda remove --name natl60 --all Delete an environment
conda list View a list of packages and versions installed in an environmentSearch for a package
conda search packagename Check to see if a package is available for conda to install
conda install packagename Install a new package
rm -rf /home/mulroy/slgentil/miniconda2 Remove conda
```

#Install a package from Anaconda.org

For packages that are not available using conda install, we can next look on Anaconda.org. Anaconda.org is a package management service for both public and private package repositories. Anaconda.org is a Continuum Analytics product, just like Anaconda and Miniconda.

In a browser, go to http://anaconda.org. We are looking for a package named “pestc4py”
There are more than a dozen copies of petsc4py available on Anaconda.org, first select your platform, then you can sort by number of downloads by clicking the “Downloads” heading.

Select the version that has the most downloads by clicking the package name. This brings you to the Anaconda.org detail page that shows the exact command to use to download it:

Check to see that the package downloaded
```
conda list
```

#Install a package with pip

For packages that are not available from conda or Anaconda.org, we can often install the package with pip (short for “pip installs packages”).
Exporting environment

```
conda env export > environment.yml on a machine
conda env create -f environment.yml -n $ENV_NAME on the new machine
```


