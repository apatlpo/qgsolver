# To build the fast modules, run
#  python setup.py build_ext --inplace

from distutils.core import setup, Extension
import numpy

# define the extension module
fast_interp3D = Extension('fast_interp3D', sources=['fast_interp3D.c'],
                          extra_compile_args = ['-fopenmp'],
                          extra_link_args = ['-fopenmp'],
                          include_dirs=[numpy.get_include()])

# run the setup
setup(ext_modules=[fast_interp3D])

