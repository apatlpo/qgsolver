#!/usr/bin/env python

from setuptools import setup

setup(name='qgsolver',
      version='0.0.0',
      description='Quasi-Geostrophic solver',
      url='https://github.com/apatlpo/qgsolver',
      author='Aurelien Ponte',
      author_email='aurelien.junk@gmail.com',
      license='',
      packages=['qgsolver'],
      install_requires=['petsc4py', 'numpy', 'scipy']
      )
