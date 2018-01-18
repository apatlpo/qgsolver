qgsolver_doc's documentation
========================================

Install
--------------------

We recommend conda for dependencies, see README on
`qgsolver github repository <https://github.com/apatlpo/qgsolver>`_

Tutorial on a desktop
--------------------

Tutorial on datarmor
--------------------

PV inversion solver
--------------------

We use `PETSc <http://www.mcs.anl.gov/petsc/documentation/index.html>`_
in order to solver PV or omega equation inversion.
The `PETSc manual <http://www.mcs.anl.gov/petsc/petsc-current/docs/manual.pdf>`_ is very useful.

The petsc4py `documentation <http://www.mcs.anl.gov/petsc/petsc4py-current/docs/>`_ and
`API <http://www.mcs.anl.gov/petsc/petsc4py-current/docs/apiref/index.html>`_ may also
be helpful.

In order to get details about the PV inversion solver, add the following options at run time:

.. code:: bash

   mpirun -n 4 python analytical.py -mf -ksp_view -ksp_monitor -ksp_converged_reason

In order to profile with `snakeviz <https://jiffyclub.github.io/snakeviz/>`_
you need first to generate a profile and then run snakeviz:

.. code: bash

   mpirun -n 4 python -m cProfile -o output.prof uniform.py
   snakeviz output.prof

Creating input files
--------------------

`input/ <https://github.com/apatlpo/qgsolver/blob/master/input/>`_ is the relevant folder.

For ROMS simulation outputs, you may be inspired to look at
`input/create_input_roms.py <https://github.com/apatlpo/qgsolver/blob/master/input/create_input_roms.py>`_

For NEMO simulation outputs, you may be inspired to look at
`input/create_input_nemo.py <https://github.com/apatlpo/qgsolver/blob/master/input/create_input_nemo.py>`_


API
--------------------

.. toctree::
   :maxdepth: 2

   api/qgsolver

Indices and tables
--------------------

* :ref:`genindex`
* :ref:`search`
