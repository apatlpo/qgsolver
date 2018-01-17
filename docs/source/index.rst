qgsolver_doc's documentation
========================================

Install
--------------------

We recommend conda for dependencies, see README on `qgsolver github repository`__

.. _qgsolvergt: https://github.com/apatlpo/qgsolver

__ qgsolvergt_


Tutorial on a desktop
--------------------

Tutorial on datarmor
--------------------

PV inversion solver
--------------------

We use `PETSc <http://www.mcs.anl.gov/petsc/documentation/index.html>`_ in order to solver PV or omega equation inversion.
The PETSc manual `here <http://www.mcs.anl.gov/petsc/petsc-current/docs/manual.pdf>`_ is very useful.

In order to get details about the PV inversion solver, add the following options at run time:

.. code:: bash

   $  mpirun -n 4 python test_analytical.py -mf -ksp_view -ksp_monitor -ksp_converged_reason

Creating input files
--------------------

`dev/input/ <https://github.com/apatlpo/qgsolver/blob/master/dev/input/>`_ is the relevant folder.

For ROMS simulations, you may be inspired to look at `/dev/input/create_input_roms.py <https://github.com/apatlpo/qgsolver/blob/master/dev/input/create_input_roms.py>`_

For ROMS simulations, you may be inspired to look at `/dev/input/create_input_nemo.py <https://github.com/apatlpo/qgsolver/blob/master/dev/input/create_input_nemo.py>`_


API
--------------------

.. toctree::
   :maxdepth: 2

   api/qgsolver

Indices and tables
--------------------

* :ref:`genindex`
* :ref:`search`
