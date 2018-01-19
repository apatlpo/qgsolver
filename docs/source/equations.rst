Equations of motions, grids
===========================

Continuous form
---------------

Central state variables are the geostrophic streamfunction :math:`\psi` and potential vorticity :math:`q` 
which are related according to:

.. math::

   q(x,y,z) = f-f_0 + \Delta \psi + \partial_z \Big ( \frac{f_0^2}{N^2} \partial_z \psi \Big ) 

where :math:`f_0` is the averaged Coriolis parameter and :math:`N(z)` is the buoyancy frequency.

Density anomalies and geostrophic currents are related to the streamfunction according to:

.. math::

   \partial_z \psi &= - \frac{g\rho}{\rho_0 f_0} \\
   (u,v) &= (-\partial_y \psi, \partial_x \psi)

The evolution of the system is governed by the advection of potential vorticity and top and bottom densities
by geostrophic currents:

.. math::

   \partial_t q + J(\psi,q) + J(\Psi,q) + J(\psi,Q) &= 0 \\
   \partial_t \partial_z \psi + J(\psi,\partial_z \psi) + J(\Psi,\partial_z \psi) + J(\psi,\partial_z \Psi) &= 0 \mathrm{at} z=0,-h

where capitals represent the large scale - slowly evolving background.

Following Arakawa and Moorthi 1988, we solve for a generalized potential vorticity :math:`\tilde{q}`:

.. math::

   \tilde{q}(x,y,z) &= f-f_0 + \Delta \psi + \partial_z \Big ( \frac{f_0^2}{N^2} \partial_z \psi \Big ) - \frac{f_0^2}{N^2} \partial_z \psi \delta(z=0) + \frac{f_0^2}{N^2} \partial_z \psi \delta(z=-h) \\
   &= f-f_0 + \Delta \psi + \partial_z \Big ( \frac{f_0^2}{N^2} \partial_z \psi \Big ) + \frac{f_0}{N^2} \frac{g\rho}{\rho_0} \delta(z=0) - \frac{f_0}{N^2} \frac{g\rho}{\rho_0} \delta(z=-h)

where :math:`\delta(z=0)=1/\delta z` at :math:`z=0` (corresponds to :math:`\rho_{kup}`, see description of the vertical grid) and :math:`\delta(z=-h)=1/\delta z` at :math:`z=-h` (corresponds to :math:`\rho_{kdown}`)  

The quasi-geostrophic evolution is then solely described by the advection of :math:`\tilde{q}`:

.. math::

   \partial_t \tilde{q} + J(\psi,\tilde{q}) + J(\Psi,\tilde{q}) + J(\psi,\tilde{Q}) &= 0 


Vertical grid
-------------

The vertical grid is Charney-Phillips type, meaning streamfunction and potential vorticity are on identical
vertical levels while density is at intermediate levels.

.. image:: vgrid.png
   :scale: 70 %


Horizontal grid
---------------





