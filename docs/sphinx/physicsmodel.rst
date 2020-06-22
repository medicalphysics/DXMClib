Physic model
============
DXMClib is a naïve Monte Carlo implementation of photon transport. Using this library comes down to creating a world and beam sources. Beam sources can be anything from a monochromatic pecilbeam to a collimated beam rotating in a spiral (CT scan). The world consist of a 3D voxalized box model, where each voxel have properties such as density and type of material. 

This section intends to give an overview of methods and techniques used in the software and a description of the physics model. 

Geometry
--------
DXMClib has a coordinate system that closely follow the DiCOM standard for, among others, CT scanners. A right handed coordinate system with the first axis going from left to right, second down to up and last from shallow to deep. 

.. image:: ./figures/coord.png
    :width: 200
    :alt: DXMClib coordinate system

The world may be rotated in the basis coordinate system. Orientation of the world is described by two orientation direction cosines: :math:`\vec{x}, \vec{y}`. The third direction cosine is given by :math:`\vec{z} = \vec{x} \times \vec{y}` since the world basis must be orthonormal. A transformation matrix from a point to a rotated world is then :math:`M=\left[ \vec x \: \vec y \: \vec z \right]^{-1}`. This is used when simulating photon transport, while computing histories throught the world the particles position and direction is given in the world basis. 

Definition of materials
-----------------------
The only material properties DXMClib cares for are mass attenuation coefficients, atomic form factors and scatter factors for chemical elements and compounds. Thankfully, Tom Schoonjans have created an excellent library providing these properties. For this reason xraylib_ is a required dependency of DXMClib. Materials can be specified by a chemical element, or a compound described by chemical symbol and number density, for example H2O or HO0.5 for water. In addition, standard NIST materials are also included. 

.. _xraylib: https://github.com/tschoonj/xraylib

Particle transport
----------------
X-ray energies above 150 keV is rareley used in diagnostic imaging, for an electron with energy 150 keV the CSDA range in soft tissue is about 0.3 mm and on par with typical voxel size in CT imaging. DXMClib assumes that all interactions creating secondary electrons positions their energy in the current voxel. 
For efficient photon transport in a voxelized volume there are two suitable algorithms; calculating the radiologic path to compute interaction point [#SUNDERMAN1998]_ or Woodcock tracking [#WOODCOCK1965]_. While Siddons path algorithm by calculating the radiologic path thorugh the whole volume to find an interaction point are suitable to track even a few photons it's quite inefficient compared to Woodcock tracking for large number of voxels. Woodcock tracking are perhaps best explained by introducing photon transport in a homogeneous volume first.
To sample a path length of a photon in a homogeneous volume, draw a random number :math:`r` in interval [0, 1). The photon path length is then :math:`l= -\ln(r)/(\mu \rho)` where :math:`\mu` and :math:`\rho` is mass attenuation coefficient and density of the material. This can be extended to a heterogeneous volume by introducing virtual interactions. The path lenght is calculated in a similar way: :math:`l= -\ln(r)/\zeta` with

..math::
    \zeta = \max\left[ \mu_i \rho_i \right]

For the material corresponding to the traversed step :math:`l`, an interaction happends if 
..math::
    \frac{\mu\rho}{\zeta} < r
where :math:`r` is a random number in interval [0, 1), otherwise the interaction is virtual and a new step is calculated.

.. [#SUNDERMAN1998] A Fast Algorithm to Calculate the Exact Radiological Path Through a Pixel Or Voxel Space, Sunderman E. et al. Journal of Computing and Information Technology 6(1). December 1998.
.. [#WOODCOCK1965] Woodcock E.R. et al. Techniques used in the GEM code for Monte Carlo neutronics calculations in reactors and other systems of complex geometry. ANL-7050. Argonne National Laboratory, 1965.