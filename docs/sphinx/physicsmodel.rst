Physic model
============
DXMClib is a naïve Monte Carlo implementation of photon transport. Using this library comes down to creating a world and beam sources. Beam sources can be anything from a monochromatic pecilbeam to a collimated beam rotating in a spiral (CT scan). The world consist of a 3D voxalized box model, where each voxel have properties such as density and type of material. 

This section intends to give an overview of methods and techniques used in the software and a description of the physics model. 

Geometry
--------
DXMClib has a coordinate system that closely follow the DiCOM standard for, among others, CT scanners. A right handed coordinate system with the first axis going from left to right, second down to up and last from shallow to deep. 

.. image:: ./figures/coord.png
    :width: 400
    :alt: DXMClib coordinate system

The world may be rotated in the basis coordinate system. Orientation of the world is described by two orientation direction cosines: :math:`\vec{x}, \vec{y}`. The third direction cosine is given by :math:`\vec{z} = \vec{x} \times \vec{y}` since the world basis must be orthonormal. A transformation matrix from a point to a rotated world is then :math:`M=\left[ \vec x \: \vec y \: \vec z \right]^{-1}`. This is used when simulating photon transport, while computing histories throught the world the particles position and direction is given in the world basis. 

Definition of materials
-----------------------
The only material properties DXMClib cares for are mass attenuation coefficients, atomic form factors and scatter factors for chemical elements and compounds. Thankfully, Tom Schoonjans have created an excellent library providing these properties. For this reason xraylib_ is a required dependency of DXMClib. Materials can be specified by a chemical element, or a compound described by chemical symbol and number density, for example H2O or HO0.5 for water. In addition, standard NIST materials are also included. 

.. _xraylib: https://github.com/tschoonj/xraylib

Particle transport
----------------
X-ray energies above 150 keV is rareley used in diagnostic imaging, for an electron with energy 150 keV the CSDA range in soft tissue is about 0.3 mm and on par with typical voxel size in CT imaging. DXMClib assumes that all interactions creating secondary electrons positions their energy in the current voxel. 
For effective photon transport in a voxelized volume there are two suitable algorithms; calculating the radiologic path to compute interaction point or Woodcock tracking. 