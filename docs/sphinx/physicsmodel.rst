Physic model
============
DXMClib is a naïve Monte Carlo implementation of photon transport. Using this library comes down to creating a world and beam sources. Beam sources can be anything from a monochromatic pecilbeam to a collimated beam rotating in a spiral (CT scan). The world consist of a 3D voxalized box model, where each voxel have properties such as density and type of material. 

This section intends to give an overview of methods and techniques used in the software and a description of the physics model. 

Geometry
--------
DXMClib has a coordinate system that closely follow the DiCOM standard for, among others, CT scanners. A right handed coordinate system with the first axis going from left to right, second down to up and last from shallow to deep.   

.. image ./figures/coord.png
    :width: 400
    :alt: DXMClib coordinate system


