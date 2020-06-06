
DXMClib
===================================
DXMClib (Diagnostic X-ray Monte Carlo) is a radiation dose scoring library for diagnostic photon energies in voxelized geometry written in C++. The main goal for this library is to provide an accurate enough physics model to describe and model x-ray sources and estimate radiation doses by the Monte-Carlo method.

If you are looking for an application with graphical user interface to perform simulations, OpenDXMC_ uses this library as Monte Carlo engine and also allows for import CT images and phantoms as scoring volumes. 

.. _OpenDXMC: https://github.com/medicalphysics/OpenDXMC/releases

Another header for me
=====================

.. toctree::
   :maxdepth: 2
   about
   license
   
Docs
====
.. doxygenclass:: CTDualSource
   :members:

Indices and tables
==================

* :ref:`genindex`
* :ref:`search`
