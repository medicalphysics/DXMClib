Physic model
============
DXMClib is a naïve Monte Carlo implementation of photon transport. Using this library comes down to creating a world and beam sources. Beam sources can be anything from a monochromatic pencil beam to a collimated beam rotating in a spiral (CT scan). The world consist of a 3D voxelized box model, where each voxel have properties such as density and type of material. 

This section intends to give an overview of methods and techniques used in the software and a description of the physics model. 

Geometry
--------
DXMClib has a coordinate system that closely follow the DiCOM standard for, among others, CT scanners. A right handed coordinate system with the first axis going from left to right, second down to up and last from shallow to deep. 

.. image:: ./figures/coord.png
    :width: 200
    :alt: DXMClib coordinate system

The world may be rotated in the basis coordinate system. Orientation of the world is described by two orientation direction cosines: :math:`\vec{x}, \vec{y}`. The third direction cosine is given by :math:`\vec{z} = \vec{x} \times \vec{y}` since the world basis must be orthonormal. A transformation matrix from a point to a rotated world is then :math:`M=\left[ \vec x \: \vec y \: \vec z \right]^{-1}`. This is used when simulating photon transport, while computing histories through the world the particles position and direction is given in the world basis. 

Definition of materials
-----------------------
The only material properties DXMClib cares for are mass attenuation coefficients, atomic form factors and scatter factors for chemical elements and compounds. Thankfully, Tom Schoonjans have created an excellent library providing these properties. For this reason xraylib_ is a required dependency of DXMClib. Materials can be specified by a chemical element, or a compound described by chemical symbol and number density, for example H2O or HO0.5 for water. In addition, standard NIST materials are also included. 

.. _xraylib: https://github.com/tschoonj/xraylib

Particle transport
----------------
X-ray energies above 150 keV is rarely used in diagnostic imaging, for an electron with energy 150 keV the CSDA range in soft tissue is about 0.3 mm and on par with typical voxel size in CT imaging. DXMClib assumes that all interactions creating secondary electrons positions their energy in the current voxel. All photons are described simply by a position, unit direction vector, energy and weight. Of these only the weight attribute should need an explanation. Ideally all photons will have a weight of unity or one. The weight is introduced to simplify variations in fluence from a source, for example modeling of a CT bow-tie filter. Instead of randomly sampling photons with a fluence profile mimicking a bow-tie filter a flat fluence profile can be used instead and assigning photon weight according to the bow-tie fluence profile. The expectation weight of a large number of photons is 1.0 with the added effect that same number density of photons are simulated on the filter edge as at the center of the filter. 

For efficient photon transport in a voxelized volume there are two suitable algorithms; calculating the radiologic path to compute interaction point [#SUNDERMAN1998]_ or Woodcock tracking [#WOODCOCK1965]_. While Siddons path algorithm by calculating the radiologic path through the whole volume to find an interaction point are suitable to track even a few photons it's quite inefficient compared to Woodcock tracking for large number of voxels. Woodcock tracking are perhaps best explained by introducing photon transport in a homogeneous volume first.
To sample a path length of a photon in a homogeneous volume, draw a random number :math:`r` in interval [0, 1). The photon path length is then :math:`l= -\ln(r)/(\mu \rho)` where :math:`\mu` and :math:`\rho` is mass attenuation coefficient and density of the material. This can be extended to a heterogeneous volume by introducing virtual interactions. The path length is calculated in a similar way: :math:`l= -\ln(r)/\zeta` with

.. math::
    \zeta = \max_i \left( \mu_i \rho_i \right)

for each volume element :math:`i`. For the material corresponding to the traversed step :math:`l`, an interaction happens if 

.. math::
    \frac{\mu \rho}{\zeta} \leqslant r

where :math:`r` is a random number in interval [0, 1), otherwise the interaction is virtual (i.e nothing happens) and a new step is calculated. This method is numerical effective but only valid when a large number of photons are simulated, and is unsuitable to showcase one particle track. 

When an interaction occurs type of interaction is sampled by drawing a new random number :math:`r` in interval [0, :math:`\mu`).

.. math::
    \mu = \mu_{photoelectric} + \mu_{incoherent} + \mu_{coherent}

if :math:`r < \mu_{photoelectric}` a photoelectric event happens, if :math:`r < \mu_{photoelectric} + \mu_{incoherent}` an incoherent scattering event happens and else a coherent scattering event. 

Forced interactions can be used as a variance reduction technique, for example calculating dose in an air chamber inside an CTDI phantom. Since air is thankfully a low density material, relative few interactions occurs compared to in water or plastic. Forcing interactions in a voxel can decrease the variance of a dose calculation in a voxel without spending more CPU cycles to simulate more histories. Forced interactions are implemented such that when a photon steps into a voxel with a forced interaction flag an photoelectric event is scored, even if the interaction is considered virtual. To balance out the scored energy only a fraction 

.. math::
    \frac{\mu_{photoelectric}\rho}{\zeta}

of the photon energy is scored. The weight of the photon is reduced accordingly by 

.. math::
    w_{after} = w_{before}(1-\frac{\mu_{photoelectric}\rho}{\zeta})

Further transport of the photon is done by the ordinary method by sampling a random number to determine if an interaction occurs and determine if a coherent or noncoherent event will happen (photoelectric effect is already dealt with).

Generating x-ray spectra
------------------------
Since most diagnostic x-ray units do not emit monochromatic photon beams this library includes a x-ray spectra generator. The implementation uses a semi-analytical model proposed by Poludniowski [#Poludniowski1]_ [#Poludniowski2]_ for simulating a spectra from a pure tungsten anode. The model is valid tube potentials from 50kVp to 150 kVp but is accurate up to 300 kVp. The implementation allows for adding filtration of any material and to freely select tube potential and anode angle proving quite flexible. Since the model requires an evaluation of a double integral for each energy bin which is quite computational expensive this implementation is multi threaded. The same formalism is also used in the SpekCalc application also by Poludniowski et al. 

Sampling photon energies from a specter is implemented by the squaring of histogram method which is quite fast after an initial generation of a lookup table. When an energy bin is sampled the photon energy is finally uniformly sampled within the bin width. 

Photon transport
----------------
Photon transport in DXMClib is implemented in a relatively simple manner. A source will set up one or multiple exposures where an exposure is emitting photons from a fixed point and a fixed beam direction.  A photon is created at the exposure (tube) position and the direction is sampled uniformly inside the collimation. The photon energy is either sampled from a specter or if the source is monochrome, given the selected monochrome energy. The weight of the photon is calculated based on direction and any selected filters, such as a CT bow tie filter or a Heel effect filter or both.

The sampled photon is first checked for intersecting the voxel volume, also called the world. If it intersects, it is transported to the world border before the Woodcock tracking starts. 

Photoelectric effect
____________





References
----------
.. [#SUNDERMAN1998] A Fast Algorithm to Calculate the Exact Radiological Path Through a Pixel Or Voxel Space, Sunderman E. et al. Journal of Computing and Information Technology 6(1). December 1998.
.. [#WOODCOCK1965] Woodcock E.R. et al. Techniques used in the GEM code for Monte Carlo neutronics calculations in reactors and other systems of complex geometry. ANL-7050. Argonne National Laboratory, 1965.
.. [#Poludniowski1] Poludniowski, G.G. and Evans, P.M. (2007), Calculation of x‐ray spectra emerging from an x‐ray tube. Part I. Electron penetration characteristics in x‐ray targets. Med. Phys., 34: 2164-2174. doi:10.1118/1.2734725
.. [#Poludniowski2] Poludniowski, G.G. (2007), Calculation of x‐ray spectra emerging from an x‐ray tube. Part II. X‐ray production and filtration in x‐ray targets. Med. Phys., 34: 2175-2186. doi:10.1118/1.2734726