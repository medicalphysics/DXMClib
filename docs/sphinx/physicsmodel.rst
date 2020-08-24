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
------------------
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
    \mu = \mu_{photoelectric} + \mu_{compton} + \mu_{rayleight}

if :math:`r < \mu_{photoelectric}` a photoelectric event happens, if :math:`r < \mu_{photoelectric} + \mu_{compton}` an Compton scattering event happens and else a Rayleight scattering event. 

Forced interactions can be used as a variance reduction technique, for example calculating dose in an air chamber inside an CTDI phantom. Since air is thankfully a low density material, relative few interactions occurs compared to in water or plastic. Forcing interactions in a voxel can decrease the variance of a dose calculation in a voxel without spending more CPU cycles to simulate more histories. Forced interactions are implemented such that when a photon steps into a voxel with a forced interaction flag an photoelectric event is scored, even if the interaction is considered virtual. To balance out the scored energy only a fraction 

.. math::
    \frac{\mu_{photoelectric}\rho}{\zeta}

of the photon energy is scored. The weight of the photon is reduced accordingly by 

.. math::
    w_{after} = w_{before}(1-\frac{\mu_{photoelectric}\rho}{\zeta})

Further transport of the photon is done by the ordinary method by sampling a random number to determine if an interaction occurs and determine if a Rayleight or Compton event will happen (photoelectric effect is already dealt with).

Generating x-ray spectra
------------------------
Since most diagnostic x-ray units do not emit monochromatic photon beams this library includes a x-ray spectra generator. The implementation uses a semi-analytical model proposed by Poludniowski [#Poludniowski1]_ [#Poludniowski2]_ for simulating a spectra from a pure tungsten anode. The model is valid tube potentials from 50kVp to 150 kVp but is accurate up to 300 kVp. The implementation allows for adding filtration of any material and to freely select tube potential and anode angle proving quite flexible. Since the model requires an evaluation of a double integral for each energy bin which is quite computational expensive this implementation is multi threaded. The same formalism is also used in the SpekCalc application also by Poludniowski et al. The Heel effect is also modelled for collimated beams along the anode cathode direction, and is equivalent to a corresponding change in anode angle in the model proposed by Poludniowski.  

Sampling photon energies from a specter is implemented by the squaring of histogram method which is quite fast after an initial generation of a lookup table. When an energy bin is sampled the photon energy is finally uniformly sampled within the bin width. 

Photon transport
----------------
Photon transport in DXMClib is implemented in a relatively simple manner. A source will set up one or multiple exposures where an exposure is emitting photons from a fixed point and a fixed beam direction.  A photon is created at the exposure (tube) position and the direction is sampled uniformly inside the collimation. The photon energy is either sampled from a specter or if the source is monochrome, given the selected monochrome energy. The weight of the photon is calculated based on direction and any selected filters, such as a CT bow tie filter or a Heel effect filter or both.

The sampled photon is first checked for intersecting the voxel volume, also called the world. If it intersects, it is transported to the world border before the Woodcock tracking starts. 

Photoelectric effect
____________
This is the simplest of three types of interactions handled by DXMClib. When a photoelectric event is triggered the photon transfers all it's energy to the voxel. The energy from a scattered electron and any photons from bremsstrahlung is assumed not to escape the voxel.

Compton scattering
__________________
Compton events are handled by sampling the Klein Nishina differential cross section for an unbound electron:

.. math::
    \frac{d\rho}{d\epsilon} = \pi r_e^2\frac{m_ec^2}{E_0}Z\left[\frac{1}{\epsilon}+\epsilon \right] \left[ 1-\frac{\epsilon \sin^2\theta}{1+\epsilon^2} \right]

with :math:`r_e` as the classical electron radius, :math:`m_ec^2`: electron mass, :math:`E_0` and :math:`E_1` as energy of incident and scattered photon respectivly, and :math:`\epsilon` as :math:`E_1/E_0`. Scatter angle :math:`\theta` is given by the Compton formula:

.. math::
    \epsilon = \frac{m_ec^2}{m_ec^2 + E_0(1-\cos\theta)}

The minimum and maximum values for :math:`\epsilon` follows from the compton formula with 

.. math::
    \epsilon_{min} = \frac{m_ec^2}{m_ec^2 +2E_0}

and 

.. math::
    \epsilon_{max} = \frac{m_ec^2}{m_ec^2} = 1

so :math:`\epsilon \in [\epsilon_{min}, 1]`. For low photon energies, i.e typical diagnostic energy levels, it's most efficient to uniformly sample :math:`\epsilon` with the rejection function: 

.. math::
    g = \frac{1}{g_{max}} \left( \frac{1}{\epsilon} + \epsilon -\sin^2\theta \right)

with

.. math::
    g_{max} = \frac{1}{\epsilon_{min}}+\epsilon_{min}

To sample the Klein-Nishina cross section an :math:`\epsilon` is uniformly sampled by 

.. math::
    \epsilon = r_1+(1-r_1)\epsilon_{min}

where :math:`r_1` is a random number in interval :math:`[0, 1]`. For the sampled :math:`\epsilon` calculate :math:`g` and :math:`\theta`. Draw a new random number :math:`r_2` in interval :math:`[0,1]`, if :math:`r_2 \leqslant g` accept the sampled :math:`\epsilon` (and :math:`\theta`) else repeat the process. 

The sampling methods described above ignores any binding effects on the electron and will overestimate forward scattering for low energy photons. DXMClib can use a simplified model (the Livermore model) for low energy correction and is enabled by default by CMake option DXMC_USE_LOW_ENERGY_COMPTON. This correction takes into account Hubbel's atomic form factor [#Hubbell]_. In this case the sampling is performed by the same procedure as a free electron except for a slighly modified rejection function:

.. math::
    g = \frac{1}{g_{max}} \left( \frac{1}{\epsilon} + \epsilon -\sin^2\theta \right) \frac{SF(q)}{Z}

Where :math:`SF(q)` is the scatter factor, :math:`Z` is the atomic number for the material and :math:`q` is the momentum transfer function:

.. math::
    q = E_0 \sin\left( \frac{\theta}{2}\right) \frac{1}{hc}

In DXMClib the scatter factor for composite materials is obtained by the independent atom approximation, simply put the scatter factor is a weighted average over the atoms in the material. A lookup table for scatter factors are generated for materials in each simulation run and involves computing of a square root thus is more computationally demanding.  


Rayleigh scattering
___________________
Differential cross section for Rayleigh scattering follows Thomson differential cross section for a free electron

.. math::
    \frac{d\rho}{d\Omega} = \frac{r_c^2}{2}\left( 1-\cos^2\theta\right)

This is valid for bound atomic electrons for energies up to 2 keV. For higher energies the photon scatter angle is decreased due to the electronic configuration of the whole atom. The Rayleight differential cross section is like the Thomson cross section but with the introduction of an atomic form factor [#Hubbell]_ :math:`F(q, Z)` where :math:`Z` is the atomic number and :math:`q` is the momentum transfer given by

.. math::
    q = E \sin\left( \frac{\theta}{2}\right) \frac{1}{hc}

for photon energy :math:`E` and :math:`hc` as Planck's constant and speed of light in vacuum. 

The differential cross section for Rayleigh scattering is then

.. math::
    \frac{d\rho}{d\Omega} = \frac{r_c^2}{2}\left( 1-\cos^2\theta\right) \left[F(q, Z)\right]^2

For sampling of scatter angle DXMClib uses a similar approach as the EGS5 monte carlo code. By defining 

.. math::
    A(q_{max}^2) = \int_0^{q_{max}^2} \left[F(q, Z)\right]^2 dq^2

with :math:`q_{max} = E/hc`. :math:`[F(q,Z)]^2/A(q_{max}^2)` can be used as a probability density function and :math:`(1-\cos^2\theta)/2` as a rejection function. To sample a scatter angle :math:`q` is first sampled by :math:`A(q^2) = r_1 A(q_{max}^2)` with :math:`r_1` as a random uniform number in interval [0,1). In DXMClib :math:`q` is found by lookup tables of the integral :math:`A(q^2)`. Scatter angle :math:`\theta` is also obtained from the sampled :math:`q` value. The sampled momentum transfer and therefore scatter angle is accepted if 

.. math::
    \frac{1+\cos^2 \theta}{2} > r_2

where :math:`r_2` is a random number in interval [0, 1). 


Radiation sources
--------------------------
DXMClib models a few radiation sources that should cover most setups in clinical x-ray imaging:

- DX: A x-ray tube source with rectangular collimation.
- CT seq: A CT source for sequental og step and shoot imaging.
- CT spiral: A CT source for spiral aqusitions.
- CT dual: A CT source for dual energy aqusitions.
- Pencil beam: A monochromatic pencil beam.
- Isotropic beam: A rectangular collimated beam, either monochromatic or a user supplied specter.

All of the radiation sources can be positioned arbitrary with the use of source direction cosines and a position vector, although most sources also implements som helper functions to make life easier. Source direction cosines are three orthonormal vectors with the first (:math:`\vec x`) perpendicular to the second vector (:math:`\vec y`) along the anode cathode direction. The third vector :math:`\vec z = \vec x \times \vec y` is along the beam direction. All vectors have basis in the world coordinate system. 

All sources in DXMClib uses the concept of an *exposure* meaning a static position and direction where a number of photon histories are emitted. This makes hardly any sense for conventional x-ray examinations, but for CT examinations an exposure is a position around the patient where a number of photons is emitted.  

.. NOTE::
    Each exposure can run in parallell for computers with multiple cores (all computers nowadays). For optimal performance use atleast twice as many exposures as cores available. Note that number of exposures for CT sources can only be controlled indirectly by setting step angle between each exposure. 

Number of histories per exposure can be set for every source. The optimal number of histories is dependent on the requirered resolution and certainties for a specific application. For example, calculating dose in a large volume of plastics requires much fewer histories compared to a detailed dose map of a CT examination. Voxel size also matters since reaching many events in a small voxel needs more histories. As a guideline, a detailed dose calculation on the voxel level for a thorax examination, either CT or DX, the total number of histories should be about :math:`10^{10}`.   

References
----------
.. [#SUNDERMAN1998] Sunderman E. et al. A Fast Algorithm to Calculate the Exact Radiological Path Through a Pixel Or Voxel Space. Journal of Computing and Information Technology 6(1). December 1998.
.. [#WOODCOCK1965] Woodcock E.R. et al. Techniques used in the GEM code for Monte Carlo neutronics calculations in reactors and other systems of complex geometry. ANL-7050. Argonne National Laboratory, 1965.
.. [#Poludniowski1] Poludniowski G.G. and Evans, P.M. Calculation of x‐ray spectra emerging from an x‐ray tube. Part I. Electron penetration characteristics in x‐ray targets. Med. Phys., 34: 2164-2174 (2007). doi:10.1118/1.2734725
.. [#Poludniowski2] Poludniowski G.G. Calculation of x‐ray spectra emerging from an x‐ray tube. Part II. X‐ray production and filtration in x‐ray targets. Med. Phys., 34: 2175-2186 (2007). doi:10.1118/1.2734726
.. [#Hubbell] Hubbell J.H. et al Atomic form factors, incoherent scattering functions, and photon scattering cross sections, J. Phys. Chem. Ref. Data, Vol.4, No. 3, 1975

