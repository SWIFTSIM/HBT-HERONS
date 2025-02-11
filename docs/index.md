# Home

Welcome to the documentation for **HBT-HERONS**! HBT-HERONS is a history-based algorithm used to find subhaloes in 
cosmological dark-matter-only and hydrodynamical simulations. It is a new 
version of the [HBT+](<https://github.com/Kambrian/HBTplus>) subhalo finder, with 
additions made to improve the tracking of subhaloes in dark-matter-only and hydrodynamical simulations. 
The main changes are described in [Forouhar Moreno et al. 2025]().

!!! question 

    Is the documentation lacking in a certain area? If so, please let us know. We
    want to ensure that the code is as accessible as possible, so your input is
    greatly appreciated. 

## About

History-based subhalo finding makes use of the fact that structure in the 
Universe forms hierarchically. In other words, any subhalo must have first formed
in relative isolation. In practice, this means that particles associated to a field
subhalo are tagged, and the memberships are propagated forward in time. Those memberships
are used if the subhalo in question has become a satellite, by grouping all particles that share
a common subhalo origin. 

This relatively simple, yet powerful approach presents several advantages relative to 
traditional subhalo finding methods.

   * The properties (and existence) of subhaloes are more robustly recovered in dense enviroments,
     as the algorithm does not rely on the often complicated instantaneous particle distribution. 
   * Merger trees are self-consistently built using the past evolution of subhaloes.
