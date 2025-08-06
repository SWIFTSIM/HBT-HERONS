# Home

Welcome to the documentation for **HBT-HERONS**! HBT-HERONS is a history-based algorithm used to find subhaloes in
cosmological dark-matter-only and hydrodynamical simulations. It is a new
version of the [HBT+](<https://github.com/Kambrian/HBTplus>) subhalo finder, with
additions made to improve the tracking of subhaloes in dark-matter-only and hydrodynamical simulations.
The main changes and methodology of HBT-HERONS are described in [Forouhar Moreno et al. 2025]().

## About HBT-HERONS

History-based subhalo finding makes use of the fact that structure in the
Universe forms hierarchically. In other words, any subhalo must have first formed
in relative isolation. In practice, this means that particles associated to a field
subhalo are tagged, and the memberships are propagated forward in time. Those memberships
are used if the subhalo in question has become a satellite, by grouping all particles that share
a common subhalo origin.

This relatively simple, yet powerful approach presents several advantages relative to
traditional subhalo finding methods.

   * The properties and the existence of subhaloes are more robustly recovered in dense enviroments, as the algorithm does not rely on the often complicated instantaneous particle distribution.
   * Merger trees are self-consistently built using the past evolution of subhaloes, resulting in more robust time tracking of subhaloes ([Chandro-Gomez et al. 2025](https://ui.adsabs.harvard.edu/abs/2025MNRAS.539..776C/abstract)).

!!! warning

    Before deciding to use HBT-HERONS, there are several crucial requirements that determine whether your simulation is suitable to be analysed by HBT-HERONS:

    * There should be a healthy amount of particle outputs ($\gtrsim 64$) sampling a variety of redshifts for the history-based subhalo finding to work.
    * These particle outputs should contain 3D FoF group information, which are used by HBT-HERONS to create new subhaloes.
    * The most massive FoF group in your simulation should be able to fit using the memory of a single node. There is no current implementation to spread the particle load of individual FoF groups over several nodes to circumvent memory issues arising from very large groups. This is something that will be solved in the future.

## About this documentation

The documentation touches upon the following aspects of the code:

* [Algorithm](./algorithm/overview.md): a general description of the philosophy and theory behind HBT-HERONS, highlighting the effect of different free parameters that arise from its code implementation.
* [Installation](./installation/requirements.md): required packages, compile-time options and how to compile.
* [Running](./running/start.md): submitting an analyis, mandatory and optional parameters, and general diagnostics to look out for during the execution of HBT-HERONS.
* [Outputs](./outputs/format.md): file structure, which subhalo properties are saved, how to select particles bound to them, and how to navigate the merger trees.

!!! question

    Is the documentation lacking in a certain area? If so, please let us know. We
    want to ensure that the code is as accessible as possible, so your input is
    greatly appreciated.