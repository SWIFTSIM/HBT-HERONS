# Overview

!!! warning

    The algorithm section is currently under construction, and will contain detailed information and visualisations about the internal HBT-HERONS algorithm. Planned sections beyond those already existing:

    * Diffuse mass accretion
    * Mass stripping

Each of the key steps, data and concepts used within the HBT-HERONS algorithm are explained in the following pages:

* [Halo catalogue](./friends_of_friends.md)
* [Host haloes](./host_haloes.md)
* [Subhalo hierarchy](./subhalo_hierarchy.md)
* [Subhalo unbinding](./unbinding.md)
* [Subhalo sinking](./subhalo_sinking.md)
* [Subhalo merger trees](./subhalo_merger_trees.md)

Several free parameters are involved in one or more of the above steps. The default values the code uses are chosen to produce sensible results. Nonetheless provide a [series of tests](./parameter_choices.md) of how changing the subhalo tracking and unbinding parameters affect subhalo-related statistics.

In addition to this documentation, there is also a [wiki](https://github.com/Kambrian/HBTplus/wiki) available for HBT+. Since HBT+ is the version previous to HBT-HERONS, you will find within the HBT+ wiki alternative explanations of how the algorithm works. However, note that some of the content is no longer relevant to HBT-HERONS (e.g. the presence of duplicate particles or the recommended compilation method).