# Population-level statistics

We provide code to make basic population-level statistics. Each of the examples we provide contain an additional "bonus" example that illustrates some of the less-trivial concepts that are part of HBT-HERONS.

## Mass functions

The bound mass of all subhaloes at a given simulation output can be directly loaded via `Mbound`, and it is given in the units specified in `Units/MassInMsunh`. If one is analysing a hydrodynamical simulation and is interested in the bound mass subdivided by type (e.g. stellar mass, gas mass), then `MboundType` is the dataset to use. Note that, because HBT-HERONS stores entries for every subhalo that has ever existed in the simulation, a non-negligible fraction of subhaloes will be orphans (zero bound mass).

In any case, it is easy to make the bound mass function at a given redshift of interest. One can also make the equivalent plot for the maximum circular velocity ($V_{\rm max}$), in which case the datasets `VmaxPhysical` and `Units/VelInKmS` should be used.

=== "After running `toolbox/catalogue_cleanup/SortCatalogues.py`"

    ```python
    import h5py
    from glob import glob
    import matplotlib.pyplot as plt

    # Get ordered list to the sorted subhalo catalogues. Create a dictionary to access its paths by output number.
    catalogue_paths = sorted(glob(f"{SORTED_CATALOGUE_BASE_PATH}/OrderedSubSnap_*.hdf5"))
    catalogue_paths = dict([(int(path[-8:-8+3]),path) for path in catalogue_paths])
    max_output_number = list(catalogue_paths)[-1]

    with h5py.File(catalogue_paths[max_output_number]) as catalogue:
        # Load the current bound mass of all subhaloes, in units of Msun/h
        mass_units = catalogue["Units/MassInMsunh"][0]
        Mbound = catalogue["Subhalos/Mbound"][()] * mass_units

    # Make the plot
    fig, ax1 = plt.subplots(1)
    ax1.plot(np.sort(Mbound)[::-1], np.arange(len(Mbound)) + 1, label="All subhaloes")
    ax1.set_yscale("log")
    ax1.set_xscale("log")
    ax1.set_xlabel(r"$M_{\rm bound} \; [\mathrm{M}_{\rm \odot}h^{-1}]$")
    ax1.set_ylabel(r"$N(\geq M_{\rm bound})$")
    ax1.legend()
    plt.show()
    ```

=== "Without running `toolbox/catalogue_cleanup/SortCatalogues.py`"

    ```python
    import sys
    sys.path.append("<HBT-HERONS_PATH>/toolbox")
    from HBTReader import HBTReader

    # The reader will parse the base folder and parameter file to identify number
    # of outputs and if any are missing.
    catalogue = HBTReader("<HBT-HERONS_CATALOGUE_BASE_PATH>")
    ```
??? abstract "Bonus: the mass function of hostless subhaloes"

    Hostless subhaloes are a HBT-HERONS-specific population of subhaloes that, as the name implies, are not associated to any host haloes. These are generally close to the resolution limit of the simulation, and arise because the host halo they were associated to in the past has fragmented. You can select hostless subhaloes with `HostHaloId = -1`.

    === "After running `toolbox/catalogue_cleanup/SortCatalogues.py`"

        ```python
        import h5py
        from glob import glob
        import matplotlib.pyplot as plt

        # Get ordered list to the sorted subhalo catalogues. Create a dictionary to access its paths by output number.
        catalogue_paths = sorted(glob(f"{SORTED_CATALOGUE_BASE_PATH}/OrderedSubSnap_*.hdf5"))
        catalogue_paths = dict([(int(path[-8:-8+3]),path) for path in catalogue_paths])
        max_output_number = list(catalogue_paths)[-1]

        with h5py.File(catalogue_paths[max_output_number]) as catalogue:
            # Load the current bound mass of all subhaloes, in units of Msun/h
            mass_units = catalogue["Units/MassInMsunh"][0]
            Mbound = catalogue["Subhalos/Mbound"][()] * mass_units

            # Identify which subhaloes are hostless and get their bound mass
            hostless_subhaloes_mask = catalogue["Subhalos/HostHaloId"][()] == -1
            Mbound_hostless = Mbound[hostless_subhaloes_mask]

        # Make the plot
        fig, ax1 = plt.subplots(1)
        ax1.plot(np.sort(Mbound)[::-1], np.arange(len(Mbound)) + 1, label="All subhaloes")
        ax1.plot(np.sort(Mbound_hostless)[::-1], np.arange(len(Mbound_hostless)) + 1, label="Hostless subhaloes")
        ax1.set_yscale("log")
        ax1.set_xscale("log")
        ax1.set_xlabel(r"$M_{\rm bound} \; [\mathrm{M}_{\rm \odot}h^{-1}]$")
        ax1.set_ylabel(r"$N(\geq M_{\rm bound})$")
        ax1.legend()
        plt.show()
        ```

    === "Without running `toolbox/catalogue_cleanup/SortCatalogues.py`"

        ```python
        import sys
        sys.path.append("<HBT-HERONS_PATH>/toolbox")
        from HBTReader import HBTReader

        # The reader will parse the base folder and parameter file to identify number
        # of outputs and if any are missing.
        catalogue = HBTReader("<HBT-HERONS_CATALOGUE_BASE_PATH>")
        ```

## Peak mass function

HBT-HERONS stores information pertaining to past evolution of subhaloes, like the maximum total bound mass it ever reached in its evolution in `LastMaxMass`. Contrary to the [instantaneous bound mass function](#mass-functions), we recommend using the latest output available for a simulation to load `LastMaxMass`, because subhaloes that keep on growing in mass only reach their peak mass at the last output of the simulation.

One may also make the peak maximum velocity function, by using `LastMaxVmaxPhysical` and `Units/VelInKmS`.


=== "After running `toolbox/catalogue_cleanup/SortCatalogues.py`"

    ```python
    import h5py
    from glob import glob
    import matplotlib.pyplot as plt

    # Get ordered list to the sorted subhalo catalogues. Create a dictionary to access its paths by output number.
    catalogue_paths = sorted(glob(f"{SORTED_CATALOGUE_BASE_PATH}/OrderedSubSnap_*.hdf5"))
    catalogue_paths = dict([(int(path[-8:-8+3]),path) for path in catalogue_paths])
    max_output_number = list(catalogue_paths)[-1]

    with h5py.File(catalogue_paths[max_output_number]) as catalogue:
        # Load the maximum-ever bound mass of all subhaloes, in units of Msun/h
        mass_units = catalogue["Units/MassInMsunh"][0]
        Mpeak = catalogue["Subhalos/LastMaxMass"][()] * mass_units

    # Make the plot
    fig, ax1 = plt.subplots(1)
    ax1.plot(np.sort(Mpeak)[::-1], np.arange(len(Mpeak)) + 1, label="All subhaloes")
    ax1.set_yscale("log")
    ax1.set_xscale("log")
    ax1.set_xlabel(r"$M_{\rm peak} \; [\mathrm{M}_{\rm \odot}h^{-1}]$")
    ax1.set_ylabel(r"$N(\geq M_{\rm peak})$")
    ax1.legend()
    plt.show()
    ```

=== "Without running `toolbox/catalogue_cleanup/SortCatalogues.py`"

    ```python
    import sys
    sys.path.append("<HBT-HERONS_PATH>/toolbox")
    from HBTReader import HBTReader

    # The reader will parse the base folder and parameter file to identify number
    # of outputs and if any are missing.
    catalogue = HBTReader("<HBT-HERONS_CATALOGUE_BASE_PATH>")
    ```

??? abstract "Bonus: the peak mass function of orphan subhaloes"

    A benefit of using the peak bound mass, instead of the instantaneous bound mass, is that we can study the population of orphan subhaloes and how massive they were in the past. Orphan subhaloes can be selected by `Nbound = 0` or `SnapshotOfDeath != -1`, since both are equivalent.

    === "After running `toolbox/catalogue_cleanup/SortCatalogues.py`"

        ```python
        import h5py
        from glob import glob
        import matplotlib.pyplot as plt

        # Get ordered list to the sorted subhalo catalogues. Create a dictionary to access its paths by output number.
        catalogue_paths = sorted(glob(f"{SORTED_CATALOGUE_BASE_PATH}/OrderedSubSnap_*.hdf5"))
        catalogue_paths = dict([(int(path[-8:-8+3]),path) for path in catalogue_paths])
        max_output_number = list(catalogue_paths)[-1]

        with h5py.File(catalogue_paths[max_output_number]) as catalogue:
            # Load the maximum-ever bound mass of all subhaloes, in units of Msun/h
            mass_units = catalogue["Units/MassInMsunh"][0]
            Mpeak = catalogue["Subhalos/LastMaxMass"][()] * mass_units

            # Identify curren orphan subhaloes and get their peak bound mass.
            orphan_subhaloes_mask = catalogue["Subhalos/SnapshotOfDeath"][()] != -1
            Mpeak_orphans = Mpeak[orphan_subhaloes_mask]

        # Make the plot
        fig, ax1 = plt.subplots(1)
        ax1.plot(np.sort(Mpeak)[::-1], np.arange(len(Mpeak)) + 1, label="All subhaloes")
        ax1.plot(np.sort(Mpeak_orphans)[::-1], np.arange(len(Mpeak_orphans)) + 1, label="Orphan subhaloes")
        ax1.set_yscale("log")
        ax1.set_xscale("log")
        ax1.set_xlabel(r"$M_{\rm peak} \; [\mathrm{M}_{\rm \odot}h^{-1}]$")
        ax1.set_ylabel(r"$N(\geq M_{\rm peak})$")
        ax1.legend()
        plt.show()
        ```

    === "Without running `toolbox/catalogue_cleanup/SortCatalogues.py`"

        ```python
        import sys
        sys.path.append("<HBT-HERONS_PATH>/toolbox")
        from HBTReader import HBTReader

        # The reader will parse the base folder and parameter file to identify number
        # of outputs and if any are missing.
        catalogue = HBTReader("<HBT-HERONS_CATALOGUE_BASE_PATH>")
        ```

## Merger rates

Another set of evolutionary milestones that HBT-HERONS saves relates to the time when subhaloes became orphans: `SnapshotOfDeath`. Loading this dataset allows the user to measure the rate at which subhaloes become orphans as a function of time, which is related to the merger rate.

=== "After running `toolbox/catalogue_cleanup/SortCatalogues.py`"

    ```python
    import h5py
    from glob import glob
    import matplotlib.pyplot as plt

    # Get ordered list to the sorted subhalo catalogues. Create a dictionary to access its paths by output number.
    catalogue_paths = sorted(glob(f"{SORTED_CATALOGUE_BASE_PATH}/OrderedSubSnap_*.hdf5"))
    catalogue_paths = dict([(int(path[-8:-8+3]),path) for path in catalogue_paths])
    max_output_number = list(catalogue_paths)[-1]

    with h5py.File(catalogue_paths[max_output_number]) as catalogue:
        # Load when subhaloes became orphans, if they have done so.
        SnapshotOfDeath = catalogue["Subhalos/SnapshotOfDeath"][()]

        # Only keep orphan subhaloes
        orphan_subhaloes_mask = SnapshotOfDeath != -1
        SnapshotOfDeath = SnapshotOfDeath[orphan_subhaloes_mask]

    # Make the plot
    fig, ax1 = plt.subplots(1)
    output_number, number_orphans_created = np.unique(SnapshotOfDeath, return_counts = 1)
    ax1.plot(output_number, number_orphans_created, label="All orphans")
    ax1.set_yscale("log")
    ax1.set_xlabel(r"Output number")
    ax1.set_ylabel(r"$N_{\rm new \; orphans}$")
    ax1.legend()
    plt.show()
    ```

=== "Without running `toolbox/catalogue_cleanup/SortCatalogues.py`"

    ```python
    import sys
    sys.path.append("<HBT-HERONS_PATH>/toolbox")
    from HBTReader import HBTReader

    # The reader will parse the base folder and parameter file to identify number
    # of outputs and if any are missing.
    catalogue = HBTReader("<HBT-HERONS_CATALOGUE_BASE_PATH>")

    ```
??? abstract "Bonus: the merger rate of disrupted and sunken subhaloes"

    Subhaloes become orphans in HBT-HERONS via [two possible mechanisms](../algorithm/subhalo_merger_trees.md#secondary-evolutionary-branch): disruption or sinking. We can subdivide the population of orphan subhaloes according to whether they disrupted or sunk by applying additional masks on both `SnapshotOfDeath` and `SnapshotOfSink`.

    === "After running `toolbox/catalogue_cleanup/SortCatalogues.py`"

        ```python
        import h5py
        from glob import glob
        import matplotlib.pyplot as plt

        # Get ordered list to the sorted subhalo catalogues. Create a dictionary to access its paths by output number.
        catalogue_paths = sorted(glob(f"{SORTED_CATALOGUE_BASE_PATH}/OrderedSubSnap_*.hdf5"))
        catalogue_paths = dict([(int(path[-8:-8+3]),path) for path in catalogue_paths])
        max_output_number = list(catalogue_paths)[-1]

        with h5py.File(catalogue_paths[max_output_number]) as catalogue:
            # Load when subhaloes became orphans and when they sunk, if they have done so.
            SnapshotOfDeath = catalogue["Subhalos/SnapshotOfDeath"][()]
            SnapshotOfSink  = catalogue["Subhalos/SnapshotOfSink"][()]

            # Only keep orphan subhaloes, but subdivide population into orphans that disrupted and those that sunk.
            sunk_orphan_subhaloes_mask = (SnapshotOfDeath != -1) & (SnapshotOfSink == SnapshotOfDeath)
            disrupted_orphan_subhaloes_mask = (SnapshotOfDeath != -1) & (SnapshotOfSink == -1)

            SnapshotOfDeath_sunk_subhaloes      = SnapshotOfDeath[sunk_orphan_subhaloes_mask]
            SnapshotOfDeath_disrupted_subhaloes = SnapshotOfDeath[disrupted_orphan_subhaloes_mask]

        # Make the plot
        fig, ax1 = plt.subplots(1)
        output_number, number_orphans_created = np.unique(SnapshotOfDeath_sunk_subhaloes, return_counts = 1)
        ax1.plot(output_number, number_orphans_created, label="Sunk orphans")
        output_number, number_orphans_created = np.unique(SnapshotOfDeath_disrupted_subhaloes, return_counts = 1)
        ax1.plot(output_number, number_orphans_created, label="Disrupted orphans")
        ax1.set_yscale("log")
        ax1.set_xlabel(r"Output number")
        ax1.set_ylabel(r"$N_{\rm new \; orphans}$")
        ax1.legend()
        plt.show()
        ```

    === "Without running `toolbox/catalogue_cleanup/SortCatalogues.py`"

        ```python
        import sys
        sys.path.append("<HBT-HERONS_PATH>/toolbox")
        from HBTReader import HBTReader

        # The reader will parse the base folder and parameter file to identify number
        # of outputs and if any are missing.
        catalogue = HBTReader("<HBT-HERONS_CATALOGUE_BASE_PATH>")
        ```