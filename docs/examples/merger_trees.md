# Merger trees

In this page we explain how to follow the evolution of a subhalo of interest ([main progenitor](#evolution-of-a-subhalo)), and how to identify subhaloes that contributed to its build-up ([secondary progenitors](#secondary-progenitors)).

## Main progenitor

The `TrackId` of a subhalo is a unique identifier that persists in time throughout the simulation. Thus, following the evolution of a subhalo from the first output when it is found in the simulation (`SnapshotOfBirth`) until the last output when it is resolved as self-bound (`SnapshotOfDeath`) only requires knowing its `TrackId`. In fact, the subhalo can still be tracked as an orphan subhalo after its "death", but only a subset of properties are computed in that case.

In the code below we follow the bound mass evolution of the most massive subhalo, identified at the last output of the simulation.

=== "After running `toolbox/catalogue_cleanup/SortCatalogues.py`"

    ```python
    import h5py
    from glob import glob
    import matplotlib.pyplot as plt

    # Get ordered list to the sorted subhalo catalogues. Create a dictionary to access its paths by output number.
    catalogue_paths = sorted(glob(f"{SORTED_CATALOGUE_BASE_PATH}/OrderedSubSnap_*.hdf5"))
    catalogue_paths = dict([(int(path[-8:-8+3]),path) for path in catalogue_paths])
    max_output_number = list(catalogue_paths)[-1]

    # Get the TrackId of the most massive subhalo at the last available output,
    # when it was first identified and when it disrupted/merged.
    with h5py.File(catalogue_paths[max_output_number]) as catalogue:
        TrackId_to_follow = catalogue["Subhalos/Mbound"][()].argmax()
        output_start = catalogue["Subhalos/SnapshotOfBirth"][TrackId_to_follow]
        output_end   = catalogue["Subhalos/SnapshotOfDeath"][TrackId_to_follow]

    # If output_end is equal to -1, that means it is still resolved at the time when the output was saved.
    output_end = output_end if output_end != -1 else max_output_number

    # Create an array to hold values we are interested in tracking.
    Mbound_evolution   = - np.ones(max_output_number)

    # Iterate over catalogues to obtain Nbound value of the entry with the TrackId we want to follow.
    for output_number in range(output_start, output_end):
        with h5py.File(catalogue_paths[output_number]) as catalogue:
            mass_units = catalogue["Units/MassInMsunh"][0]
            Mbound_evolution[output_number] = catalogue["Subhalos/Mbound"][TrackId_to_follow] * mass_units

    fig, ax1 = plt.subplots(1)
    ax1.plot(Mbound_evolution)
    ax1.set_xlabel('Output Number')
    ax1.set_ylabel(r"$M_{\rm bound} \; [\mathrm{M}_{\rm \odot}h^{-1}]$")
    ax1.set_yscale('log')
    plt.show()
    ```

=== "Without running `toolbox/catalogue_cleanup/SortCatalogues.py`"

    ```python
    import sys
    sys.path.append(f"{HBT_HERONS_PATH}/toolbox")
    from HBTReader import HBTReader
    import matplotlib.pyplot as plt

    # The reader will parse the base folder and parameter file to identify number
    # of outputs and if any are missing.
    catalogue = HBTReader(f"{RAW_CATALOGUE_BASE_PATH}")

    # Get the TrackId of the most massive subhalo. The reader loads the latest
    # available snapshot by default and all subhalo properties.
    subhaloes = catalogue.LoadSubhalos()
    TrackId_to_follow = subhaloes["TrackId"][subhaloes["Mbound"].argmax()]

    # Get its bound mass evolution, which returns by default all properties and
    # the associated scale factors and snashot output numbers.
    subhalo_evolution  = catalogue.GetTrackEvolution(TrackId_to_follow)
    Mbound_evolution   = subhalo_evolution["Mbound"] * catalogue.get_mass_units()

    fig, ax1 = plt.subplots(1)
    ax1.plot(Mbound_evolution)
    ax1.set_xlabel('Output Number')
    ax1.set_ylabel(r"$M_{\rm bound} \; [\mathrm{M}_{\rm \odot}h^{-1}]$")
    ax1.set_yscale('log')
    plt.show()
    ```

## Secondary progenitors

The catalogues also provide sufficient infomation to identify which subhaloes merged together, helping to establish links between different evolutionary branches through so-called secondary progenitors. Obtaining the secondary progenitors of a given subhalo is more involved that following its main branch. One complication is that subhaloes disappear from the HBT-HERONS catalogues in two different ways:

* [Subhalo disruption](#disrupted-progenitors) as a consequence of subhaloes no longer being self-bound.
* [Subhalo sinking](#sunk-progenitors) when the subhalo is still self-bound but its core coalesces in phase-space with that of another subhalo.

HBT-HERONS provides different information depending on which of the two processes led to the removal of a subhalo from the simulation. Commonly used merger trees do not provide this two-category classification, as following the entire process of
subhalo sinking is difficult and is inadecuately followed in traditional subhalo finders.

### Disrupted progenitors

Disruption occurs when the subhalo is no longer considered to be self-bound. To identify the descendant of a disrupted subhalo, HBT-HERONS uses the `NumTracersForDescendants` most bound tracer particles when the subhalo was last self-bound. The descendant subhalo (`DescendantTrackId`) is the one that is self-bound and contains the majority of the aforementioned particles.

Here we show how to identify all subhaloes that disrupted and merged with the subhalo used in the [first example](#main-progenitor).

=== "After running `toolbox/catalogue_cleanup/SortCatalogues.py`"

    ```python
    import h5py
    from glob import glob
    import matplotlib.pyplot as plt

    # Get ordered list to the sorted subhalo catalogues. Create a dictionary to access its paths by output number.
    catalogue_paths = sorted(glob(f"{SORTED_CATALOGUE_BASE_PATH}/OrderedSubSnap_*.hdf5"))
    catalogue_paths = dict([(int(path[-8:-8+3]),path) for path in catalogue_paths])
    max_output_number = list(catalogue_paths)[-1]

    # Get the TrackId of the most massive subhalo at the last available output,
    # when it was first identified and when it disrupted/merged.
    with h5py.File(catalogue_paths[max_output_number]) as catalogue:
        TrackId_to_follow = catalogue["Subhalos/Mbound"][()].argmax()

        # Select orphan subhaloes that disrupted or underwent unresolved sinking.
        orphaned_subhalo_mask  = catalogue["Subhalos/SnapshotOfDeath"][()] != -1
        disrupted_subhalo_mask = (catalogue["Subhalos/SnapshotOfSink"][()] == -1) | (catalogue["Subhalos/SnapshotOfDeath"][()] < catalogue["Subhalos/SnapshotOfSink"][()])

        # For these cases, we always use the DescendantTrackId to find the secondary progenitors of the object of interest.
        disrupted_subhalo_progenitors = np.flatnonzero((catalogue["Subhalos/DescendantTrackId"][()] ==  TrackId_to_follow) & orphaned_subhalo_mask & disrupted_subhalo_mask)

    # Plot the evolution of a random subset of progenitors
    Mbound_evolution = - np.ones((len(disrupted_subhalo_progenitors), max_output_number))

    # Iterate over catalogues to obtain bound mass evolution of all disrupted progenitors.
    for output_number in range(max_output_number):
        with h5py.File(catalogue_paths[output_number]) as catalogue:

            # Identify which progenitors exist in the current snapshot
            number_subhaloes = catalogue["NumberOfSubhalosInAllFiles"][0]
            progenitors_in_current_snapshot = disrupted_subhalo_progenitors < number_subhaloes

            # Load the mass of existing progenitors
            mass_units = catalogue["Units/MassInMsunh"][0]
            Mbound_evolution[progenitors_in_current_snapshot, output_number] = catalogue["Subhalos/Mbound"][()][disrupted_subhalo_progenitors[progenitors_in_current_snapshot]] * mass_units

    fig, ax1 = plt.subplots(1)
    ax1.plot(Mbound_evolution.T)
    ax1.set_xlabel('Output Number')
    ax1.set_ylabel(r"$M_{\rm bound} \; [\mathrm{M}_{\rm \odot}h^{-1}]$")
    ax1.set_yscale('log')
    plt.show()
    ```

=== "Without running `toolbox/catalogue_cleanup/SortCatalogues.py`"

    ``` python
    import sys
    sys.path.append("<HBT-HERONS_PATH>/toolbox")
    from HBTReader import HBTReader
    import matplotlib.pyplot as plt

    # The reader will parse the base folder and parameter file to identify number
    # of outputs and if any are missing.
    catalogue = HBTReader(f"{RAW_CATALOGUE_BASE_PATH}")

    # Get the TrackId of the most massive subhalo. The reader loads the latest
    # available snapshot by default and all subhalo properties.
    subhaloes = catalogue.LoadSubhalos()
    TrackId_to_follow = subhaloes["TrackId"][subhaloes["Mbound"].argmax()]

    # Select orphan subhaloes that disrupted or underwent unresolved sinking.
    orphaned_subhalo_mask  = subhaloes["SnapshotOfDeath"] != -1
    disrupted_subhalo_mask = (subhaloes["SnapshotOfSink"][()] == -1) | (subhaloes["SnapshotOfDeath"][()] < subhaloes["SnapshotOfSink"][()])
    disrupted_progenitors  = subhaloes["TrackId"][(subhaloes["DescendantTrackId"] == TrackId_to_follow) & orphaned_subhalo_mask & disrupted_subhalo_mask]

    # Plot the evolution of a random subset of progenitors
    Mbound_evolution = - np.ones((len(disrupted_progenitors),catalogue.MaximumSnapshotIndex + 1))
    for i, disrupted_progenitor in enumerate(disrupted_progenitors):
        subhalo_evolution = catalogue.GetTrackEvolution(disrupted_progenitor)
        Mbound_evolution[i,subhalo_evolution["Snapshot"]] = subhalo_evolution["Mbound"] * catalogue.get_mass_units()

    fig, ax1 = plt.subplots(1)
    ax1.plot(Mbound_evolution.T)
    ax1.set_xlabel('Output Number')
    ax1.set_ylabel(r"$M_{\rm bound} \; [\mathrm{M}_{\rm \odot}h^{-1}]$")
    ax1.set_yscale('log')
    plt.show()
    ```

### Sunk progenitors

The sinking of a subhalo refers to process of its core becoming indistinguishable in phase-space from the core of another subhalo, whilst both remain self-bound. HBT-HERONS identifies a sinking even when two subhaloes are within a threshold distance in phase-space, flagging the least massive of the two as being sunk and converting it to an orphan subhalo. The `TrackId` of the most massive of the two is assigned as the `SinkTrackId` of the subhalo that sunk.

Here we show how to identify all subhaloes that sunk and merged with the subhalo used in the [first example](#main-progenitor).

=== "After running `toolbox/catalogue_cleanup/SortCatalogues.py`"

    ```python
    import h5py
    from glob import glob
    import matplotlib.pyplot as plt

    # Get ordered list to the sorted subhalo catalogues. Create a dictionary to access its paths by output number.
    catalogue_paths = sorted(glob(f"{SORTED_CATALOGUE_BASE_PATH}/OrderedSubSnap_*.hdf5"))
    catalogue_paths = dict([(int(path[-8:-8+3]),path) for path in catalogue_paths])
    max_output_number = list(catalogue_paths)[-1]

    # Get the TrackId of the most massive subhalo at the last available output,
    # when it was first identified and when it disrupted/merged.
    with h5py.File(catalogue_paths[max_output_number]) as catalogue:
        TrackId_to_follow = catalogue["Subhalos/Mbound"][()].argmax()

        # Select orphan subhaloes that disrupted or underwent unresolved sinking.
        orphaned_subhalo_mask  = catalogue["Subhalos/SnapshotOfDeath"][()] != -1
        sunk_subhalo_mask = catalogue["Subhalos/SnapshotOfSink"][()] == catalogue["Subhalos/SnapshotOfDeath"][()]

        # For these cases, we always use the SinkTrackId to find the secondary progenitors of the object of interest.
        sunk_subhalo_progenitors = np.flatnonzero((catalogue["Subhalos/SinkTrackId"][()] ==  TrackId_to_follow) & orphaned_subhalo_mask & sunk_subhalo_mask)

    # Plot the evolution of a random subset of progenitors
    Mbound_evolution = - np.ones((len(sunk_subhalo_progenitors), max_output_number))

    # Iterate over catalogues to obtain bound mass evolution of all disrupted progenitors.
    for output_number in range(max_output_number):
        with h5py.File(catalogue_paths[output_number]) as catalogue:

            # Identify which progenitors exist in the current snapshot
            number_subhaloes = catalogue["NumberOfdisrupted_subhalo_progenitors_in_current_snapshotSubhalosInAllFiles"][0]
            progenitors_in_current_snapshot = sunk_subhalo_progenitors < number_subhaloes

            # Load the mass of existing progenitors
            mass_units = catalogue["Units/MassInMsunh"][0]
            Mbound_evolution[progenitors_in_current_snapshot, output_number] = catalogue["Subhalos/Mbound"][()][sunk_subhalo_progenitors[progenitors_in_current_snapshot]] * mass_units

    fig, ax1 = plt.subplots(1)
    ax1.plot(Mbound_evolution.T)
    ax1.set_xlabel('Output Number')
    ax1.set_ylabel(r"$M_{\rm bound} \; [\mathrm{M}_{\rm \odot}h^{-1}]$")
    ax1.set_yscale('log')
    plt.show()
    ```

=== "Without running `toolbox/catalogue_cleanup/SortCatalogues.py`"

    ``` python
    import sys
    sys.path.append("<HBT-HERONS_PATH>/toolbox")
    from HBTReader import HBTReader
    import matplotlib.pyplot as plt

    # The reader will parse the base folder and parameter file to identify number
    # of outputs and if any are missing.
    catalogue = HBTReader(f"{RAW_CATALOGUE_BASE_PATH}")

    # Get the TrackId of the most massive subhalo. The reader loads the latest
    # available snapshot by default and all subhalo properties.
    subhaloes = catalogue.LoadSubhalos()
    TrackId_to_follow = subhaloes["TrackId"][subhaloes["Mbound"].argmax()]

    # Select orphan subhaloes that underwent resolved sinking
    orphaned_subhalo_mask  = subhaloes["SnapshotOfDeath"] != -1
    sunk_subhalo_mask = (subhaloes["SnapshotOfSink"][()] == subhaloes["SnapshotOfDeath"] )
    sunk_progenitors  = subhaloes["TrackId"][(subhaloes["DescendantTrackId"] == TrackId_to_follow) & orphaned_subhalo_mask & sunk_subhalo_mask]

    # Plot the evolution of a random subset of progenitors
    Mbound_evolution = - np.ones((len(sunk_progenitors),catalogue.MaximumSnapshotIndex + 1))
    for i, sunk_progenitor in enumerate(sunk_progenitors):
        subhalo_evolution = catalogue.GetTrackEvolution(sunk_progenitor)
        Mbound_evolution[i,subhalo_evolution["Snapshot"]] = subhalo_evolution["Mbound"] * catalogue.get_mass_units()

    fig, ax1 = plt.subplots(1)
    ax1.plot(Mbound_evolution.T)
    ax1.set_xlabel('Output Number')
    ax1.set_ylabel(r"$M_{\rm bound} \; [\mathrm{M}_{\rm \odot}h^{-1}]$")
    ax1.set_yscale('log')
    plt.show()
    ```
