# Internal halo population

In this page we showcase how to study the population of subhaloes within a single halo, with particular emphasis to the concept of subhalo hierarchy and subhalo parents.

## Satellite peak mass function

Subhaloes that are in the same halo share the value of `HostHaloId`, allowing us to select satellite systems. Each satellite has an assigned hierarchical depth, which encodes how many connections away a subhalo is from the central. Subhaloes that have <br> `Depth = 0`  are centrals, `Depth= 1` are satellites, `Depth= 2` are satellites-of-satellites, etc

In the provided code, we show how to select satellite subhaloes that are in the same halo as the most massive subhalo at the last output of the simulation. We will use this information to plot the satellite peak mass function, which is further subdivided according to the `Depth` of satellites.

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

        # Load bound masses of all subhaloes in the simulation
        mass_units = catalogue["Units/MassInMsunh"][0]
        Mbound = catalogue["Subhalos/Mbound"][()] * mass_units
        Mpeak  = catalogue["Subhalos/LastMaxMass"][()] * mass_units
        Depth  = catalogue["Subhalos/Depth"][()]

        # Select the most massive subhalo in the box and its host halo
        TargetTrackId = Mbound.argmax()
        HostHaloId = catalogue["Subhalos/HostHaloId"][TargetTrackId]

        # Select all satellite subhaloes that share the same host halo
        satellite_subhaloes = np.flatnonzero((catalogue["Subhalos/HostHaloId"][()] == HostHaloId) & (Depth != 0))

        Mpeak_satellites = Mpeak[satellite_subhaloes]
        depth_satellites = Depth[satellite_subhaloes]

    # Make the plot
    fig, ax1 = plt.subplots(1)

    # All satellites together
    ax1.plot(np.sort(Mpeak_satellites)[::-1], np.arange(len(Mpeak_satellites)) + 1, label="All satellite subhaloes")

    # Do it according to hierarchical depth
    for depth_level in range(1, depth_satellites.max() + 1):
        ax1.plot(np.sort(Mpeak_satellites[depth_satellites == depth_level])[::-1],
                np.arange(len(Mpeak_satellites[depth_satellites == depth_level])) + 1, label=f"Depth = {depth_level:d}  subhaloes")

    ax1.set_yscale("log")
    ax1.set_xscale("log")
    ax1.set_xlabel(r"$M_{\rm peak} \; [\mathrm{M}_{\rm \odot}h^{-1}]$")
    ax1.set_ylabel(r"$N(\geq M_{\rm peak})$")
    ax1.legend()
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

## Hierarchical connections

Every satellite subhalo has `NestedParentTrackId != -1`, which tells us the `TrackId` of the *immediate* parent subhalo of a satellite subhalo. In other words, the `NestedParentTrackId` tells us which `Depth = l` subhalo is the parent of a subhalo of <br> `Depth = l + 1`. Note that the value of `NestedParentTrackId` reflects whether a previously-associated collection of subhaloes were accreted at the same time by a halo, but the subhalo systems may be far apart due to post-infall orbital scattering.

In the following code, we plot the spatial position all subhaloes in the halo that contains the most massive subhalo at the last output of the simulation. We draw connections between subhaloes following the value of `NestedParentTrackId` to indicate connections between subhaloes.

!!! warning

    We have not taken into account periodic boundary conditions when plotting the positions of subhaloes.

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

        # Required unit information
        mass_units = catalogue["Units/MassInMsunh"][0]
        length_units = catalogue["Units/LengthInMpch"][0]

        # Load current and peak bound masses of all subhaloes in the simulation
        Mbound = catalogue["Subhalos/Mbound"][()] * mass_units
        Mpeak  = catalogue["Subhalos/LastMaxMass"][()] * mass_units

        # Load the positions of subhaloes
        Centres = catalogue["Subhalos/ComovingMostBoundPosition"][()] * length_units

        # Information relating to subhalo hierarchy
        Depth  = catalogue["Subhalos/Depth"][()]
        ParentTrackId = catalogue["Subhalos/NestedParentTrackId"][()]

        # Select the most massive subhalo in the box and its host halo
        TargetTrackId = Mbound.argmax()
        HostHaloId = catalogue["Subhalos/HostHaloId"][TargetTrackId]

        # Select all subhaloes that share the same host halo (including central)
        subhaloes_in_halo = np.flatnonzero((catalogue["Subhalos/HostHaloId"][()] == HostHaloId))

        # Information needed about the subhaloes in this halo
        Mpeak_satellites = Mpeak[subhaloes_in_halo]
        Depth_satellites = Depth[subhaloes_in_halo]
        ParentTrackId_satellites = ParentTrackId[subhaloes_in_halo]
        CentreSatellites = Centres[subhaloes_in_halo]

    # Make the plot
    fig, ax1 = plt.subplots(1)

    ax1.scatter(CentreSatellites[:,0],CentreSatellites[:,1],
                lw=0.3, s=0.01 * Mpeak_satellites**(1/3), edgecolor='k',c=Depth_satellites / Depth_satellites.max(),alpha=0.7)

    # Now connect the subhaloes according to hierarchy
    for (current_track_id, parent) in zip(subhaloes_in_halo,ParentTrackId_satellites):
        if parent == -1:
            continue
        ax1.plot([CentreSatellites[subhaloes_in_halo == current_track_id][0][0], CentreSatellites[subhaloes_in_halo == parent][0][0]],
                [CentreSatellites[subhaloes_in_halo == current_track_id][0][1], CentreSatellites[subhaloes_in_halo == parent][0][1]],
                    lw=0.05,color='k')

    ax1.set_xlabel(r"$x \; [\mathrm{Mpc}h^{-1}]$")
    ax1.set_ylabel(r"$y \; [\mathrm{Mpc}h^{-1}]$")

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

    # Required unit information
    mass_units = catalogue.get_mass_units()
    length_units = catalogue.get_length_units()

    # The reader loads the latest available snapshot by default and all subhalo properties.
    subhaloes = catalogue.LoadSubhalos()

    # Load current and peak bound masses of all subhaloes in the simulation
    Mbound = subhaloes["Mbound"] * mass_units
    Mpeak  = subhaloes["LastMaxMass"] * mass_units

    # Load the positions of subhaloes
    Centres = subhaloes["ComovingMostBoundPosition"] * length_units

    # Information relating to subhalo hierarchy
    Depth  = subhaloes["Depth"]
    ParentTrackId = subhaloes["NestedParentTrackId"]

    # Select the most massive subhalo in the box and its host halo
    TargetTrackId = Mbound.argmax()
    HostHaloId = subhaloes["HostHaloId"][TargetTrackId]

    # Select all subhaloes that share the same host halo (including central)
    subhaloes_in_halo = subhaloes["HostHaloId"] == HostHaloId

    # Information needed about the subhaloes in this halo
    TrackId_subhaloes_in_halo = subhaloes["TrackId"][subhaloes_in_halo]
    Mpeak_subhaloes_in_halo = Mpeak[subhaloes_in_halo]
    Depth_subhaloes_in_halo = Depth[subhaloes_in_halo]
    ParentTrackId_subhaloes_in_halo = ParentTrackId[subhaloes_in_halo]
    CentreSatellites = Centres[subhaloes_in_halo]

    # Make the plot
    fig, ax1 = plt.subplots(1)

    ax1.scatter(CentreSatellites[:,0],CentreSatellites[:,1],
                lw=0.3, s=0.01 * Mpeak_subhaloes_in_halo**(1/3), edgecolor='k',c=Depth_subhaloes_in_halo / Depth_subhaloes_in_halo.max(),alpha=0.7)

    # Now connect the subhaloes according to hierarchy
    for (current_track_id, parent) in zip(TrackId_subhaloes_in_halo,ParentTrackId_subhaloes_in_halo):
        if parent == -1:
            continue
        ax1.plot([CentreSatellites[TrackId_subhaloes_in_halo == current_track_id][0][0], CentreSatellites[TrackId_subhaloes_in_halo == parent][0][0]],
                [CentreSatellites[TrackId_subhaloes_in_halo == current_track_id][0][1], CentreSatellites[TrackId_subhaloes_in_halo == parent][0][1]],
                    lw=0.05,color='k')

    ax1.set_xlabel(r"$x \; [\mathrm{Mpc}h^{-1}]$")
    ax1.set_ylabel(r"$y \; [\mathrm{Mpc}h^{-1}]$")

    plt.show()
    ```