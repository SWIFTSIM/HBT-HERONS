# Merger trees

The history-based approach of HBT-HERONS means that self-consistent
merger trees are constructed simultaneously when finding subhaloes. No additional algorithm is
therefore required to link subhaloes across time, as the catalogues already contain all
the information required to traverse the merger trees.

The merger tree of a single subhalo is traditionally visualised as being made up of 
distinct branches. There are two types of branches:

* [Main progenitor branch](#main-progenitor), which is the subhalo structure
* [Secondary progenitor branches](#secondary-progenitors), which are subhaloes that grew independently from the main progenitor but subsequently merged with it. 

We explain below how to use the information saved by HBT-HERONS to follow the evolution of main progenitors, as  
well as to identify their secondary progenitors. Example code is also provided for this purpose.

## Main progenitor

The `TrackId` of a subhalo is a unique identifier that persists in time throughout
the simulation. This means that following the main progenitor branch of a given
subhalo simply relies on obtaining the properties of the entry with the same `TrackId` 
that we are interested in tracking.

<h4>Code example</h4>

In this snipet we show how to follow the mass time evolution for a subhalo that has piked our interest. 
Its `TrackId` equals `23`, which we will use to find the relevant information across different outputs.

=== "After running `toolbox/SortCatalogues.py`"

    ```python
    import h5py

    # Identifier for the subhalo whose evolution we want to follow.
    TrackId_to_follow = 23

    # We get the times at which the subhalo was identified and disrupted/merged.
    with h5py.File(<PATH_TO_CATALOGUE>) as catalogue:
        output_start = catalogue['Subhalos']['SnapshotIndexOfBirth'][TrackId_to_follow]
        output_end   = catalogue['Subhalos']['SnapshotIndexOfDeath'][TrackId_to_follow]

    # If output_end is equal to -1, that means it is still resolved by the output 
    # we just opened.
    output_end = output_end if output_end != -1 else <MAX_SNAPSHOT>

    # Create an array to hold values we are interested in tracking (Mass vs Time)
    mass_evolution     = - np.ones(output_end - output_start)
    redshift_evolution = - np.ones(output_end - output_start)

    # Iterate over catalogues to obtain Mbound value of the entry with the TrackId we want to follow
    for i, output_number in enumerate(range(output_start, output_end + 1)):
        with h5py.File(f"<PATH_TO_CATALOGUE>/SortedOutput_{output_number:03d}") as catalogue:
          mass_evolution    [i] = catalogue['Subhalo/Mbound'][TrackId_to_follow]
          redshift_evolution[i] = catalogue['Subhalo/Mbound'][TrackId_to_follow]
    ```

=== "Without running `toolbox/SortCatalogues.py`"


    ``` python
    # TODO: finish writing.
    ```

## Secondary progenitors

The catalogues also provide sufficient infomation to identify the secondary merger tree branches
of a given subhalo, although it is somewhat more involved that following the main branch. One 
complication is that subhaloes disappear from the HBT-HERONS catalogues in two different ways:

* [Subhalo disruptions](#subhalo-disruption) as a consequence of subhaloes no longer being self-bound. 
* [Subhalo mergers](#subhalo-mergers) when the subhalo is still self-bound but coalesces in phase-space with another subhalo. 

Commonly used merger trees do not provide this two-category classification, as following the entire process of 
subhalo merging is difficult and is inadecuately followed in traditional subhalo finders. HBT-HERONS provides 
different information depending on which process led to the removal of a subhalo from the simulation.

### Subhalo disruption

Disruption occurs when the subhalo is no longer considered to be self-bound. There are two conditions within HBT-HERONS used to determine whether a subhalo is self-bound:

* The total number of bound particles has to be equal or greater than `MinNumPartOfSub`.
* The total number of bound tracer particles has to be equal or greater than `MinNumTracerPartOfSub`.

If either of these conditions is not satisfied, the subhalo is considered as disrupted (`Nbound = 0`) and is 
subsequently tracked as an orphan subhalo. HBT-HERONS uses the `NumTracersForDescendants` most bound tracer particles identified when the subhalo was last self-bound. The descendant subhalo of the disrupted subhalo is the one that is self-bound and contains the majority of particles used. The `TrackId` of the descendant subhalo is saved in `DescendantTrackId`.

<h4>Code example</h4>

In this snipet we identify all of the subhaloes that disrupted and whose most bound
tracer particles end up primarly bound to the same subhalo of the previous example (`TrackId = 23`).

=== "After running `toolbox/SortCatalogues.py`"

    ```python
    import h5py

    # Subhalo whose (disrupted) secondary progenitors we want to identify.
    TrackId_to_follow = 23

    # We get the time when the subhalo was last identified as self-bound.
    with h5py.File(<PATH_TO_CATALOGUE>) as catalogue:
        output_end   = catalogue['Subhalos']['SnapshotIndexOfDeath'][TrackId_to_follow]

    # The disrupted subhalo is done 
    # Open that output, since it contains the most up to date information for its merger tree.
    with h5py.File(f"<PATH_TO_CATALOGUE>/SortedOutput_{output_end:03d}") as catalogue:
        merger_progenitors = np.where(catalogue["Subhalos/SinkTrackId"][()] == TrackId_to_follow)[0]

    # We iterate over all outputs and find subhaloes who are connected to 
    secondary_progenitor_trackids = []
    for i, output_number in enumerate(range(output_start, output_end + 1)):
        with h5py.File(f"<PATH_TO_CATALOGUE>/SortedOutput_{output_number:03d}") as catalogue:
          mass_evolution    [i] = catalogue['Subhalo/Mbound'][TrackId_to_follow]
          redshift_evolution[i] = catalogue['Subhalo/Mbound'][TrackId_to_follow]

    # All the TrackIds that merged directly with TrackId_to_follow are now in merger_progenitors. If you also require
    # knowing which subhaloes merged before with those that merged with TrackId_to_follow, the above process should
    # be recursively done for each TrackId in merger_progenitor.
    ```

=== "Without running `toolbox/SortCatalogues.py`"


    ``` python
    # TODO: finish writing.
    ```

### Subhalo mergers

Mergers differ from disruptions in the sense that they are driven by the effects of dynamical friction 
instead of tidal stripping. As dynamical friction becomes more efficient as the mass ratio of the two objects becomes
more comparable, this process generally happens between subhaloes with comparable mass ratios.

Formally, the phase-space offset is the distance between the centre of mass 
position and velocity of the 10 most bound particles of two subhaloes. The offset is normalised 
by the position and velocity dispersion of the 10 most bound particles of the most massive subhalo of the two.
If the check is triggered, the most massive subhalo of the two accretes all of the particles of the least massive one.
The merger gets registered through the `SinkTrackId`.

<h4>Code example</h4>

In this snipet we identify all of the subhaloes whose cores physically overlapped in phase-space
with the same subhalo of the previous example (`TrackId = 23`). 

=== "After running `toolbox/SortCatalogues.py`"

    ```python
    import h5py

    # Subhalo whose (merged) secondary progenitors we want to identify.
    TrackId_to_follow = 23

    # We get the time when the subhalo was last identified as self-bound.
    with h5py.File(<PATH_TO_CATALOGUE>) as catalogue:
        output_end   = catalogue['Subhalos']['SnapshotIndexOfDeath'][TrackId_to_follow]

    # Open that output, since it contains the most up to date information for its merger tree.
    with h5py.File(f"<PATH_TO_CATALOGUE>/SortedOutput_{output_end:03d}") as catalogue:
        merger_progenitors = np.where(catalogue["Subhalos/SinkTrackId"][()] == TrackId_to_follow)[0]

    # All the TrackIds that merged directly with TrackId_to_follow are now in merger_progenitors. If you also require
    # knowing which subhaloes merged before with those that merged with TrackId_to_follow, the above process should
    # be recursively done for each TrackId in merger_progenitor.
    ```

=== "Without running `toolbox/SortCatalogues.py`"


    ``` python
    # TODO: finish writing.
    ```