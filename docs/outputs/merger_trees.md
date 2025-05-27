# Time evolution

The history-based approach of HBT-HERONS means that the evolution of subhaloes is tracked
simultaneously with their identification in each simulation output. No additional algorithm is therefore required to link subhaloes forward in time, so the information required to build merger trees is already contained within the catalogues.

Since HBT-HERONS analyses a simulation from early to late times, the information used to follow the evolution of subhaloes is primarily represented through the notion of descendant subhaloes. There are two types of fate a given subhalo can have at any given output, and the dataset used to find the descendant of a subhalo varies depending on which of the two it experienced:

* **It remains self-bound**. The subhalo with the same `TrackId`, which is unique and time-persistent for each subhalo, is its descendant.
* **It merges with another subhalo**. Depending on the way in which this happened, i.e. through disruption or because their cores overlapped in phase-space, the `DescendantTrackId` or `SinkTrackId` values identify which subhalo it merged with.

We provide examples on how to use these datasets to follow the [complete evolution of a subhalo of interest](#evolution-of-a-single-subhalo), and to [identify subhaloes that merged with it](#identifying-subhalo-mergers).

## Evolution of a single subhalo

The `TrackId` of a subhalo is a unique identifier that persists in time throughout the simulation. Thus, following the evolution of a subhalo from the first output when it is found in the simulation (`SnapshotIndexOfBirth`) until the last output when it is resolved as self-bound (`SnapshotIndexOfDeath`) only requires knowing its `TrackId`. In fact, the subhalo can still be tracked as an orphan subhalo after its "death", but only a subset of properties are computed in that case. 

<h4>Code example</h4>

Here we show how to follow the mass evolution of the most massive subhalo identified in the last output of a simulation.
We assume the simulation has 64 outputs throughout this example.

=== "After running `toolbox/SortCatalogues.py`"

    ```python
    import h5py

    # Path to where the sorted catalogues are located.
    catalogue_path = "<SORTED_CATALOGUE_BASE_PATH>/OrderedSubSnap_{output_number:03d}.hdf5"

    # Maximum output number for this example (64 outputs but SnapshotIndex uses 0-indexing)
    max_output_number = 63

    # Get the TrackId of the most massive subhalo, when it was first identified and when it disrupted/merged.
    with h5py.File(catalogue_path.format(output_number = max_output_number)) as catalogue:
        TrackId_to_follow = catalogue['Subhalos']['Mbound'][()].argmax()
        output_start = catalogue['Subhalos']['SnapshotIndexOfBirth'][TrackId_to_follow]
        output_end   = catalogue['Subhalos']['SnapshotIndexOfDeath'][TrackId_to_follow]

    # If output_end is equal to -1, that means it is still resolved at the time when the output was saved. 
    output_end = output_end if output_end != -1 else max_output_number

    # Create an array to hold values we are interested in tracking (number of bound particles)
    Nbound_evolution = - np.ones(output_end - output_start + 1)

    # Iterate over catalogues to obtain Nbound value of the entry with the TrackId we want to follow. 
    for i, output_number in enumerate(range(output_start, output_end + 1)):
        with h5py.File(catalogue_path.format(output_number = output_number)) as catalogue:
          Nbound_evolution[i] = catalogue['Subhalos/Nbound'][TrackId_to_follow]
    ```

=== "Without running `toolbox/SortCatalogues.py`"


    ``` python
    EXAMPLE INCOMING SOON!
    ```

We now have an array that contains the number of bound particles for that subhalo across time. We can plot it using
the code below:

```bash
import matplotlib.pyplot as plt

fig, ax1 = plt.subplots(1)
ax1.plot(np.arange(output_start, output_end + 1), Nbound_evolution, 'k-')
ax1.set_xlabel('Output Number')
ax1.set_ylabel('Total number of bound particles')
ax1.set_yscale('log')
plt.show()
```

## Identifying subhalo mergers

The catalogues also provide sufficient infomation to identify which subhaloes merged together, helping to establish links between different evolutionary branches through so-called secondary progenitors. Obtaining the secondary progenitors of a given subhalo is more involved that following its main branch. One complication is that subhaloes disappear from the HBT-HERONS catalogues in two different ways:

* [Subhalo disruption](#subhalo-disruption) as a consequence of subhaloes no longer being self-bound. 
* [Subhalo sinking](#subhalo-sinking) when the subhalo is still self-bound but its core coalesces in phase-space with the core of another subhalo. 

HBT-HERONS provides different information depending on which of the two processes led to the removal of a subhalo from the simulation. Commonly used merger trees do not provide this two-category classification, as following the entire process of 
subhalo sinking is difficult and is inadecuately followed in traditional subhalo finders. 

### Subhalo disruption

Disruption occurs when the subhalo is no longer considered to be self-bound. There are two conditions within HBT-HERONS used to determine whether a subhalo is self-bound:

* The total number of bound particles has to be equal or greater than `MinNumPartOfSub`.
* The total number of bound tracer particles has to be equal or greater than `MinNumTracerPartOfSub`.

If either of these conditions is not satisfied, the subhalo is considered as disrupted (`Nbound = 0`) and is 
subsequently tracked as an orphan subhalo. 

To identify the descendant of a disrupted subhalo, HBT-HERONS uses the `NumTracersForDescendants` most bound tracer particles when the subhalo was last self-bound. The descendant subhalo (`DescendantTrackId`) is the one that is self-bound and contains the majority of the aforementioned particles.

<h4>Code example</h4>

Here we show how to identify all subhaloes that disrupted and hence merged with the same subhalo of the [first example](#evolution-of-a-single-subhalo). We assume the simulation has 64 outputs throughout this example.

!!! note

    At the moment, `DescendantTrackId` always equals -1 for every subhalo, except for the output when a subhalo disrupts. This means that we need to iterate over every output to which subhaloes merged with another subhalo by disrupting. We are likely to change this in the future so that the information is available after the disruption takes place, to remove the requirement to load every output. 

=== "After running `toolbox/SortCatalogues.py`"

    ```python
    import h5py

    # Path to where the sorted catalogues are located.
    catalogue_path = "<SORTED_CATALOGUE_BASE_PATH>/OrderedSubSnap_{output_number:03d}.hdf5"

    # Maximum output number for this example (64 outputs but SnapshotIndex uses 0-indexing)
    max_output_number = 63

    # Get the TrackId of the most massive subhalo, when it was first identified and when it disrupted/merged.
    with h5py.File(catalogue_path.format(output_number = max_output_number)) as catalogue:
        TrackId_to_follow = catalogue['Subhalos']['Mbound'][()].argmax()
        disrupted_progenitors  = np.where(catalogue["Subhalos/DescendantTrackId"][()] == TrackId_to_follow)[0]
    ```

=== "Without running `toolbox/SortCatalogues.py`"


    ``` python
    EXAMPLE INCOMING SOON!
    ```

The `disrupted_progenitors` array contains a `TrackID` for each subhalo that disrupted and whose core ended bound to the example subhalo (`TrackId_to_follow`). If you want to find *every* subhalo that is connected to the example subhalo through disruption, e.g. to build its full merger tree, you will need to repeat the above code for each of the `disrupted_progenitors` until a complete list is built. The evolution of each subhalo prior to disrupting can be followed in the same way as for an [individual subhalo](#evolution-of-a-single-subhalo).

### Subhalo sinking

The sinking of a subhalo refers to process of its core becoming indistinguishable in phase-space from the core of another subhalo, whilst both remain self-bound. The coalescence occurs as dynamical friction makes the core of the least massive of the pair to "sink" towards the core of more massive one. As the efficiency of dynamical friction increases when the masses of both subhaloes becomes comparable, sinking only happens for non-negligible mass ratios.

Formally, HBT-HERONS checks the phase-space offset between the cores of subhaloes that have a hierarchical connection. If it identifies them to be within a threshold distance, HBT-HERONS flags the least massive of the two as being sunk and removes it from the simulation. The `TrackId` of the most massive of the two is assigned as the `SinkTrackId` of the subhalo that sunk.

<h4>Code example</h4>

Here we show how to identify all subhaloes that directly sunk with the same subhalo of the [first example](#evolution-of-a-single-subhalo). We assume the simulation has 64 outputs throughout this example.

=== "After running `toolbox/SortCatalogues.py`"

    ```python
    import h5py

    # Path to where the sorted catalogues are located.
    catalogue_path = "<SORTED_CATALOGUE_BASE_PATH>/OrderedSubSnap_{output_number:03d}.hdf5"

    # Maximum output number for this example (64 outputs but SnapshotIndex uses 0-indexing)
    max_output_number = 63

    # Get the TrackId of the most massive subhalo, and use it to find the TrackIds that directly sunk to it.
    with h5py.File(catalogue_path.format(output_number = max_output_number)) as catalogue:
        TrackId_to_follow = catalogue['Subhalos']['Mbound'][()].argmax()
        sink_progenitors  = np.where(catalogue["Subhalos/SinkTrackId"][()] == TrackId_to_follow)[0]
    ```

=== "Without running `toolbox/SortCatalogues.py`"


    ``` python
    EXAMPLE INCOMING SOON!
    ```

The `sink_progenitors` array contains a `TrackID` for each subhalo that sunk with the example subhalo (`TrackId_to_follow`). If you want to find *every* subhalo that is connected to the example subhalo through sinking, e.g. to build its full merger tree, you will need to repeat the above code for each of the `sink_progenitors` until a complete list is built. The evolution of each subhalo prior to sinking can be followed in the same way as for an [individual subhalo](#evolution-of-a-single-subhalo).