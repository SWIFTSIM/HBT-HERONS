# Subhalo merger trees

The history-based approach of HBT-HERONS means that the evolution of subhaloes is tracked simultaneously with their identification in each simulation output. No additional algorithm is therefore required to find subhalo progenitors or descendants, because the information required to build merger trees is already contained within the catalogues.

## Main evolutionary branch

For many applications, the only information required to follow the evolution of a subhalo is its **main evolutionary branch**. The main evolutionary branch of a subhalo can be identified by selecting the entries across all simulation outputs that share the same value of `TrackId`. Note that the main evolutionary branch will also provide the evolution of the orphan subhalo, which can be removed by discarding all outputs greater than `SnapshotOfDeath`.

We provide an example of how this works in practice [in the following page](../outputs/merger_trees.md#evolution-of-a-single-subhalo).

## Secondary evolutionary branch

If one is interested in the evolution of all of the subhaloes that contributed to the build-up of a given subhalo, then its **secondary evolutionary branches** also need to be considered. Identifying secondary evolutionary branches requires connecting disjoint main evolutionary branches, identified by their `TrackId`, at the time when their associated subhaloes first became orphans. HBT-HERONS identifies the descendants of subhaloes that have just become orphans in two different ways, depending on whether [sinking](./subhalo_sinking.md) or [disruption](./unbinding.md) lead to its conversion to an orphan subhalo.

We provide an example of how to use the information that HBT-HERONS outputs to find secondary evolutionary branches [in the following page](../outputs/merger_trees.md#identifying-subhalo-mergers).

### Disruption descendants

At the beginning of the analysis of each simulation output, HBT-HERONS stores for each resolved subhalo the particle IDs of the `NumTracersForDescendants` most bound tracer particles from the last analysed output. If a subhalo becomes an orphan, HBT-HERONS finds which self-bound subhaloes this set of particles is now bound to.

The descendant subhalo, stored in `DescendantTrackId`, is identified as the `TrackId` that contains the largest share of the tagged tracer particles of the now-orphan subhalo. Note that this may result in `DescendantTrackId = -1` if the largest share of particles are unbound.

The subhaloes where this descendant entry should be used can be identified by `SnapshotOfDeath != SnapshotOfSink == -1`, as well as `SnapshotOfSink > SnapshotOfDeath != -1` (see [unresolved sinking](#unresolved-sinking)).

### Sinking descendants

Subhaloes that are found to overlap in phase-space with the core of another resolved subhalo of `TrackId` store this value as `SinkTrackId`. Contrary to `DescendantTrackId`, the value of `SinkTrackId` can never be `-1`, because sinking needs another existing subhalo to serve as a reference for the phase-space overlap. This means that subhaloes that have sunk can be selected via `SinkTrackId != -1` or `SnapshotOfSink != -1`.

The subhaloes where this descendant entry should be used can be identified by `SnapshotOfDeath = SnapshotOfSink != -1`. Note that HBT-HERONS also computes a `DescendantTrackId` at this time because the subhalo becomes an orphan. The `DescendantTrackId` is the same as `SinkTrackId` in $\approx 99.9\%$ of sinking events, with the discrepant values explained in [unusual descendants](#unusual-descendants).

## Important considerations

Although HBT-HERONS provides robust merger trees and a clean way to navigate them, there will always be more complicated and non-trivial scenarios to bear in mind when analysing a large cosmological simulation.

### Unresolved sinking

Beyond disruption and sinking, subhaloes may also become orphans through a third pathway that is a combination of those two cases: **unresolved sinking**. Unresolved sinking happens when a subhalo first disrupts and is later found to be sunk. These cases can be identified via `SnapshotOfSink > SnapshotOfDeath != -1`. Based on tests using cosmological simulations, around $10\%$ of orphan subhaloes have undergone unresolved sinking, with the fraction decreasing slightly as the resolution is increased.

The "unresolved" qualifier reflects the fact that, if the resolution of the simulation were to be increased, these subhaloes would sink before disrupting. This expectation stems from the fact that dynamical friction in the associated orbital system must be efficient enough for the sinking process to occur, as otherwise it is unlikely for the phase-space overlap to happen. Hence, increasing the resolution will make the subhalo more resilient to disruption before the sinking is complete.

For these instances, we recommend using as the subhalo descendant the value provided in `DescendantTrackId`, and not the one in `SinkTrackId` (they have the same value for $\approx90\%$ of unresolved sinking events).

### Unusual descendants

HBT-HERONS tries to assign `DescendantTrackId` when a subhalo first becomes an orphan, including when the orphan is created because of sinking. There are certain cases where no descendant is found (`DescendantTrackId = -1`) or when the descendant is not the same as the `SinkTrackId`.

Although we make specific recommendations of which entry to use depending on how the orphan was created (`DescendantTrackId` for disruption and unresolved sinking, and `SinkTrackId` for sinking), we provide a table that gives a rough order of magnitude of how common some of these cases are. Note that the subhalo population used to populate the table is entirely made up of orphans, i.e. `SnapshotOfDeath != -1`.

| <div style="width:75px">Description</div> | <div style="width:70">Mask</div> |<div style="width:100px">Statistics</div> |
| :-------------------------------------- | :---- | :-------------------------------------------------------------------------------------------------------------- |
| The subhalo **disrupts** but the majority of its core is not bound to any subhalo. | `DescendantTrackId = -1 &` <br>`SinkTrackId = -1`                                                      | $20\%$ of all orphan subhaloes. |
| The subhalo **sinks** but the majority of its core is bound to a different subhalo from the one it sunk to. | `DescendantTrackId != SinkTrackId &` <br> `DescendantTrackId != -1 & SinkTrackId != -1`                                                      | $0.01\%$ of all orphan subhaloes. |
| The subhalo **sinks** but the majority of its core is not bound to any subhalo. | `DescendantTrackId = -1 &`  <br>  `SinkTrackId != -1`                                                      |$0.0001\%$ of all orphan subhaloes.  |

### Re-resolving orphans

Orphan subhaloes can re-appear as resolved subhaloes in HBT-HERONS. This happens if the only subhalo in a FoF group is an orphan subhalo, because it is designated as the central subhalo of the halo. All particles in the FoF are consequently added to the source of the orphan subhalo, and if they are found to be self-bound, the orphan is re-classified as a resolved subhalo.

This means that some subhaloes, typically those close to the resolution limit of the simulation, may re-appear (and disrupt) several times throughout the simulation. When a orphan subhalo becomes a resolved subhalo again, the entries associated to the disappearance of the subhalo (e.g. `SnapshotOfDeath`, `SnapshotOfSink`) are reset to a value of `-1`.

The choice of allowing orphan subhaloes to re-appear under these conditions is a preferable approach to the alternative of spawning a completely new subhalo in the halo. Doing so would lead to the creation of many short-lived subhaloes with disjoint evolutionary branches, even if the underlying overdense region is the same.