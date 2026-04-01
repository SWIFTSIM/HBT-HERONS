# Host haloes

The first step that HBT-HERONS does when analysing a simulation output is to assign to every pre-existing subhalo a host halo. This is the first step in the algorithm because:

*   Haloes that contain no subhaloes are used as seeds for new subhaloes in the simulation.
*   The subhalo hierarchy, which determines mass transfer between subhaloes, is based on subhaloes within the same host halo.

## Assignment of host haloes

For a subhalo that was resolved in the previous simulation output, HBT-HERONS uses its `NumTracerHostFinding` most bound tracer particles as identified from that output. It then checks which FoF group each of those particles belongs to in the current simulation output, assigning a particle weight that reflects the bound ranking it had in the last output ($r_{i}$):
$$
w_{i} = \dfrac{1}{1 + \sqrt{r_{i}}}.
$$
The weighting gives greater importance to particles that had larger binding energies when determining the host halo of the subhalo.

HBT-HERONS scores all potential host halo candidates using the total sum of the weights of tracer particles associated to each halo candidate. The assigned host halo of a subgroup is the highest-scoring halo. Note that HBT-HERONS also allows subhaloes to have no assigned host halo, in which case subhaloes are termed [*hostless subhaloes*](#hostless-subhaloes).

<h4>Orphan subhaloes</h4>

The method is the same as for resolved subhaloes, but the difference is that it only uses the most bound tracer particle when the subhalo was last resolved (`SnapshotOfDeath - 1`).

## Hostless subhaloes

There is a population of subhaloes that do not have an assigned host halo in HBT-HERONS. These cases occur when a significant fraction of the most bound tracer particles used to identify host haloes does not belong to a FoF group. This is often caused by the fragmentation of poorly-resolved FoF groups, because a missing particle link can make the number of particles drop down below the minimum number required to identify a FoF group.

However, the lack of a FoF group does not mean that there is no self-bound subhalo in that region. Whereas other subhalo finders do not consider this possibility (their starting point is often a 3D FoF group), HBT-HERONS uses past subhalo memberships to check if there is still a (hostless) subhalo.
