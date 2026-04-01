# Subhalo sinking

Subhaloes grow in mass either by accreting diffuse mass (e.g. by stripping mass from less massive subhaloes) or by merging (e.g. by coalescing in phase-space with neighbouring subhaloes). HBT-HERONS identifies which of these two mass growth modes a subhalo has experienced throughout its evolution. This page describes the algorithm used to determine whether two subhaloes have coalesced in phase-space, which is termed as **subhalo sinking**.

The systems where subhalo sinking happens are those in which dynamical friction is non-neglibible, i.e. during major mergers, because one of the subhaloes needs to lose sufficient orbital kinetic energy to coalesce with the neighbouring subhalo. In the frame of reference of the most massive subhalo, it appears as if the lower mass subhalo "sinks" towards its centre, hence the name "subhalo sinking".

In this page, we define how the phase-space location of subhaloes is defined, the metric used to identify whether subhaloes overlap in phase-space, and what happens if so.

## Core position

To determine if subhaloes overlap with each other, we need to assign each of them a phase-space position. There are different choices available, like using the position and velocity of the most bound particle or of the centre of mass of all bound particles. Using a single particle is a noisy estimate, but using the centre of mass of all bound particles may result in a systematically offset position for disturbed subhaloes (as one would expect for those involved in a sinking event).

HBT-HERONS thus defines the phase-space location of a subhalo using the centre of mass position and velocity of the $N_{\rm core}$ most bound particles:

$$
\vec{x}_{\rm core} = \dfrac{\sum \vec{x}_{i}m_{i}}{\sum m_{i}} \; \;\; \;\mathrm{and}  \; \;\; \; \vec{v}_{\rm core} = \dfrac{\sum \vec{v}_{i}m_{i}}{\sum m_{i}}.
$$

The idea is that this estimate represents the phase-space position of the *subhalo core*. However, since we want to evaluate whether two subhaloes overlap in phase-space, we need to associate an "uncertainty" to the measured position and velocity of the subhalo core. Is it the gravitational softening length? A multiple of the half-mass radius? HBT-HERONS uses the position and velocity dispersion of the $N_{\rm core}$ most bound particles as "uncertainty", i.e.

$$
\sigma_{x} = \sqrt{ \dfrac{\sum m_{i}|\vec{x}_{i} - \vec{x}_{\rm core} |^{2}}{\sum m_{i}}} \; \;\; \;\mathrm{and}  \; \;\; \; \sigma_{v} = \sqrt{ \dfrac{\sum m_{i}|\vec{v}_{i} - \vec{v}_{\rm core} |^{2}}{\sum m_{i}}} .
$$

The number of most bound particles that contribute to the sums defined above is determined by a combination of the user-defined values of $N^{\rm min}_{\rm core}$ ([`SubCoreSizeMin`](../running/parameter_file.md#sinking)) and $f_{\rm core}$ ([`SubCoreSizeFactor`](../running/parameter_file.md#sinking)):

$$
N_{\rm core} = \max(N^{\rm min}_{\rm core}, f_{\rm core}N_{\rm bound}).
$$

For the reasons provided at the beginning of this section, making $N_{\rm core} = 1$ or $N_{\rm core} \approx N_{\rm bound}$ is discouraged. HBT-HERONS uses a default of $N^{\rm min}_{\rm core} = 20$ and $f_{\rm core} = 0$, but see [here](./parameter_choices.md/#mincoresize) to see the effect of changing $N^{\rm min}_{\rm core}$.

## Flagging overlaps

Two subhaloes ($i$, $j$) with masses $m_{i} \geq m_{j}$ overlap if their centres are separated in phase-space by less than the phase-space dispersion of their cores:

$$
\Delta \equiv \min(\Delta_{ij},\Delta_{ji}) \leq 2,
$$

where:

$$
\Delta_{ij} = \dfrac{|\vec{x}_{\mathrm{core}, i} - \vec{x}_{\mathrm{core},j}|}{\sigma_{x,i}} + \dfrac{|\vec{v}_{\mathrm{core},i} - \vec{v}_{\mathrm{core},j}|}{\sigma_{v,i}} \; \;\; \;\mathrm{and}  \; \;\; \; \Delta_{ji} = \dfrac{|\vec{x}_{\mathrm{core}, i} - \vec{x}_{\mathrm{core},j}|}{\sigma_{x,j}} + \dfrac{|\vec{v}_{\mathrm{core},i} - \vec{v}_{\mathrm{core},j}|}{\sigma_{v,j}}.
$$

The difference between $\Delta_{ij}$ and $\Delta_{ji}$ is that the first expression uses the core dispersion of the most massive subhalo of the pair, whereas the second expression uses the core dispersion of the least massive subhalo of the pair. One may expect the more massive subhalo to have a larger phase-space dispersion that the lower mass one, e.g. the gravitational potential is deeper and the particle velocities are therefore larger. This is indeed typically the case, but there are a small number of cases in which the lower mass subhalo has a larger dispersion, like when it is close to the resolution limit of the simulation. To not miss these cases and delay the merging of overlapping subhaloes, we choose the smallest phase-space offset.

## Merging subhaloes

HBT-HERONS checks whether subhaloes overlap in phase-space inmediately after subjecting a subhalo to the iterative unbinding procedure. It does so by computing the phase-space offset between the current subhalo and all of the subhaloes contained deeper within its hierarchy tree. This means that subhalo sinking can only happen between subhaloes that fell in as a group, since hierarchical connections need to exist between them.

Any subhalo that satisfies $\Delta \leq 2$ is flagged as sunk (`SnapshotOfSink != -1`), and all of its bound particles (if any) are acquired by the most massive subhalo of the pair. To account for the potential accretion of new mass from the sunk subhalo, HBT-HERONS triggers a re-unbinding of the subhalo with which the subhalo sunk to.

!!! warning

    Poorly resolved systems can also undergo unresolved sinking events. This happens if the subhalo that is sinking towards the centre disrupts before the sinking is complete. These "unresolved" sinking events can be identified by `SnapshotOfSink > SnapshotOfDeath != -1`.