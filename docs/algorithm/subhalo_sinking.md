# Subhalo sinking

In the paradigm of hierarchical structure formation, subhaloes grow in mass either by accreting diffuse mass (e.g. by stripping mass from satellite subhaloes) or by merging with neighbouring subhaloes. Merging in HBT-HERONS precisely refers to the coalesence in phase-space between two or more *resolved* subhaloes. This definition thus excludes mass growth by subhaloes that are stripped below the resolution level of the simulation.

The subhalo systems where this definition of coalescence applies are those in which dynamical friction is non-neglibible, e.g. during major mergers. This is because at least one of the subhaloes needs to lose orbital kinetic energy so that its core coalesces with the core of the neighbouring subhalo. In the frame of reference of the most massive subhalo, this appears as if the lower mass subhalo "sinks" towards its centre, hence the name "subhalo sinking".

In this page, we define how the phase-space location of subhaloes is defined, the metric used to identify whether subhaloes overlap in phase-space, and what happens if so. We further explain why it is an essential step within HBT-HERONS, and how its omission can lead to double-counting subhaloes.

## Subhalo cores

To determine if subhaloes overlap with each other, we first need to locate them in phase-space. There are alternative ways of doing so, like using the position and velocity of the most bound particle or of the centre of mass of all bound particles. For the purpose of identifying subhalo cores, using a single particle is noisy. Using the centre of mass of all bound particles may be offset from the true centre, particularly for disturbed subhaloes (as one would expect for those involved in a major merger).

HBT-HERONS defines the phase-space location of a subhalo core using the centre of mass position and velocity of the $N_{\rm core}$ most bound particles:

$$
\vec{x}_{\rm core} = \dfrac{\sum \vec{x}_{i}m_{i}}{\sum m_{i}} \; \;\; \;\mathrm{and}  \; \;\; \; \vec{v}_{\rm core} = \dfrac{\sum \vec{v}_{i}m_{i}}{\sum m_{i}}.
$$

However, locating the subhalo core is one part of the problem. Since we want to evaluate whether two subhaloes overlap in phase-space, we need to associate an "uncertainty" to the measured position and velocity of the subhalo. Is it the gravitational softening length? A multiple of the half-mass radius? For the "uncertainty" of the core phase-space location, HBT-HERONS uses the position and velocity dispersion of the $N_{\rm core}$ most bound particles:

$$
\sigma_{x} = \sqrt{ \dfrac{\sum m_{i}|\vec{x}_{i} - \vec{x}_{\rm core} |^{2}}{\sum m_{i}}} \; \;\; \;\mathrm{and}  \; \;\; \; \sigma_{v} = \sqrt{ \dfrac{\sum m_{i}|\vec{v}_{i} - \vec{v}_{\rm core} |^{2}}{\sum m_{i}}} .
$$

The number of most bound particles that contribute to the sums defined above is determined by a combination of the user-defined values of $N^{\rm min}_{\rm core}$ ([`SubCoreSizeMin`](../running/parameter_file.md#sinking)) and $f_{\rm core}$ ([`SubCoreSizeFactor`](../running/parameter_file.md#sinking)):

$$
N_{\rm core} = \max(N^{\rm min}_{\rm core}, f_{\rm core}N_{\rm bound}).
$$

For the reasons provided at the beginning of this section, making $N_{\rm core} = 1$ or $N_{\rm core} \approx N_{\rm bound}$ is discouraged. Using small values yields noisy estimates of the core location and its dispersion, whereas large values overestimate its dispersion. By default, HBT-HERONS uses $N^{\rm min}_{\rm core} = 20$ and $f_{\rm core} = 0$.

## Overlap metric

Two subhaloes ($i$, $j$) with masses $m_{i} \geq m_{j}$ overlap if their centres are separated in phase-space by less than the phase-space dispersion of their cores:

$$
\Delta \equiv \min(\Delta_{ij},\Delta_{ji}) \leq 2,
$$

where:

$$
\Delta_{ij} = \dfrac{|\vec{x}_{\mathrm{core}, i} - \vec{x}_{\mathrm{core},j}|}{\sigma_{x,i}} + \dfrac{|\vec{v}_{\mathrm{core},i} - \vec{v}_{\mathrm{core},j}|}{\sigma_{v,i}} \; \;\; \;\mathrm{and}  \; \;\; \; \Delta_{ji} = \dfrac{|\vec{x}_{\mathrm{core}, i} - \vec{x}_{\mathrm{core},j}|}{\sigma_{x,j}} + \dfrac{|\vec{v}_{\mathrm{core},i} - \vec{v}_{\mathrm{core},j}|}{\sigma_{v,j}}.
$$

The difference between $\Delta_{ij}$ and $\Delta_{ji}$ is that the first expression uses the core dispersion of the most massive subhalo of the pair, whereas the second expression uses the core dispersion of the least massive subhalo of the pair. One may expect the more massive subhalo to have a larger core dispersion that the lower mass one, e.g. the gravitational potential is deeper and the particle velocities are therefore larger. This is indeed typically the case, but there are a small number of cases in which the lower mass subhalo has a larger dispersion. These cases typically occur close to the resolution limit of the simulation.
