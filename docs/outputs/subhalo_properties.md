# Subhalo properties

There are three types of information datasets that are saved in every HBT-HERONS
output:

- [Bound properties](#bound-properties) based on particles currently bound to the subhalo, like its current mass.
- [Evolutionary milestones](#evolutionary-milestones), such as when the subhalo was initially identified in the simulation.
- [Hierarchical relationships](#hierarchical-relationships) between subhaloes, like which Friends of Friends group a subhalo belongs to.

### Bound properties

These properties are computed for each subhalo within HBT-HERONS using _only_ the
particles that are bound to it at the current output time. As orphan subhaloes do not have
any formally bound particles to them, only a subset of these properties are computed.
We describe what each property means below.

!!! info

    We recommend using the [SOAP](https://github.com/SWIFTSIM/SOAP) post-processing package.
    It is a Python package capable of computing a wide range of subhalo properties using spherical
    and projected apertures, whose sizes are specified based on spherical overdensities or physical
    apertures. It can also use the bound membership information to include or exclude unbound particles.
    It is also unit aware, making each property carry over the correct units.

#### Mass metrics

Different ways to quantify how massive a subhalo is.

| <div style="width:210px">Property</div> | <div style="width:750px">Description</div>                             |
| :------------------------------------- | :-------------------------------------------------------------------------------------------------------------------------------------------------- |
|                `Nbound`                 | The total number of particles bound to the subhalo.                                                                                                 |
|                `Mbound`                 | The total mass of particles bound to the subhalo.                                                                                                   |
|              `NboundType`               | The total number of particles bound to the subhalo, classifed according to their particle type.                                                     |
|              `MboundType`               | The total mass of particles bound to the subhalo, classifed according to their particle type.                                                       |
|             `VmaxPhysical`              | The maximum value of the circularised rotation curve of the subhalo. The centre of the aperture used to compute this quantity corresponds to `ComovingMostBoundPosition`.                                                                            |
|             `BoundM200Crit`             | Mass of a region with a spherical overdensity of 200 times the critical density of the universe. Only bound mass is included and only for centrals. The centre of the aperture used to compute this quantity corresponds to `ComovingMostBoundPosition`. |

!!! warning

    Since the HBT-HERONS definition of spherical overdensity is based on the enclosed mass of **bound particles** only, it does not follow the common convention in the literature, which includes the mass contribution of **all particles**.

    The values provided for $R_{\rm 200}$ and $M_{\rm 200c}$ are therefore lower than the commonly used definition. We only recommend using it as a quick reference, and the user to calculate the spherical overdensity quantities using the accepted definition.

#### Size metrics

Different ways to quantify the size of a subhalo.

| <div style="width:210px">Property</div> | <div style="width:750px">Description</div>                             |
| :-------------------------------------  | :--------------------------------------------------------------------------------------------------------------------------------------------------------- |
|             `RmaxComoving`              | The radius at which the circularised rotation curve of the subhalo reaches its maximum value (`VmaxPhysical`). The centre of the aperture used to compute this quantity corresponds to `ComovingMostBoundPosition`.                                                                             |
|             `RHalfComoving`             | The smallest radius that encloses 50% of the total bound mass. The centre of the aperture used to compute this quantity corresponds to `ComovingMostBoundPosition`.                                                                                                     |
|           `REncloseComoving`            | The smallest radius that encloses 100% of the total bound mass. Useful when interested in doing spatial masking. The centre of the aperture used to compute this quantity corresponds to `ComovingMostBoundPosition`.                                                                |
|         `BoundR200CritComoving`         | The radius of a sphere enclosing a mean density that is 200 times the critical density of the Universe. Only bound mass is included and it is only computed for centrals. The centre of the aperture used to compute this quantity corresponds to `ComovingMostBoundPosition`.  |

#### Position metrics

Different ways to locate the subhalo in 6D phase-space.

| <div style="width:210px">Property</div> | <div style="width:750px">Description</div>                             |
| :-------------------------------------  | :----------------------------------------------------- |
|        `ComovingAveragePosition`        | Mass-weighted average position of all bound particles. |
|        `PhysicalAverageVelocity`        | Mass-weighted average velocity of all bound particles. |
|       `ComovingMostBoundPosition`       | Position of the most bound particle.                   |
|       `PhysicalMostBoundVelocity`       | Velocity of the most bound particle.                   |

#### Shape metrics

Different ways to measure the shape of the subhalo.

| <div style="width:210px">Property</div> | <div style="width:750px">Description</div>                             |
| :-------------------------------------- | :------------------------------------------------ |
| `InertialTensor`                        | Flattened representation of the inertia tensor of the subhalo.                                         |
| `InertialTensorWeighted`                | Flattened representation of the inertia tensor of the subhalo, weighted by particle mass and 3D distance to `ComovingMostBoundPosition`. |

#### Dynamical metrics

Estimates of the internal dynamics and energetics of the subhalo.

| <div style="width:210px">Property</div> | <div style="width:750px">Description</div>                             |
| :-------------------------------------- | :-------------------------------------------------------------------------------------------------------------- |
| `SpecificSelfPotentialEnergy`           | Mass-weighted average potential energy of bound particles.                                                      |
| `SpecificSelfKineticEnergy`             | Mass-weighted average kinetic energy of bound particles in the centre of mass reference frame of the subhalo.   |
| `SpecificAngularMomentum`               | Mass-weighted average angular momentum of bound particles in the centre of mass reference frame of the subhalo. |

#### Miscellaneous

Quantities used for more detailed tracking, restarting runs or debugging.

| <div style="width:210px">Property</div> | <div style="width:750px">Description</div>                             |
| :-------------------------------------- | :----------------------------------------------- |
| `TracerIndex`                           | Bound ranking of the most bound tracer particle. |
| `MostBoundParticleId`                   | ID of the most bound particle.                   |

### Evolutionary milestones

The second type of property corresponds to key events in the evolution of a subhalo, such as when it was first identified or when and if it disrupted. It also includes the peak values of a small subset of properties, which can be useful to have when selecting subhaloes based on their past evolution.

!!! warning

    The milestones get updated as the simulaton is analysed from early to late times. This entails that information from earliest outputs are not reflective of the whole evolution of the subhalo. For example, a subhalo that continously grows in mass will keep updating its `LastMaxMass`. Therefore, the only values that one should use are those stored in the last output of the simulation that you intend to analyse.

| <div style="width:210px">Property</div> | <div style="width:750px">Description</div>                             |
| :-------------------------------------  | :--------------------------------------------------------------------- |
|         `SnapshotOfBirth`          | The output when the subhalo was first identified.                      |
|          `SnapshotOfSink`          | If the subhalo has [sunk](../outputs/merger_trees.md/#subhalo-sinking), the output when that happened. If it has not sunk, it equals -1. |
|         `SnapshotOfDeath`          | If the subhalo has sunk or disrupted, the output when that happened. If neither has happened, it equals -1. |
|     `SnapshotOfLastIsolation`      | If the subhalo has ever been a satellite, the output before it ever became a satellite for the first time. If it has always been a central subhalo, it equals -1.    |
|      `SnapshotOfLastMaxVmax`       | The output when the subhalo reached its maximum value of VmaxPhysical. |
|      `SnapshotOfLastMaxMass`       | The output when the subhalo reached its maximum value of Mbound        |
|              `LastMaxMass`              | The maximum mass that the subhalo has reached so far.                  |
|          `LastMaxVmaxPhysical`          | The maximum VmaxPhysical that the subhalo has reached so far.          |

### Hierarchical relationships

The last type of property is used to connect subhaloes between each other at a fixed time or across time.

| <div style="width:210px">Property</div> | <div style="width:750px">Description</div>                             |
| :-------------------------------------  | :--------------------------------------------------------------------- |
|                `TrackId`                | Unique identifier for the subhalo, which persists across time.                      |
|              `SinkTrackId`              | If the subhalo has sunk, the `TrackId` of the subhalo that accreted it. If it has not happened, it equals -1.              |
|           `DescendantTrackId`           | If the subhalo has sunk or disrupted, the `TrackId` of the subhalo that accreted its most bound particles. If neither has happened, it equals -1. |
|          `NestedParentTrackId`          | The `TrackId` of the parent subhalo. Central subhaloes have -1.       |
|              `HostHaloId`               | The Friends of Friends group that this subhalo is a part of. |
|                 `Rank`                  | The mass ranking of the subhalo compared to all of the subhaloes that have the same `HostHaloId`. |
|                 `Depth`                 | The number of hierarchical connections that the subhalo is away from the central, e.g. 0 for centrals, 1 for satellites, 2 for satellites of satellites.                 |

### Analysis cost statistics

Information used to reconstruct the cost to analyse individual subhaloes and how well parallelised their analysis is. Only saved if `HBT_MEASURE_UNBINDING_TIME` is enabled.

| <div style="width:210px">Property</div> | <div style="width:750px">Description</div>                             |
| :-------------------------------------- | :----------------------------------------------- |
| `MPIRank`                           | MPI rank which analysed the subhalo. |
| `NumberUnbindingIterations`                           | How many unbinding iterations were done before the subhalo disrupted or satisfied the convergence threshold. |

<h5>Timestamps</h5>

All timestamps are provided in nanoseconds since the beginning of the analysis time of the current output. All time measurements occur within [`Subhalo_t::Unbind`](https://github.com/SWIFTSIM/HBT-HERONS/blob/7fbc6515c2e0d4fba7c7f8e427b396f3d5a6f1d5/src/subhalo_unbind.cpp#L463).

| <div style="width:210px">Property</div> | <div style="width:750px">Description</div>                             |
| :-------------------------------------- | :----------------------------------------------- |
| `StartSubhalo`                           | Time when `Subhalo_t::Unbind` was entered. |
| `EndSubhalo`                           | Time when `Subhalo_t::Unbind` was exited. |
| `StartUnbinding`                           | Time when the first unbinding iteration started. It equals `-1` if no unbinding iterations were done, e.g. for pre-existing orphans. |
| `EndUnbinding`                           | Time when the last unbinding iteration finished. It equals `-1` if no unbinding iterations were done, e.g. for pre-existing orphans. |
| `StartCentreRefinement`                           | Time when the subhalo centre refinement started. It equals `-1` if no refinement is done, e.g. for any subhaloes with $N_{\rm bound} \leq N_{\rm subsample}$ or for runs with no subsampling (`MaxSampleSizeOfPotentialEstimate = -1`). |
| `EndCentreRefinement`                           |  Time when the subhalo centre refinement finished. It equals `-1` if no refinement is done, e.g. for any subhaloes with $N_{\rm bound} \leq N_{\rm subsample}$ or for runs with no subsampling (`MaxSampleSizeOfPotentialEstimate = -1`). |
| `StartPhaseSpace`                   | Time when the calculation of the phase-space location of the subhalo started. It equals `-1` if no phase-space calculations are done, e.g. for any orphan subhaloes. |
| `EndPhaseSpace`                   | Time when the calculation of the phase-space location of the subhalo fnished. It equals `-1` if no phase-space calculations are done, e.g. for any orphan subhaloes. |