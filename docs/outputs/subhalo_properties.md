# Subhaloes

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
|             `VmaxPhysical`              | The maximum value of the circularised rotation curve of the subhalo.                                                                             |
|             `BoundM200Crit`             | Mass of a region with a spherical overdensity of 200 times the critical density of the universe. Only bound mass is included and only for centrals. |

#### Size metrics

Different ways to quantify the size of a subhalo.

| <div style="width:210px">Property</div> | <div style="width:750px">Description</div>                             |
| :-------------------------------------  | :--------------------------------------------------------------------------------------------------------------------------------------------------------- |
|             `RmaxComoving`              | The radius at which the circularised rotation curve of the subhalo reaches its maximum value (`VmaxPhysical`).                                                                             |
|             `RHalfComoving`             | The smallest radius that encloses 50% of the total bound mass.                                                                                                     |
|           `REncloseComoving`            | The smallest radius that encloses 100% of the total bound mass. Useful when interested in doing spatial masking.                                                               |
|         `BoundR200CritComoving`         | The radius of a sphere enclosing a mean density that is 200 times the critical density of the Universe. Only bound mass is included and only for centrals. |

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
| `InertialEigenVector`                   | Only computed if GSL was enabled at compile time. |
| `InertialEigenVectorWeighted`           | Only computed if GSL was enabled at compile time. |
| `InertialTensor`                        | centrals.                                         |
| `InertialTensorWeighted`                | centrals.                                         |

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
|         `SnapshotIndexOfBirth`          | The output when the subhalo was first identified.                      |
|          `SnapshotIndexOfSink`          | If the subhalo has merged, the output when that happened.              |
|         `SnapshotIndexOfDeath`          | If the subhalo has merged or disrupted, the output when that happened  |
|     `SnapshotIndexOfLastIsolation`      | If the subhalo has ever been a satellite, the output before it ever became a satellite for the first time.    |
|      `SnapshotIndexOfLastMaxVmax`       | The output when the subhalo reached its maximum value of VmaxPhysical. |
|      `SnapshotIndexOfLastMaxMass`       | The output when the subhalo reached its maximum value of Mbound        |
|              `LastMaxMass`              | The maximum mass that the subhalo has reached so far.                  |
|          `LastMaxVmaxPhysical`          | The maximum VmaxPhysical that the subhalo has reached so far.          |

### Hierarchical relationships

The last type of property that is saved corresponds ot. This helps analyse which
FoF groups, which subhalo is the central, etc.

| <div style="width:210px">Property</div> | <div style="width:750px">Description</div>                             |
| :-------------------------------------  | :--------------------------------------------------------------------- |
|                `TrackId`                | Unique identifier for the subhalo, which persists across time.                      |
|              `SinkTrackId`              | If the subhalo has merged, the `TrackId` of the subhalo that accreted it.              |
|           `DescendantTrackId`           | If the subhalo has merged or disrupted, the `TrackId` of the subhalo that accreted its most bound particles.  |
|          `NestedParentTrackId`          | The `TrackId` of the parent subhalo.         |
|              `HostHaloId`               | The output when the subhalo reached its maximum value of VmaxPhysical. |
|                 `Rank`                  | The output when the subhalo reached its maximum value of Mbound        |
|                 `Depth`                 | The number of hierarchical connections that the subhalo is away from the central, e.g. 0 for centrals, 1 for satellites, 2 for satellites of satellites.                 |