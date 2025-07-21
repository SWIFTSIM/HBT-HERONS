# Parameters

HBT-HERONS accepts several parameters, which tells the code about the simulation [format and file structure](#file-io), how subhaloes are [tracked and unbound](#subhaloes), the gravity [tree and softening](#gravity), the [unit system](#units) and whether the box is [cosmological](#simulation-box). There are also [miscellaneous](#miscellaneous) options provided for very specific simulation setups.

Some of these parameters are **mandatory**, which means that the code will not work if no information is provided. Other parameters are optional, either because they are not used depending on the simulation setup, or because they default to predefined value. The default values are specified in [`src/config_parser.h`](https://github.com/SWIFTSIM/HBT-HERONS/blob/master/src/config_parser.h) and we provide them when relevant in the tables below.

!!! tip

    When adding a new simulation format, we always recommend relying on the output metadata to load important parameters e.g. the unit system, box size, gravitational softening etc. This is implemented for `Swift` simulations at the moment.

## I/O

These parameters set the format of the simulation, how the output of the simulation is saved and how many outputs to analyse.

### File paths

| <div style="width:260px">Property</div> | <div style="width:50px">Default</div>       | <div style="width:100px">Description</div>       |
| :-------------------------------------- | :----------------------------------------------- | :----------------------------------------------- |
| `SnapshotPath`                   | - | Path to the files containing the particle information.                   |
| `HaloPath`                           | - | Path to the files containing the particle FoF group memberships. |
| `SnapshotFileBase`                           | - | Base name for the snapshot files. |
| `SnapshotDirBase`                           | - | Base name for the directories where snapshot subfiles are saved, if applicable. If there is no subdirectory,  omit this option.|
| `SubhaloPath`                           | - | Base directory where the HBT-HERON catalogues will be saved. |
| `MinSnapshotIndex`                           | 0 | The minimum output number for the HBT-HERONS analysis. |
| `MaxSnapshotIndex`                           | - | The maximum output number for the HBT-HERONS analysis. This does not equal the number of simulation outputs if only a subset of outputs is being analysed, as indicated by `SnapshotIdList`.  |
| `SnapshotIdList`                           | - | Space-separated list indicating a subset of simulation outputs to analyse, if applicable. For example, if we only want to analyse every second snapshot of a ten snapshot run: 1 3 5 7 9. |

### File I/O

!!! warning

    File formats that are not based on Swift outputs have not been explicitly tested using HBT-HERONS. Compatibility with
    alternative formats are there based on what was implemented at the time of forking HBT+. We do not expect the changes that HBT-HERONS makes will affect I/O significantly, but tread carefully.

These options specify the format of the simulation particle and FoF group data. It also determines whether additional information for the particles bound to a subhalo is saved.

| <div style="width:260px">Property</div> | <div style="width:50px">Default</div>       | <div style="width:100px">Description</div>       |
| :-------------------------------------- | :----------------------------------------------- | :----------------------------------------------- |
| `SnapshotFormat`                   | `gadget3` | The format used by files containing the particle information.                   |
| `GroupFileFormat`                           | `gadget3_int` | The format used by the files containing the FoF memberships of particles. |
| `ParticleNullGroupId`                           | - | The value in the FoF particle membership catalogue that corresponds to belonging to no FoF group. |
| `ParticlesSplit`                           | 0 | If particles split in the simulation, in which case their splitting history needs to be provided. Only implemented for `Swift` hydrodynamical  simulations|
| `SaveBoundParticleProperties`                           | 0 | If the position, velocity, mass and type of particles bound to subhaloes should be saved alongside the subhalo catalogues. It also includes thermal energy if thermal unbinding is enabled. |
| `SaveBoundParticleBindingEnergies`                           | 0 | If the internally calculated binding energies of particles bound to subhaloes should be saved alongside the subhalo catalogues. |

## Units

The unit system used by the simulation and used in the outputs. For `Swift` outputs,
the values are set using the snapshot metadata.

| <div style="width:260px">Property</div> | <div style="width:50px">Default</div>       | <div style="width:100px">Description</div>       |
| :-------------------------------------- | :----------------------------------------------- | :----------------------------------------------- |
| `MassInMsunh`                   | 1e10 | Mass of particles in units of Msun / h.                   |
| `LengthInMpch`                           | 1 | Lengths in units of Mpc / h. |
| `VelInKmS`                           | 1| Velocity in units of km  / s. |

## Simulation box

If the box is periodic, and if so, the box size.

<!-- | <div style="width:260px">Property</div> | <div style="width:787px">Description</div>       | -->
| <div style="width:260px">Property</div> | <div style="width:50px">Default</div>       | <div style="width:100px">Description</div>       |
| :-------------------------------------- | :----------------------------------------------- | :----------------------------------------------- |
| `PeriodicBoundaryOn`                    | `1` | Whether the box is periodic.                   |
| `BoxSize`                               | - | Side length of each side of the box. Only cubes are supported. |

## Gravity

### Tree

These options specify how the tree is built and the geometric condition used to open a tree node.

| <div style="width:260px">Property</div> | <div style="width:50px">Default</div>       | <div style="width:100px">Description</div>       |
| :-------------------------------------- | :----------------------------------------------- | :----------------------------------------------- |
| `TreeNodeOpenAngle`                   | 0.45 | Geometric criterion used to determine whether a tree node should be opened or can be used as is. Smaller values make the code slower but more accurate.             |
| `TreeAllocFactor`                           | 0.8 | The maximum number of cells used in the gravity tree, relative to the number of particles. |
| `TreeMinNumOfCells`                           | 10 | Minimum number of tree cells to use. |

### Softening

The value of the softening used for the gravity kernel. Should always be provided, unless the I/O can load it directly from the particle output metadata (e.g. `Swift`).

| <div style="width:260px">Property</div> | <div style="width:50px">Default</div>       | <div style="width:100px">Description</div>       |
| :-------------------------------------- | :----------------------------------------------- | :----------------------------------------------- |
| `SofteningHalo`                   | - | The comoving gravitational softening value. Assumed to be the same for all particle types.                   |
| `MaxPhysicalSofteningHalo`                           | - | The maximum physical gravitational softening value. If no value is provided, the gravitational softening length will always be based on `SofteningHalo`. Assumed to be the same for all particle types. |

## Subhaloes

Various options that tell the code how to do the unbinding, the subhalo self-boundness conditions and the tracking of subhaloes.

### Unbinding

| <div style="width:260px">Property</div> | <div style="width:50px">Default</div>       | <div style="width:100px">Description</div>       |
| :-------------------------------------- | :----------------------------------------------- | :----------------------------------------------- |
| `BoundMassPrecision`                                | 0.995 | Convergence threshold of the fractional difference in the number of bound particles between two consecutive unbinding iterations.                   |
| `SourceSubRelaxFactor`                              | 3 | How many particles a subhalo can have associated to it, relative to its number of bound particles. |
| `MaxSampleSizeOfPotentialEstimate`                  | 1000 | Maximum number of particles used to estimate the gravitational potential of particles. If the subhalo has more particles than this value, a randomly selected `MaxSampleSizeOfPotentialEstimate` particles are used as gravity sources. To disable subsampling, set its value equal to 0.|
| `RefineMostBoundParticle`                           | 1 | If the self-binding energy of the most bound subset of particles should be computed after unbinding the subhalo. This step is done to better identify the most bound particle if potential subsampling was enabled. |
| `BoundFractionCenterRefinement`                     | 0.1 | Fraction of the most bound particles whose self-binding energies are computed to better estimate the most bound particle. The actual value is `max(MaxSampleSizeOfPotentialEstimate, BoundFractionCenterRefinement)` |

### Tracking

Parameters relating to the tracking and merging of subhaloes, as well as the criteria used to determine if a subhalo is self-bound or not.

| <div style="width:260px">Property</div> | <div style="width:50px">Default</div>       | <div style="width:100px">Description</div>       |
| :-------------------------------------- | :----------------------------------------------- | :----------------------------------------------- |
| `TracerParticleTypes`                   | 1 4 | Particle types that can be used as subhalo tracers, given as a space-separated list of particle type numbers. We recommend using time-persistent, collisionless particles (e.g. DM & stars).                   |
| `MinNumPartOfSub`                       | 20 | Minimum number of bound particles required for a subhalo to be resolved, regardless of their particle type.        |
| `MinNumTracerPartOfSub`                  | 10 |   Minimum number of bound particles required for a subhalo to be resolved, of the chosen tracer particle type. |
| `NumTracerHostFinding`                  | 10 |   How many tracer particles are used to identify the host FoF group of subhaloes. The weighting of each particle reflects their binding energy ordering in the previous output. |
| `NumTracersForDescendants`              | 10 |  How many tracer particles are used to identify which subhalo has accreted the core of subhaloes that have disrupted in the current snapshot. |
| `MergeTrappedSubhalos`                  | 1 |  Whether to manually merge self-bound subhaloes that overlap in phase-space. Recommended to always leave on. |
| `MajorProgenitorMassRatio`              | 0.8 |  The threshold used to identify which subhaloes are central candidates within a Friends of Friends group. Expressed relative to the previous bound mass of the most massive subhalo in said host FoF group. |

## Miscellaneous

| <div style="width:260px">Property</div> | <div style="width:50px">Default</div>       | <div style="width:100px">Description</div>       |
| :-------------------------------------- | :----------------------------------------------- | :----------------------------------------------- |
| `SnapshotHasIdBlock`                   | 1 | Whether the simulation outputs have a dataset containing the particle IDs.|
| `ParticleIdNeedHash`                           | 1 | Whether to create a hash map to retrieve particle properties given an ID. |
| `GroupParticleIdMask`                           | - | Hexadecimal mask indicating which digits of the particle IDs are significant. For example, if only the first 4 digits are required (i.e. "1111" in binary), then this would correspond to "F". Only required for peculiar Gadget formats. |
| `SnapshotIdUnsigned`                           | 0 | If the particle IDs are saved as an unsigned integer (4 bytes). |
