# Parameters

HBT-HERONS accepts several parameters, which tells the code about the simulation [format and file structure](#file-io), how subhaloes are [tracked and unbound](#subhaloes), the gravity [tree and softening](#gravity), the [unit system](#units) and whether the box is [cosmological](#simulation-box). There are also [miscellaneous](#miscellaneous) options provided for very specific simulation setups.

Some of these parameters are **mandatory**, which means that the code will not work if no information is provided. Other parameters are optional either because they are not used depending on the simulation setup, or because they default to given value.

!!! tip

    When adding a new simulation format, we always recommend relying on the output metadata to load important parameters e.g. the unit system, box size, gravitational softening etc. This is implemented for `Swift` simulations at the moment.

## I/O

### File paths

| <div style="width:220px">Property</div> | <div style="width:827px">Description</div>       |
| :-------------------------------------- | :----------------------------------------------- |
| `SnapshotPath`                   | Path to the files containing the particle information.                   |
| `HaloPath`                           | Path to the files containing information used to reconstruct the FoF group membership of particles. |
| `SnapshotFileBase`                           | Base name for the snapshot files, e.g. `<SnapshotFileBase>_<snap_nr>`. |
| `SnapshotDirBase`                           | Base name for the directories where snapshot subfiles are saved, if applicable. If there is only one file per snapshot, omit this option.|
| `SubhaloPath`                           | Base directory of where the HBT catalogues will be saved. |
| `MinSnapshotIndex`                           | The minimum output number of the simulation being analysed. |
| `MaxSnapshotIndex`                           | The maximum output number of the simulation being analysed. |
| `SnapshotIdList`                           | Which subset of snapshots to analyse, if applicable. For example, if we only  want to analyse every second snapshot of a ten snapshot run: 1 3 5 7 9.. |

### File I/O

| <div style="width:260px">Property</div> | <div style="width:787px">Description</div>       |
| :-------------------------------------- | :----------------------------------------------- |
| `SnapshotFormat`                   | The format used by files containing the particle information.                   |
| `GroupFileFormat`                           | The format used by the files containing the FoF memberships of particles. |
| `ParticleNullGroupId`                           |  The value in the FoF particle membership catalogue that corresponds to belonging to no FoF group. |
| `ParticlesSplit`                           | If particles split in the simulation, in which case their splitting history needs to be provided. Only currently implemented for SWIFT hydrodynamical  simulations|
| `SaveBoundParticleProperties`                           | If the position, velocity, mass and type of particles bound to subhaloes should be saved alongside the subhalo catalogues. It also includes thermal energy if thermal unbinding is enabled. |
| `SaveBoundParticleBindingEnergies`                           |  If the internally calculated binding energies of particles bound to subhaloes should be saved alongside the subhalo catalogues. |
<!-- | `MaxSnapshotIndex`                           | The maximum output number of the simulation being analysed. | -->
<!-- | `SnapshotIdList`                           | Which subset of snapshots to analyse, if applicable. For example, if we only  want to analyse every second snapshot of a ten snapshot run: 1 3 5 7 9.. | -->

 0


 0

## Units
 
The unit system used by the simulation and used in the outputs. For `Swift` outputs,
they are set using the snapshot metadata.

| <div style="width:260px">Property</div> | <div style="width:787px">Description</div>       |
| :-------------------------------------- | :----------------------------------------------- |
| `MassInMsunh`                   | Mass of particles in units of Msun / h.                   |
| `LengthInMpch`                           | Lengths in units of Mpc / h. |
| `VelInKmS`                           | Velocity in units of km  / s. |


## Simulation box

| <div style="width:260px">Property</div> | <div style="width:787px">Description</div>       |
| :-------------------------------------- | :----------------------------------------------- |
| `PeriodicBoundaryOn`                   | Whether the box is periodic.                   |
| `BoxSize`                           | Side length of each side of the box. Only cubes are supported. |
<!-- # .
 1

# Side length of each side of the box. Note that only cubes are currently
# supported.
 -1 -->

## Gravity

### Tree

| <div style="width:260px">Property</div> | <div style="width:787px">Description</div>       |
| :-------------------------------------- | :----------------------------------------------- |
| `TreeNodeOpenAngle`                   | Geometric criteria to determine whether a tree node should be opened or can be used as is.                   |
| `TreeAllocFactor`                           | The maximum number of cells used in the gravity tree, relative to the number of particles. |
| `TreeMinNumOfCells`                           | Minimum number of tree cells to use. |


### Softening

| <div style="width:260px">Property</div> | <div style="width:787px">Description</div>       |
| :-------------------------------------- | :----------------------------------------------- |
| `SofteningHalo`                   | The comoving gravitational softening value. Assumed to be the same for all particle types.                   |
| `MaxPhysicalSofteningHalo`                           | The maximum physical gravitational softening value. Assumed to be the same for all particle types. |



<!-- # The comoving softening length of dark matter particles.
 -1

# The maximum physical softening length of dark matter particles.
 -1 -->


<!-- # Geometric criteria to determine whether a tree node should be opened or can be
# used as is.
 0.45

# The maximum number of cells used in the gravity tree, relative to the number
# of particles.
TreeAllocFactor 0.8

# Minimum number of tree cells to use.
 10 -->

## Subhaloes

### Unbinding

| <div style="width:260px">Property</div> | <div style="width:787px">Description</div>       |
| :-------------------------------------- | :----------------------------------------------- |
| `BoundMassPrecision`                   | Convergence threshold of the fractional difference in the number of bound particles between unbinding iterations.                   |
| `SourceSubRelaxFactor`                           | How many particles a subhalo can have associated to its source subhalo, relative to its number of bound particles. |
| `MaxSampleSizeOfPotentialEstimate`                           | Maximum number of particles used to estimate the gravitational potential of particles. If the subhalo has more than this value, a randomly selected `MaxSampleSizeOfPotentialEstimate` particles are used as gravity sources. |
| `RefineMostBoundParticle`                           | If the self-binding energy of the most bound subset of particles should be computed after unbinding the subhalo. This step is done to better identify the most bound particle if potential subsampling was enabled. |
| `BoundFractionCenterRefinement`                           | Fraction of the most bound particles whose self-binding energies are computed to better estimate the most bound particle. The actual value is `max(MaxSampleSizeOfPotentialEstimate, BoundFractionCenterRefinement)` |


<!-- # 
 0.995

# 
 3

# 
 1000

#  Not enabling this
# could make HBT incorrectly identify which particle is the most bound.
 1

# . -->
<!-- BoundFractionCenterRefinement 0.1 -->


### Tracking

| <div style="width:260px">Property</div> | <div style="width:787px">Description</div>       |
| :-------------------------------------- | :----------------------------------------------- |
| `TracerParticleTypes`                   | Particle types that can be used as subhalo tracers. Tracers are used to identify which FoF group hosts a given subhalo, the offsets in phase-space between subhaloes, and the position/velocity of orphan (disrupted) subhaloes. We recommend using time-persistent, collisionless particles (e.g. DM & stars).                   |
| `MinNumPartOfSub`                           | Minimum number of bound particles required for a subhalo to be resolved, regardless of their particle type. |
| `NumTracerHostFinding`                           |  How many tracer particles are used to identify the host FoF group of subhaloes. The weighting of each particle reflects their binding energy ordering in the previous output. |
| `NumTracersForDescendants`                           | How many tracer particles are used to identify which subhalo has accreted the core of subhaloes that have disrupted in the current snapshot. |
| `MergeTrappedSubhalos`                           | Whether to manually merge resolved subhaloes that overlap in phase-space. Recommended to always leave on. |
| `MajorProgenitorMassRatio`                           | The threshold used to identify which subhaloes are central candidates within a FoF. Expressed relative to the most massive subhalo in said host FoF group. |

<!-- # 
 1 4

# 
 20

# 
 10

#
 10

# 
 10

# 
 1

# 
 0.8

# Keeping here as a reference as it is defined, but it is not used internally at
# the moment. It will be good to add as an additional option when computing the
# phase-space distance between overlapping subhaloes. The current behaviour is
# equivalent to SubCoreSizeMin = 10 and SubCoreSizeFactor=0.
# SubCoreSizeMin 20
# SubCoreSizeFactor 0.25 -->

## Miscellaneous

| <div style="width:260px">Property</div> | <div style="width:787px">Description</div>       |
| :-------------------------------------- | :----------------------------------------------- |
| `SnapshotHasIdBlock`                   | Whether particles are already sorted by Particle IDs, and hence no ParticleID dataset is specified within the snapshot file.                  |
| `ParticleIdNeedHash`                           | Whether to create a hash map to retrieve particle properties given an ID type. |
| `GroupParticleIdMask`                           |  Hexadecimal mask indicating which digits of the particle IDs are significant. For example, if only the first 4 digits are required (i.e. "1111" in binary), then this would correspond to "F". Only required for peculiar Gadget formats. |
| `SnapshotIdUnsigned`                           | If the particle IDs are saved as an unsigned integer. |

<!-- 
# 
 1

# 
 1

# 


# 
 -->
