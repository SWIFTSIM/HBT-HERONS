/* This macro contains the definition of each user-defined input that HBT-HERONS
 * accepts. The column format is as follows: Type, name, default value. */
#define AVAILABLE_PARAMETERS(PARAMETER)                                        \
  /* Path where the particle data are saved. */                              \
  PARAMETER(std::string, SnapshotPath)                                           \
  /* Path where the information used to assign a FoF to each particle are 
   * saved. */                                                               \
  PARAMETER(std::string, HaloPath)                                           \
  /* Path where the catalogues will be saved. */                             \
  PARAMETER(std::string, SubhaloPath)                                           \
  /* Base string used to generate the names of each particle snapshot. */    \
  PARAMETER(std::string, SnapshotFileBase)                                           \
  /* Maximum physical gravitational softening length (used for all 
   * particles). */                                                          \
  PARAMETER(std::string, SnapshotDirBase, "")                                    \
  /* Format of the particle snapshot, generally dependent on the code used 
   * to generate them. */                                                    \
  PARAMETER(std::string, SnapshotFormat, "swiftsim")                                     \
  /* Format of the FoF catalogues, generally dependent on the code used to
   * generate them. */                                                       \
  PARAMETER(std::string, GroupFileFormat, "swiftsim")                                    \
  /* Total number of snapshots that will be analysed. */                     \
  PARAMETER(int,         MaxSnapshotIndex)                                           \
  /* Box size of the simulation, if periodic boundary conditions are used. */\
  PARAMETER(HBTReal,      BoxSize)                                          \
  /* Comoving gravitational softening length (used for all particles). */    \
  PARAMETER(HBTReal,      SofteningHalo)                                          \
  /* Maximum physical gravitational softening length (used for all 
   * particles). */                                                          \
  PARAMETER(HBTReal,      MaxPhysicalSofteningHalo, -1)                               \
  PARAMETER(int,          MaxConcurrentIO, 10)                                        \
  /* Snapshot to start the analysis from. */                                 \
  PARAMETER(int, MinSnapshotIndex, 0)                                       \
  /* Minimum number of bound particles for a subhalo to be resolved. */      \
  PARAMETER(int, MinNumPartOfSub, 20)                                        \
  /* Minimum number of bound tracer particles for a subhalo to be 
   * resolved. */                                                            \
  PARAMETER(int, MinNumTracerPartOfSub, 10)                                  \
  /* Number of tracer particles used to assign FoF hosts to subhaloes */     \
  PARAMETER(int, NumTracerHostFinding, MinNumTracerPartOfSub)                                   \
  /* Number of tracer particles used to assign descendants to subhaloes */   \
  PARAMETER(int, NumTracersForDescendants,MinNumTracerPartOfSub)                               \
  PARAMETER(long, GroupParticleIdMask, 0)                                    \
  /* Whether the simulation has periodic boundary conditions. */             \
  PARAMETER(bool, PeriodicBoundaryOn, true)                                     \
  PARAMETER(bool, SnapshotHasIdBlock, true)                                     \
  PARAMETER(bool, ParticleIdNeedHash, true)                                                \
  PARAMETER(bool, SnapshotIdUnsigned, false)                                                \
  /* Whether to also save mass, position, velocity (and internal energy) of 
   * all particles bound to subhaloes. */                                    \
  PARAMETER(bool, SaveBoundParticleProperties, false)                                       \
  /* Whether to also  save mass, position, velocity (and internal energy) of 
   * all particles bound to subhaloes. */                                    \
  PARAMETER(bool, SaveBoundParticleBindingEnergies, false)                                  \
  /* Whether to also  save mass, position, velocity (and internal energy) of 
   * all particles bound to subhaloes. */                                    \
  PARAMETER(bool, SaveBoundParticlePotentialEnergies, false)                                \
  /* Whether to manually merge subhaloes that overlap in phase-space */      \
  PARAMETER(bool, MergeTrappedSubhalos, true)                                              \
  /* If snapshots will not be analysed consecutively, this vector tells the 
   * code which numbers will be read. */                                     \
  PARAMETER(vector<int>, SnapshotIdList)                                             \
  /* Particle types used to trace the host FoF of subhaloes (and orphans) */ \
  PARAMETER(vector<int>, TracerParticleTypes,  (std::vector<int>{1, 4}))                                        \
  /* Particle types that cannot be subsampled when unbinding */              \
  PARAMETER(vector<int>, DoNotSubsampleParticleTypes, (std::vector<int>{5}))                                \
  /* If snapshots are named differently within the same run. Otherwise better 
   * to use SnapshotFileBase. */                                             \
  PARAMETER(vector<string>, SnapshotNameList)                                        \
  /* Mass threshold above which subhaloes are candidates to become centrals
   * of their FoF group. Expressed relative to the bound mass of the most
   * massive subhalo. */                                                     \
  PARAMETER(HBTReal, MajorProgenitorMassRatio, 0.8)                                       \
  /* Threshold ratio in the number of bound particles between consecutive
   * unbinding iterations to deem the unbinding as converged. */             \
  PARAMETER(HBTReal, BoundMassPrecision, 0.995)                                             \
  /* Maximum number of particles a subhalo can keep track of, relative to its
   * number of bound particles. */                                           \
  PARAMETER(HBTReal, SourceSubRelaxFactor, 3)                                           \
  PARAMETER(HBTReal, TreeAllocFactor, 0.8)                                                \
  PARAMETER(HBTReal, TreeNodeOpenAngle, 0.45)                                              \
  PARAMETER(HBTInt, TreeMinNumOfCells, 10)                                               \
  /* Value used to identify which particles do not belong to any FoF group */\
  PARAMETER(HBTInt, ParticleNullGroupId)                                             \
  /* Maximum number of particles that can be subsampled used during unbinding.
   * Triggers subsampling if MaxSampleSizeOfPotentialEstimate is less than 
   * the number of particles that can be subsampled. */                      \
  PARAMETER(HBTInt, MaxSampleSizeOfPotentialEstimate, 1000)                                \
  /* Compute the self-binding energy of the most bound subset of particles if     
   * subsampling is enabled, to refine the subhalo coordinates. */            \
  PARAMETER(bool, RefineMostBoundParticle, true)                                           \
  /* Minimum fraction of most bound particles used to refine the subhalo 
   * coordinates */                                                          \
  PARAMETER(float, BoundFractionCenterRefinement, 0.1)                                    \
  /* If baryonic particles can split. Relevant to swift simulations */       \
  PARAMETER(int, ParticlesSplit, -1)                                                     \
  /* Fraction of most bound particles used to estimate position and velocity
   * when doing merger checks. */                                            \
  PARAMETER(HBTReal, SubCoreSizeFactor, 0.25)                                              \
  /* Minimum number of most bound particles used to estimate position and 
   * velocity when doing merger checks. */                                   \
  PARAMETER(HBTInt, SubCoreSizeMin, 20)                                                  \

/* This macro contains the definition of configuration-related parameters that
 * are derived from user-defined values. The column format is as follows: 
 * Type, name, default value. */
#define DERIVED_PARAMETERS(PARAMETER)                                                     \
  /* Used to determine if mandatory parameters have been set. */               \
  PARAMETER(vector<bool>,      IsSet)                                               \
  /* Base units of the simulation */                                           \
  PARAMETER(HBTReal, MassInMsunh, 1e10)                                              \
  PARAMETER(HBTReal, LengthInMpch, 1)                                             \
  PARAMETER(HBTReal, VelInKmS, 1)                                                 \
  /* Bitmask used to identify particle types used to trace the host FoF of       
   * subhaloes (and orphans) */                                                \
  PARAMETER(int, TracerParticleBitMask)                                                \
  /* Bitmask used to identify particle types that cannot be subsampled during 
   * unbinding */                                                              \
  PARAMETER(int, DoNotSubsampleParticleBitMask)                                        \
  /*derived parameters; do not require user input*/                            \
  PARAMETER(HBTReal,TreeNodeOpenAngleSquare)                                           \
  PARAMETER(HBTReal,TreeNodeResolution)                                                \
  PARAMETER(HBTReal,TreeNodeResolutionHalf)                                            \
  PARAMETER(HBTReal,BoxHalf)                                                           \
  PARAMETER(bool, GroupLoadedFullParticle, false)

/* This block of defines handles parameter declarations and their initial 
 * default values. We omit the setting of parameters without 3 elements 
 * (i.e. without specified default column) */
#define DECLARE3(type, name, default_value) type name = default_value;
#define DECLARE2(type, name) type name; /* no default assignment */
#define GET_MACRO(_1, _2, _3, NAME, ...) NAME
#define DECLARE(...) GET_MACRO(__VA_ARGS__, DECLARE3, DECLARE2)(__VA_ARGS__)