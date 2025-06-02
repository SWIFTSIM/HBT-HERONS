/* This macro contains the definition of each user-defined input that HBT-HERONS
 * accepts. The column format is as follows: Type, name, default value. */
#define AVAILABLE_PARAMETERS                                                   \
    /* Path where the particle data are saved. */                              \
    X(std::string, SnapshotPath    )                                           \
    /* Path where the information used to assign a FoF to each particle are 
     * saved. */                                                               \
    X(std::string, HaloPath        )                                           \
    /* Path where the catalogues will be saved. */                             \
    X(std::string, SubhaloPath     )                                           \
    /* Base string used to generate the names of each particle snapshot. */    \
    X(std::string, SnapshotFileBase)                                           \
    /* Maximum physical gravitational softening length (used for all 
     * particles). */                                                          \
    X(std::string,      SnapshotDirBase   )                                    \
    /* Format of the particle snapshot, generally dependent on the code used 
     * to generate them. */                                                    \
    X(std::string,      SnapshotFormat   )                                     \
    /* Format of the FoF catalogues, generally dependent on the code used to
     * generate them. */                                                       \
    X(std::string,      GroupFileFormat   )                                    \
    /* Total number of snapshots that will be analysed. */                     \
    X(int,         MaxSnapshotIndex)                                           \
    /* Box size of the simulation, if periodic boundary conditions are used. */\
    X(HBTReal,      BoxSize         )                                          \
    /* Comoving gravitational softening length (used for all particles). */    \
    X(HBTReal,      SofteningHalo   )                                          \
    /* Maximum physical gravitational softening length (used for all 
     * particles). */                                                          \
    X(HBTReal,      MaxPhysicalSofteningHalo   )                               \
    X(int,          MaxConcurrentIO   )                                        \
    /* Snapshot to start the analysis from. */                                 \
    X(int,          MinSnapshotIndex   )                                       \
    /* Minimum number of bound particles for a subhalo to be resolved. */      \
    X(int,          MinNumPartOfSub   )                                        \
    /* Minimum number of bound tracer particles for a subhalo to be 
     * resolved. */                                                            \
    X(int,          MinNumTracerPartOfSub   )                                  \
    /* Number of tracer particles used to assign FoF hosts to subhaloes */     \
    X(int,          NumTracerHostFinding   )                                   \
    /* Number of tracer particles used to assign descendants to subhaloes */   \
    X(int,          NumTracersForDescendants   )                               \
    X(long,         GroupParticleIdMask   )                                    \
    /* Whether the simulation has periodic boundary conditions. */             \
    X(bool,         PeriodicBoundaryOn   )                                     \
    X(bool,         SnapshotHasIdBlock   )                                     \
    X(bool, ParticleIdNeedHash)                                                \
    X(bool, SnapshotIdUnsigned)                                                \
    /* Whether to also save mass, position, velocity (and internal energy) of 
     * all particles bound to subhaloes. */                                    \
    X(bool, SaveBoundParticleProperties)                                       \
    /* Whether to also  save mass, position, velocity (and internal energy) of 
     * all particles bound to subhaloes. */                                    \
    X(bool, SaveBoundParticleBindingEnergies)                                  \
    /* Whether to also  save mass, position, velocity (and internal energy) of 
     * all particles bound to subhaloes. */                                    \
    X(bool, SaveBoundParticlePotentialEnergies)                                \
    /* Whether to manually merge subhaloes that overlap in phase-space */      \
    X(bool, MergeTrappedSubhalos)                                              \
    /* If snapshots will not be analysed consecutively, this vector tells the 
     * code which numbers will be read. */                                     \
    X(vector<int>, SnapshotIdList)                                             \
    /* Particle types used to trace the host FoF of subhaloes (and orphans) */ \
    X(vector<int>, TracerParticleTypes)                                        \
    /* Particle types that cannot be subsampled when unbinding */              \
    X(vector<int>, DoNotSubsampleParticleTypes)                                \
    /* If snapshots are named differently within the same run. Otherwise better 
     * to use SnapshotFileBase. */                                             \
    X(vector<string>, SnapshotNameList)                                        \
    /* Mass threshold above which subhaloes are candidates to become centrals
     * of their FoF group. Expressed relative to the bound mass of the most
     * massive subhalo. */                                                     \
    X(HBTReal, MajorProgenitorMassRatio)                                       \
    /* Threshold ratio in the number of bound particles between consecutive
     * unbinding iterations to deem the unbinding as converged. */             \
    X(HBTReal, BoundMassPrecision)                                             \
    /* Maximum number of particles a subhalo can keep track of, relative to its
     * number of bound particles. */                                           \
    X(HBTReal, SourceSubRelaxFactor)                                           \
    X(HBTReal, TreeAllocFactor)                                                \
    X(HBTReal, TreeNodeOpenAngle)                                              \
    X(HBTInt, TreeMinNumOfCells)                                               \
    /* Value used to identify which particles do not belong to any FoF group */\
    X(HBTInt, ParticleNullGroupId)                                             \
    /* Maximum number of particles that can be subsampled used during unbinding.
     * Triggers subsampling if MaxSampleSizeOfPotentialEstimate is less than 
     * the number of particles that can be subsampled. */                      \
    X(HBTInt, MaxSampleSizeOfPotentialEstimate)                                \
    /* Compute the self-binding energy of the most bound subset of particles if     
      subsampling is enabled, to refine the subhalo coordinates. */            \
    X(bool, RefineMostBoundParticle)                                           \
    /* Minimum fraction of most bound particles used to refine the subhalo 
     * coordinates */                                                          \
    X(float, BoundFractionCenterRefinement)                                    \
    /* If baryonic particles can split. Relevant to swift simulations */       \
    X(int, ParticlesSplit)                                                     \
    /* Fraction of most bound particles used to estimate position and velocity
     * when doing merger checks. */                                            \
    X(HBTReal, SubCoreSizeFactor)                                              \
    /* Minimum number of most bound particles used to estimate position and 
     * velocity when doing merger checks. */                                   \
    X(HBTInt, SubCoreSizeMin)                                                  \

/* This macro contains the definition of configuration-related parameters that
 * are derived from user-defined values. The column format is as follows: 
 * Type, name, default value. */
#define DERIVED_PARAMETERS                                                     \
  /* Used to determine if mandatory parameters have been set. */               \
  X(vector<bool>,      IsSet   )                                               \
  /* Base units of the simulation */                                           \
  X(HBTReal,      MassInMsunh   )                                              \
  X(HBTReal,      LengthInMpch   )                                             \
  X(HBTReal,      VelInKmS   )                                                 \
  /* Bitmask used to identify particle types used to trace the host FoF of       
   * subhaloes (and orphans) */                                                \
  X(int, TracerParticleBitMask)                                                \
  /* Bitmask used to identify particle types that cannot be subsampled during 
   * unbinding */                                                              \
  X(int, DoNotSubsampleParticleBitMask)                                        \
  /*derived parameters; do not require user input*/                            \
  X(HBTReal,TreeNodeOpenAngleSquare)                                           \
  X(HBTReal,TreeNodeResolution)                                                \
  X(HBTReal,TreeNodeResolutionHalf)                                            \
  X(HBTReal,BoxHalf)                                                           \
  X(bool, GroupLoadedFullParticle)                                             \