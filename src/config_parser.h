#ifndef CONFIG_PARSER_H_INCLUDED
#define CONFIG_PARSER_H_INCLUDED

#include "datatypes.h"
#include "mpi_wrapper.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#define HBT_VERSION "1.16.1.MPI"

namespace PhysicalConst
{ // initialized after reading parameter file.
extern HBTReal G;
extern HBTReal H0;
} // namespace PhysicalConst

#define NumberOfCompulsaryConfigEntries 7
class Parameter_t
{ /*!remember to register members in BroadCast() and SetParameterValue() functions if you change them!*/
public:
  // remember to update SetParameterValue() and DumpParameters() accordingly if you change any parameter definition.
  /*compulsory parameters*/
  string SnapshotPath;
  string HaloPath;
  string SubhaloPath;
  string SnapshotFileBase;
  int MaxSnapshotIndex;
  HBTReal BoxSize; // to check the unit of snapshot according to the BoxSize in header
  HBTReal SofteningHalo;
  HBTReal MaxPhysicalSofteningHalo;
  vector<bool> IsSet;

  /*optional*/
  string SnapshotDirBase;
  string SnapshotFormat;
  string GroupFileFormat;
  int MaxConcurrentIO;
  int MinSnapshotIndex;
  int MinNumPartOfSub;
  int MinNumTracerPartOfSub;
  int NumTracerHostFinding;
  int NumTracersForDescendants;
  long GroupParticleIdMask; // only used for a peculiar gadget format.
  HBTReal MassInMsunh;
  HBTReal LengthInMpch;
  HBTReal VelInKmS;
  bool PeriodicBoundaryOn;
  bool SnapshotHasIdBlock; // set to False when your snapshot is sorted according to particle id so that no id block is
                           // present.
  bool ParticleIdNeedHash;
  bool SnapshotIdUnsigned;
  bool SaveBoundParticleProperties;
  bool SaveBoundParticleBindingEnergies;
  bool SaveBoundParticlePotentialEnergies;
  bool MergeTrappedSubhalos; // whether to MergeTrappedSubhalos, see code paper for more info.
  bool PotentialEstimateUpscaleMassesPerType;

  vector<int> SnapshotIdList;
  vector<int> TracerParticleTypes;
  vector<int> DoNotSubsampleParticleTypes; /* Which particles types cannot be subsampled. */
  vector<string> SnapshotNameList;

  HBTReal MajorProgenitorMassRatio;
  HBTReal BoundMassPrecision;
  HBTReal SourceSubRelaxFactor;
  HBTReal SubCoreSizeFactor; // coresize=Nbound*CoreSizeFactor, to get center coordinates for the KineticDistance test.
  HBTInt SubCoreSizeMin;     // Minimum coresize

  HBTReal TreeAllocFactor;
  HBTReal TreeNodeOpenAngle;
  HBTInt TreeMinNumOfCells;
  HBTInt ParticleNullGroupId;

  HBTInt MaxSampleSizeOfPotentialEstimate;
  bool RefineMostBoundParticle; // whether to further improve mostbound particle accuracy in case a
                                // MaxSampleSizeOfPotentialEstimate is used. this introduces some overhead if true, but
                                // leads to more accuracy mostbound particle
  float BoundFractionCenterRefinement; /* The fraction of most bound particles to use if the most bound particle will be
                                          refined */

  int TracerParticleBitMask; /* Bitmask used to identify which particle type can be used as tracer */
  int DoNotSubsampleParticleBitMask; /* Bitmask used to identify which particle type cannot be subsampled */
  int ParticlesSplit;        /* Whether baryonic particles are able to split. Relevant to swift simulations */

  /*derived parameters; do not require user input*/
  HBTReal TreeNodeOpenAngleSquare;
  HBTReal TreeNodeResolution;
  HBTReal TreeNodeResolutionHalf;
  HBTReal BoxHalf;
  bool GroupLoadedFullParticle; // whether group particles are loaded with full particle properties or just ids.
  int ReassignParticles; // Move particles between subhalos based on neighbour search
  int NumNeighboursForReassignment; // how many neighbours to search for

  Parameter_t() : IsSet(NumberOfCompulsaryConfigEntries, false), SnapshotIdList(), SnapshotNameList()
  {
    SnapshotDirBase = "";
    SnapshotFormat = "gadget"; // see example config file for alternative formats
    GroupFileFormat = "gadget3_int";
    MaxConcurrentIO = 10;
    MinSnapshotIndex = 0;
    MinNumPartOfSub = 20;
    MinNumTracerPartOfSub = 10;
    NumTracerHostFinding = MinNumTracerPartOfSub;
    NumTracersForDescendants = MinNumTracerPartOfSub;
    GroupParticleIdMask = 0;
    MassInMsunh = 1e10;
    LengthInMpch = 1;
    VelInKmS = 1.;
    PeriodicBoundaryOn = true;
    SnapshotHasIdBlock = true;
    ParticleIdNeedHash = true;
    SnapshotIdUnsigned = false;
    SaveBoundParticleProperties = false;
    SaveBoundParticleBindingEnergies = false;
    SaveBoundParticlePotentialEnergies = false;
#ifdef NO_STRIPPING
    MergeTrappedSubhalos = false;
#else
    MergeTrappedSubhalos = true;
#endif
    MajorProgenitorMassRatio = 0.8;
    BoundMassPrecision = 0.995;
    SourceSubRelaxFactor = 3.;
    SubCoreSizeFactor = 0.25;
    SubCoreSizeMin = 20;
    TreeAllocFactor = 0.8; /* a value of 2 should be more than sufficient*/
    TreeNodeOpenAngle = 0.45;
    TreeMinNumOfCells = 10;
    MaxSampleSizeOfPotentialEstimate = 1000; // set to 0 to disable sampling
    PotentialEstimateUpscaleMassesPerType = 1;
    RefineMostBoundParticle = true;
    BoundFractionCenterRefinement = 0.1; /* Default values chosen based on tests */
    GroupLoadedFullParticle = false;
    MaxPhysicalSofteningHalo = -1; // Indicates no max. physical softening is used.

    ParticleNullGroupId = -1; /* Value of FoF group corresponding to no FoF */

    /* Tracer-related parameters. If unset, only use collisionless particles (DM
     * + Stars) as tracer. Here we assume they correspond to particle types 1
     * and 4, respectively. */
    TracerParticleTypes = vector<int>{1, 4};
    TracerParticleBitMask = 0;
    for (int i : TracerParticleTypes)
      TracerParticleBitMask += 1 << i;

    /* We default to not subsampling black holes, since they are generally very
     * massive relative to all other particle types. */
    DoNotSubsampleParticleTypes = vector<int>{5};
    DoNotSubsampleParticleBitMask = 0;
    for (int i : DoNotSubsampleParticleTypes)
      DoNotSubsampleParticleBitMask += 1 << i;
 
    /* The value is negative to indicate whether the parameter has been set in the. If not,
     * we will default to a value of 1 if this is a swift HYDRO run. This way we reminder the
     * user to pre-process snapshots (toolbox/swiftsim/generate_splitting_information.py) */
    ParticlesSplit = -1;

    /* Whether we move particles between subhalos based on nearest neighbours:
       0 = do nothing
       1 = non-tracer type particles are moved based on nearest tracer type particles
     */
#ifdef DM_ONLY
    ReassignParticles = 0;
#else
    ReassignParticles = 1; // Enabled by default in hydro runs
#endif
    NumNeighboursForReassignment = 10;
  }
  void ReadSnapshotNameList();
  void ParseConfigFile(const char *param_file);
  void SetParameterValue(const string &line);

  /* Functions that will check if the input parameter file contains all required
   * parameters and that they have valid values. */
  void CheckParameters();
  void CheckRequiredParameters();
  void CheckValidityParameters();

  void BroadCast(MpiWorker_t &world, int root);
  void BroadCast(MpiWorker_t &world, int root, int &snapshot_start, int &snapshot_end)
  {
    BroadCast(world, root);
    world.SyncAtom(snapshot_start, MPI_INT, root);
    world.SyncAtom(snapshot_end, MPI_INT, root);
  }
  void DumpParameters();

  // Decides whether to use comoving or max. physical softening (if defined).
  HBTReal GetCurrentSoftening(HBTReal ScaleFactor); // NOTE: perhaps makes more sense to include in gravity tree...

private:
  bool TryCompulsoryParameterValue(string ParameterName, stringstream &ParameterValue);
  bool TrySingleValueParameter(string ParameterName, stringstream &ParameterValue);
  bool TryMultipleValueParameter(string ParameterName, stringstream &ParameterValues);
};

extern Parameter_t HBTConfig;
extern void ParseHBTParams(int argc, char **argv, Parameter_t &config, int &snapshot_start, int &snapshot_end);
inline void trim_leading_garbage(string &s, const string &garbage_list)
{
  int pos = s.find_first_not_of(garbage_list); // look for any good staff
  if (string::npos != pos)
    s.erase(0, pos); // s=s.substr(pos);
  else               // no good staff, clear everything
    s.clear();
}
inline void trim_trailing_garbage(string &s, const string &garbage_list)
{
  int pos = s.find_first_of(garbage_list);
  if (string::npos != pos)
    s.erase(pos);
}

#define NEAREST(x)                                                                                                     \
  (((x) > HBTConfig.BoxHalf) ? ((x)-HBTConfig.BoxSize) : (((x) < -HBTConfig.BoxHalf) ? ((x) + HBTConfig.BoxSize) : (x)))
inline HBTReal PeriodicDistance(const HBTxyz &x, const HBTxyz &y)
{
  HBTxyz dx;
  dx[0] = x[0] - y[0];
  dx[1] = x[1] - y[1];
  dx[2] = x[2] - y[2];
  if (HBTConfig.PeriodicBoundaryOn)
  {
    dx[0] = NEAREST(dx[0]);
    dx[1] = NEAREST(dx[1]);
    dx[2] = NEAREST(dx[2]);
  }
  return sqrt(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);
}
#endif
