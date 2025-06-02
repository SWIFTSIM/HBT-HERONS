#ifndef CONFIG_PARSER_H_INCLUDED
#define CONFIG_PARSER_H_INCLUDED

#include "../datatypes.h"
#include "../mpi_wrapper.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include "./parameters.h"

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

  /* Automatic member declaration based on what is in parameters.h */
#define DECLARE(type, name) type name;
  AVAILABLE_PARAMETERS(DECLARE)
  DERIVED_PARAMETERS(DECLARE)
#undef DECLARE

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
    RefineMostBoundParticle = true;
    BoundFractionCenterRefinement = 0.1; /* Default values chosen based on tests */
    GroupLoadedFullParticle = false;
    MaxPhysicalSofteningHalo = -1; // Indicates no max. physical softening is used.

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
