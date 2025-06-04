#include <cstdlib>
#include <unordered_set>
#include "config_parser.h"

namespace PhysicalConst
{
HBTReal G;
HBTReal H0;
} // namespace PhysicalConst

Parameter_t HBTConfig;

bool Parameter_t::TryCompulsoryParameterValue(string ParameterName, stringstream &ParameterValue)
{
#define TrySetPar(var, i)                                                                                              \
  if (ParameterName == #var)                                                                                           \
  {                                                                                                                    \
    ParameterValue >> var;                                                                                             \
    IsSet[i] = !ParameterValue.fail();                                                                                 \
    return true; /* Signals to not continue looking for matching parameter names */                                    \
  }

  TrySetPar(SnapshotPath, 0);
  TrySetPar(HaloPath, 1);
  TrySetPar(SubhaloPath, 2);
  TrySetPar(SnapshotFileBase, 3);
  TrySetPar(MaxSnapshotIndex, 4);
  TrySetPar(BoxSize, 5);
  TrySetPar(SofteningHalo, 6);

#undef TrySetPar

  return false; // This signals to continue looking for valid parameter names
}

bool Parameter_t::TrySingleValueParameter(string ParameterName, stringstream &ParameterValue)
{
#define TrySetPar(var)                                                                                                 \
  if (ParameterName == #var)                                                                                           \
  {                                                                                                                    \
    ParameterValue >> var;                                                                                             \
    return true; /* Signals to not continue looking for matching parameter names */                                    \
  }

  TrySetPar(SnapshotDirBase);
  TrySetPar(SnapshotFormat);
  TrySetPar(GroupFileFormat);
  TrySetPar(MaxConcurrentIO);
  TrySetPar(MinSnapshotIndex);
  TrySetPar(MinNumPartOfSub);
  TrySetPar(MinNumTracerPartOfSub);
  TrySetPar(NumTracerHostFinding);
  TrySetPar(NumTracersForDescendants);
  TrySetPar(ParticleIdNeedHash);
  TrySetPar(SnapshotIdUnsigned);
  TrySetPar(SaveBoundParticleProperties);
  TrySetPar(SaveBoundParticleBindingEnergies);
  TrySetPar(SaveBoundParticlePotentialEnergies);
  TrySetPar(MergeTrappedSubhalos);
  TrySetPar(PotentialEstimateUpscaleMassesPerType);
  TrySetPar(MajorProgenitorMassRatio);
  TrySetPar(BoundMassPrecision);
  TrySetPar(SourceSubRelaxFactor);
  TrySetPar(SubCoreSizeFactor);
  TrySetPar(SubCoreSizeMin);
  TrySetPar(TreeAllocFactor);
  TrySetPar(TreeNodeOpenAngle);
  TrySetPar(TreeMinNumOfCells);
  TrySetPar(MaxSampleSizeOfPotentialEstimate);
  TrySetPar(RefineMostBoundParticle);
  TrySetPar(BoundFractionCenterRefinement);
  TrySetPar(MassInMsunh);
  TrySetPar(LengthInMpch);
  TrySetPar(VelInKmS);
  TrySetPar(PeriodicBoundaryOn);
  TrySetPar(SnapshotHasIdBlock);
  TrySetPar(MaxPhysicalSofteningHalo);
  TrySetPar(ParticlesSplit);

#undef TrySetPar

  if (ParameterName == "GroupParticleIdMask")
  {
    ParameterValue >> hex >> GroupParticleIdMask >> dec;
    cout << "GroupParticleIdMask = " << hex << GroupParticleIdMask << dec << endl;
    return true;
  }

  return false; // This signals to continue looking for valid parameter names
}

bool Parameter_t::TryMultipleValueParameter(string ParameterName, stringstream &ParameterValues)
{

  if (ParameterName == "SnapshotIdList")
  {
    for (int i; ParameterValues >> i;)
      SnapshotIdList.push_back(i);
    return true;
  }
  if (ParameterName == "TracerParticleTypes")
  {
    /* Store as a vector first to output in Params.log in a human-readable
     * format */
    TracerParticleTypes.clear(); // To remove default values
    for (int i; ParameterValues >> i;)
      TracerParticleTypes.push_back(i);

    /* Create a bitmask, to be used internally by the code to identify valid
     * tracer particle types*/
    TracerParticleBitMask = 0;
    for (int i : TracerParticleTypes)
      TracerParticleBitMask += 1 << i;
    return true;
  }
  if (ParameterName == "DoNotSubsampleParticleTypes")
  {
    /* Store as a vector first to output in Params.log in a human-readable
     * format */
    DoNotSubsampleParticleTypes.clear(); // To remove default values
    for (int i; ParameterValues >> i;)
      DoNotSubsampleParticleTypes.push_back(i);

    /* Create a bitmask, to be used internally by the code to identify which 
     * particles it should not subsample during unbinding. */
    DoNotSubsampleParticleBitMask = 0;

    /* We specified -1, we allow any particle type to be subsampled. */
    auto it = std::find(DoNotSubsampleParticleTypes.begin(), DoNotSubsampleParticleTypes.end(), -1);
    if (it != DoNotSubsampleParticleTypes.end())
    {
      /* Sanity check: we can only specify a single number if -1 is present. */
      if(DoNotSubsampleParticleTypes.size() > 1)
      {
        stringstream error_message;
        error_message << "DoNotSubsampleParticleTypes only takes a single number if -1 is present. " << endl;
        throw runtime_error(error_message.str());
      }

      /* We need a right bit shift, as it is otherwise undefined. */
      DoNotSubsampleParticleBitMask = -1 >> 1;
      return true;
    }

    for (int i : DoNotSubsampleParticleTypes)
      DoNotSubsampleParticleBitMask += 1 << i;
    return true;
  }
  return false; // This signals to continue looking for valid parameter names
}

void Parameter_t::SetParameterValue(const string &line)
{
  /* Get the name of the parameter */
  stringstream ss(line);
  string name;
  ss >> name;
  //   transform(name.begin(),name.end(),name.begin(),::tolower);

  /* We will try matching the name of the input parameter to those which have
   * been defined in the code. If we find a match, we assign its value. First
   * we try compulsory parameters, and then optional ones */
  if (TryCompulsoryParameterValue(name, ss))
    return;
  if (TrySingleValueParameter(name, ss))
    return;
  if (TryMultipleValueParameter(name, ss))
    return;

  /* No matching parameter name has been found, throw an error message */
  stringstream error_message;
  error_message << "Unrecognized configuration entry: " << name << endl;
  throw runtime_error(error_message.str());
}

void Parameter_t::ParseConfigFile(const char *param_file)
{
  ifstream ifs;
  ifs.open(param_file);
  if (!ifs.is_open()) // or ifs.fail()
  {
    cerr << "Error: failed to open configuration: " << param_file << endl;
    exit(1);
  }
  vector<string> lines;
  string line;

  cout << "Reading configuration file " << param_file << endl;

  while (getline(ifs, line))
  {
    trim_trailing_garbage(line, "#[");
    trim_leading_garbage(line, " \t");
    if (!line.empty())
      SetParameterValue(line);
  }

  CheckParameters();

  PhysicalConst::G = 43.0071 * (MassInMsunh / 1e10) / VelInKmS / VelInKmS / LengthInMpch;
  PhysicalConst::H0 = 100. * (1. / VelInKmS) / (1. / LengthInMpch);

  /* Make particles split by default in swift (only relevant for hydro runs) */
  if ((SnapshotFormat == "swiftsim") & (ParticlesSplit == -1))
    ParticlesSplit = 1;
  else
    ParticlesSplit = 0;

  BoxHalf = BoxSize / 2.;
  TreeNodeResolution = SofteningHalo * 0.1;
  TreeNodeResolutionHalf = TreeNodeResolution / 2.;
  TreeNodeOpenAngleSquare = TreeNodeOpenAngle * TreeNodeOpenAngle;

  if (GroupFileFormat == "apostle_particle_index" || GroupFileFormat == "swiftsim_particle_index")
    GroupLoadedFullParticle = true;

  ReadSnapshotNameList();
}
void Parameter_t::ReadSnapshotNameList()
{ // to specify snapshotnamelist, create a file "snapshotlist.txt" under SubhaloPath, and list the filenames inside, one
  // per line.
  string snaplist_filename = SubhaloPath + "/snapshotlist.txt";
  ifstream ifs;
  ifs.open(snaplist_filename);
  if (ifs.is_open())
  {
    cout << "Found SnapshotNameList file " << snaplist_filename << endl;

    string line;
    while (getline(ifs, line))
    {
      trim_trailing_garbage(line, "#");
      istringstream ss(line);
      string name;
      ss >> name;
      if (!name.empty())
      {
        // 		cout<<name<<endl;
        SnapshotNameList.push_back(name);
      }
    }
  }
  if (SnapshotNameList.size())
    assert(SnapshotNameList.size() == MaxSnapshotIndex + 1);
}

/* Checks whether the all mandatory parameters have been specified */
void Parameter_t::CheckRequiredParameters()
{
  stringstream error_message; /* In case we need to raise an error. */

  if(!IsSet[0])
  {
    error_message << "SnapshotPath: Mandatory parameter not found." << endl;
    throw invalid_argument(error_message.str());
  }
  if(!IsSet[1])
  {
    error_message << "HaloPath: Mandatory parameter not found." << endl;
    throw invalid_argument(error_message.str());
  }
  if(!IsSet[2])
  {
    error_message << "SubhaloPath: Mandatory parameter not found." << endl;
    throw invalid_argument(error_message.str());
  }
  if(!IsSet[3])
  {
    error_message << "SnapshotFileBase: Mandatory parameter not found." << endl;
    throw invalid_argument(error_message.str());
  }
  if(!IsSet[4])
  {
    error_message << "MaxSnapshotIndex: Mandatory parameter not found." << endl;
    throw invalid_argument(error_message.str());
  }
  if(!IsSet[5])
  {
    error_message << "BoxSize: Mandatory parameter not found." << endl;
    throw invalid_argument(error_message.str());
  }
  if(!IsSet[6])
  {
    error_message << "SofteningHalo: Mandatory parameter not found." << endl;
    throw invalid_argument(error_message.str());
  }
}

/* Checks if we have duplicate values by creating a set, which will have the same
 * size as the input vector only if all its elements are unique. */
bool ContainsDuplicateValues(const std::vector<int> &InputVectorParameter)
{
  std::unordered_set<int> UniqueValues(InputVectorParameter.begin(), InputVectorParameter.end());
  return UniqueValues.size() != InputVectorParameter.size();
}

/* Checks whether the input parameters are valid */
void Parameter_t::CheckValidityParameters()
{
  stringstream error_message; /* In case we need to raise an error. */

  /* Make sure we have not provided duplicate entries for vector input types */
  {
    /* SnapshotIdList will be empty if no values are passed at runtime. */
    if(SnapshotIdList.size() && ContainsDuplicateValues(SnapshotIdList))
    {
      error_message << "SnapshotIdList: There are duplicate values, please remove repeated elements." << endl;
      throw invalid_argument(error_message.str());
    }

    if(ContainsDuplicateValues(TracerParticleTypes))
    {
      error_message << "TracerParticleTypes: There are duplicate values, please remove repeated elements." << endl;
      throw invalid_argument(error_message.str());
    }

    if(ContainsDuplicateValues(DoNotSubsampleParticleTypes))
    {
      error_message << "DoNotSubsampleParticleTypes: There are duplicate values, please remove repeated elements." << endl;
      throw invalid_argument(error_message.str());
    }
  }

  /* Checks relating to SnapshotIdList. */
  if (SnapshotIdList.size())
  {
    /* The snapshot ID list should be in ascending order. */
    if(!std::is_sorted(SnapshotIdList.begin(), SnapshotIdList.end()))
    {
      error_message << "SnapshotIdList: Values are not in ascending order." << endl;
      throw invalid_argument(error_message.str());
    }

    /* MaxSnapshotIndex should be consistent with SnapshotIdList */
    if (MaxSnapshotIndex != (SnapshotIdList.size() - 1))
    {
      error_message << "MaxSnapshotIndex (" << MaxSnapshotIndex << "): inconsistent with SnapshotIdList (SnapshotIdList.size() - 1 = " << (SnapshotIdList.size() - 1) << ")." << endl;
      throw invalid_argument(error_message.str());      
    }
  }

  /* No negative values for MaxSampleSizeOfPotentialEstimate allowed. */
  if(MaxSampleSizeOfPotentialEstimate < 0)
  {
    error_message << "MaxSampleSizeOfPotentialEstimate: Negative values are not allowed. Use a value of 0 (no subsampling) or larger (subsampling)." << endl;
    throw invalid_argument(error_message.str());
  }

  /* Cannot say we subsample in DMO and then disable DM particle subsampling. 
   * Conversely, we cannot subsample in hydro and then disable subsampling for all
   * relevant particle types */
  {
#ifdef DM_ONLY
    vector<int> TestParticleTypes = {1};
#else
    vector<int> TestParticleTypes = {0, 1, 4, 5};
#endif

    /* Create a bit mask corresponding to disabling subsampling for all relevant types */    
    int TestBitMask = 0;
    for (int i : TestParticleTypes)
      TestBitMask += 1 << i;

    /* Is this the same as the mask generated from the user input? */
    if (TestBitMask == DoNotSubsampleParticleBitMask)
    {
#ifdef DM_ONLY
      error_message << "Cannot enable subsampling (MaxSampleSizeOfPotentialEstimate > 0) and then disable subsampling (DoNotSubsampleParticleTypes 1) in DMO simulations. Disable either of the two options." << endl;
#else
      error_message << "Cannot enable subsampling (MaxSampleSizeOfPotentialEstimate > 0) and then disable subsampling (DoNotSubsampleParticleTypes 0 1 4 5) in HYDRO simulations. Disable either of the two options." << endl;
#endif
      throw invalid_argument(error_message.str());
    }
  }
}

/* Checks if the input parameter file contains all required parameters and that 
 * they have valid values. */
void Parameter_t::CheckParameters()
{
  CheckRequiredParameters();
  CheckValidityParameters();
}

void ParseHBTParams(int argc, char **argv, Parameter_t &config, int &snapshot_start, int &snapshot_end)
{
  if (argc < 2)
  {
    cerr << "Usage: " << argv[0] << " [param_file] <snapshot_start> <snapshot_end>\n";
    exit(1);
  }
  config.ParseConfigFile(argv[1]);
  if (2 == argc)
  {
    snapshot_start = config.MinSnapshotIndex;
    snapshot_end = config.MaxSnapshotIndex;
  }
  else
  {
    snapshot_start = atoi(argv[2]);
    if (argc > 3)
      snapshot_end = atoi(argv[3]);
    else
      snapshot_end = snapshot_start;
  }
  cout << "Running " << argv[0] << " from snapshot " << snapshot_start << " to " << snapshot_end
       << " using configuration file " << argv[1] << endl;
}

void Parameter_t::BroadCast(MpiWorker_t &world, int root)
/*sync parameters and physical consts across*/
{
#define _SyncVec(x, t) world.SyncContainer(x, t, root)
#define _SyncAtom(x, t) world.SyncAtom(x, t, root)
#define _SyncBool(x) world.SyncAtomBool(x, root)
#define _SyncVecBool(x) world.SyncVectorBool(x, root)
#define _SyncReal(x) _SyncAtom(x, MPI_HBT_REAL)

  _SyncVec(SnapshotPath, MPI_CHAR);
  _SyncVec(HaloPath, MPI_CHAR);
  _SyncVec(SubhaloPath, MPI_CHAR);
  _SyncVec(SnapshotFileBase, MPI_CHAR);
  _SyncAtom(MaxSnapshotIndex, MPI_INT);
  _SyncReal(BoxSize);
  _SyncReal(SofteningHalo);
  _SyncReal(MaxPhysicalSofteningHalo);
  _SyncVecBool(IsSet);

  _SyncVec(SnapshotDirBase, MPI_CHAR);
  _SyncVec(SnapshotFormat, MPI_CHAR);
  _SyncVec(GroupFileFormat, MPI_CHAR);
  _SyncAtom(MaxConcurrentIO, MPI_INT);
  _SyncAtom(MinSnapshotIndex, MPI_INT);
  _SyncAtom(MinNumPartOfSub, MPI_INT);
  _SyncAtom(MinNumTracerPartOfSub, MPI_INT);
  _SyncAtom(NumTracerHostFinding, MPI_INT);
  _SyncAtom(NumTracersForDescendants, MPI_INT);
  _SyncAtom(GroupParticleIdMask, MPI_LONG);
  _SyncReal(MassInMsunh);
  _SyncReal(LengthInMpch);
  _SyncReal(VelInKmS);
  _SyncBool(PeriodicBoundaryOn);
  _SyncBool(SnapshotHasIdBlock);
  _SyncBool(ParticleIdNeedHash);
  _SyncBool(SnapshotIdUnsigned);
  _SyncBool(SaveBoundParticleProperties);
  _SyncBool(SaveBoundParticleBindingEnergies);
  _SyncBool(SaveBoundParticlePotentialEnergies);
  _SyncBool(MergeTrappedSubhalos);
  _SyncBool(PotentialEstimateUpscaleMassesPerType);
  _SyncVec(SnapshotIdList, MPI_INT);
  world.SyncVectorString(SnapshotNameList, root);

  _SyncReal(MajorProgenitorMassRatio);
#ifdef ALLOW_BINARY_SYSTEM
  _SyncReal(BinaryMassRatioLimit);
#endif
  _SyncReal(BoundMassPrecision);
  _SyncReal(SourceSubRelaxFactor);
  _SyncReal(SubCoreSizeFactor);
  _SyncAtom(SubCoreSizeMin, MPI_HBT_INT);

  _SyncReal(TreeAllocFactor);
  _SyncReal(TreeNodeOpenAngle);
  _SyncAtom(TreeMinNumOfCells, MPI_HBT_INT);

  _SyncAtom(MaxSampleSizeOfPotentialEstimate, MPI_HBT_INT);
  _SyncBool(RefineMostBoundParticle);
  _SyncReal(BoundFractionCenterRefinement);

  _SyncReal(TreeNodeOpenAngleSquare);
  _SyncReal(TreeNodeResolution);
  _SyncReal(TreeNodeResolutionHalf);
  _SyncReal(BoxHalf);

  _SyncAtom(ParticleNullGroupId,
            MPI_HBT_INT); // Sync here for consistency, but uninitialised until read from snapshots.

  _SyncBool(GroupLoadedFullParticle);
  _SyncAtom(TracerParticleBitMask, MPI_INT);
  _SyncAtom(DoNotSubsampleParticleBitMask, MPI_INT);
  _SyncAtom(ParticlesSplit, MPI_INT);
  //---------------end sync params-------------------------//

  _SyncReal(PhysicalConst::G);
  _SyncReal(PhysicalConst::H0);

#undef _SyncVec
#undef _SyncAtom
#undef _SyncBool
#undef _SyncVecBool
#undef _SyncReal
}
void Parameter_t::DumpParameters()
{
#ifndef HBT_VERSION
  cout << "Warning: HBT_VERSION unknown.\n";
#define HBT_VERSION "unknown"
#endif
  string filename = SubhaloPath + "/Parameters.log";
  ofstream version_file(filename, ios::out | ios::trunc);
  if (!version_file.is_open())
  {
    cerr << "Error opening " << filename << " for parameter dump." << endl;
    exit(1);
  }
  version_file << "#VERSION " << HBT_VERSION << endl;

#define DumpPar(var) version_file << #var << "  " << var << endl;
#define DumpComment(var)                                                                                               \
  {                                                                                                                    \
    version_file << "#";                                                                                               \
    DumpPar(var);                                                                                                      \
  }
#define DumpHeader(title)                                                                                              \
  {                                                                                                                    \
    version_file << endl;                                                                                              \
    version_file << "[" << title << "]" << endl;                                                                       \
  }

  DumpHeader("File paths");
  DumpPar(SnapshotPath);
  DumpPar(HaloPath);
  DumpPar(SubhaloPath);
  DumpPar(SnapshotFileBase);
  DumpPar(SnapshotDirBase);

  DumpHeader("File IO");
  DumpPar(SnapshotFormat);
  DumpPar(GroupFileFormat);
  DumpPar(MaxConcurrentIO);
  DumpPar(MinSnapshotIndex);
  DumpPar(MaxSnapshotIndex);
  if (SnapshotIdList.size())
  {
    version_file << "SnapshotIdList";
    for (auto &&i : SnapshotIdList)
      version_file << " " << i;
    version_file << endl;
  }

  if (SnapshotNameList.size())
  {
    version_file << "#SnapshotNameList";
    for (auto &&i : SnapshotNameList)
      version_file << " " << i;
    version_file << endl;
  }

  DumpHeader("Particle Properties");
  DumpComment(GroupLoadedFullParticle);
  DumpPar(SnapshotHasIdBlock);
  DumpPar(ParticleIdNeedHash);
  DumpPar(SnapshotIdUnsigned);
  DumpPar(SaveBoundParticleProperties);
  DumpPar(SaveBoundParticleBindingEnergies);
  DumpPar(SaveBoundParticlePotentialEnergies);
#ifndef DM_ONLY
  DumpPar(ParticlesSplit);
#endif
  if (GroupParticleIdMask)
    version_file << "GroupParticleIdMask " << hex << GroupParticleIdMask << dec << endl;

  DumpHeader("Units");
  DumpPar(MassInMsunh);
  DumpPar(LengthInMpch);
  DumpPar(VelInKmS);

  DumpHeader("Gravity tree");
  DumpPar(TreeAllocFactor);
  DumpPar(TreeNodeOpenAngle);
  DumpPar(TreeMinNumOfCells);

  DumpHeader("Gravity softening");
  DumpPar(SofteningHalo);
  DumpPar(MaxPhysicalSofteningHalo);

  DumpHeader("Simulation Box");
  DumpPar(BoxSize);
  DumpPar(PeriodicBoundaryOn);

  DumpHeader("Subhalo Unbinding");
  DumpPar(SubCoreSizeMin);
  DumpPar(SubCoreSizeFactor);
  DumpPar(BoundMassPrecision);
  DumpPar(SourceSubRelaxFactor);
  DumpPar(RefineMostBoundParticle);
  if (RefineMostBoundParticle)
    DumpPar(BoundFractionCenterRefinement);
  DumpPar(MaxSampleSizeOfPotentialEstimate);
  DumpPar(PotentialEstimateUpscaleMassesPerType);

  DumpHeader("Subhalo Tracking");
  DumpPar(MinNumPartOfSub);
  DumpPar(MinNumTracerPartOfSub);
  DumpPar(NumTracerHostFinding);
  DumpPar(NumTracersForDescendants);
  if (TracerParticleTypes.size())
  {
    version_file << "TracerParticleTypes";
    for (auto &&i : TracerParticleTypes)
      version_file << " " << i;
    version_file << endl;
  }
  if (DoNotSubsampleParticleTypes.size())
  {
    version_file << "DoNotSubsampleParticleTypes";
    for (auto &&i : DoNotSubsampleParticleTypes)
      version_file << " " << i;
    version_file << endl;
  }
  DumpPar(MergeTrappedSubhalos);
  DumpPar(MajorProgenitorMassRatio);

#undef DumpPar
#undef DumpHeader
#undef DumpComment
  version_file.close();
}

HBTReal Parameter_t::GetCurrentSoftening(HBTReal ScaleFactor)
{
  // Only one softening defined, use comoving
  if (MaxPhysicalSofteningHalo == -1)
    return SofteningHalo;

  // If two are defined, choose the one with the smallest value
  return min(SofteningHalo, MaxPhysicalSofteningHalo / ScaleFactor);
}
