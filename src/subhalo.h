#ifndef SUBHALO_HEADER_INCLUDED
#define SUBHALO_HEADER_INCLUDED

#include <iostream>
#include <new>
#include <vector>
#include <cassert>
#include <cstdlib>

#include "datatypes.h"
#include "halo.h"
#include "hdf_wrapper.h"
#include "mpi_wrapper.h"
#include "snapshot_number.h"

enum class SubReaderDepth_t
{
  SubTable,     // read only the properties of subhalos
  SubParticles, // read sub particle list, plus properties of sub
  SrcParticles  // read src particle list instead of sub particles, plus properties of sub
};

class Subhalo_t;
typedef vector<Subhalo_t> SubhaloList_t;

class Subhalo_t
{
public:
  typedef vector<Particle_t> ParticleList_t;
  typedef vector<HBTInt> SubIdList_t;
  HBTInt TrackId;
  HBTInt Nbound;
  float Mbound;
#ifndef DM_ONLY
  HBTInt NboundType[TypeMax];
  float MboundType[TypeMax];
#endif
  HBTInt TracerIndex; // Index of the most bound collisionless particle, used as tracer.
#ifdef CHECK_TRACER_INDEX
  HBTInt TracerId; // Id of the particle with index TracerIndex
#endif
  HBTInt HostHaloId;
  HBTInt Rank; // 0 for central and field subs, >0 for satellites
  int Depth;   // depth of the subhalo: central=0, sub=1, sub-sub=2, ...
  float LastMaxMass;
  int SnapshotIndexOfLastMaxMass; // the snapshot when it has the maximum subhalo mass, only considering past snapshots.
  int SnapshotIndexOfLastIsolation; // the last snapshot when it was a central, only considering past snapshots. -1 if
                                    // the subhalo has always been a central

  int SnapshotIndexOfBirth; // when the subhalo first becomes resolved
  int SnapshotIndexOfDeath; // when the subhalo first becomes un-resolved; only set if
                            // currentsnapshot>=SnapshotIndexOfDeath.
  int SnapshotIndexOfSink;  // when the subhalo sinks

  // profile properties
  float RmaxComoving;
  float RmaxComovingOfLastMaxVmax;
  float VmaxPhysical;
  float LastMaxVmaxPhysical;
  int SnapshotIndexOfLastMaxVmax; // the snapshot when it has the maximum Vmax, only considering past snapshots.

  float REncloseComoving; // Radius of minimum sphere which contains all bound particles
  float RHalfComoving;

  // SO properties using subhalo particles alone, meant for quick and dirty calculations
  float BoundR200CritComoving;
  //   float R200MeanComoving;
  //   float RVirComoving;
  float BoundM200Crit;
  //   float M200Mean;
  //   float MVir;

  // kinetic properties
  float SpecificSelfPotentialEnergy; // average specific potential energy of each particle, <phi/m>. the total potential
                                     // energy of the system is 0.5*<phi/m>*M, where the 0.5 factor corrects for double
                                     // counting of mutual potential.
  float SpecificSelfKineticEnergy;   //<0.5*v^2>
  float SpecificAngularMomentum[3];  //<Rphysical x Vphysical>

  // shapes
#ifdef HAS_GSL
  float InertialEigenVector[3][3]; // three float[3] vectors, with amplitude equal to eigenvalue.
  float InertialEigenVectorWeighted[3][3];
#endif
  float InertialTensor[6]; //{Ixx, Ixy, Ixz, Iyy, Iyz, Izz}
  float InertialTensorWeighted[6];

  HBTxyz ComovingAveragePosition;
  HBTxyz PhysicalAverageVelocity;   // default vel of sub
  HBTxyz ComovingMostBoundPosition; // default pos of sub
  HBTxyz PhysicalMostBoundVelocity;
  HBTInt MostBoundParticleId; // a shortcut to the first particle in Particles, for easier IO in GALFORM

  // for merging
  HBTInt SinkTrackId;         // the trackId it sinked to, -1 if it hasn't sunk.
  HBTInt DescendantTrackId;   // the trackId containing a subset of the particles most bound to this object in the
                              // previous output
  HBTInt NestedParentTrackId; // the trackID of the subhalo containing this subhalo, or -1 for top level subhalos

  ParticleList_t Particles;

  /* Binding/potential energies of the particles bound to the subhalo. They are a separate instance from ParticleList_t
   * because they do not need communicating before unbinding. */
  vector<float> ParticleBindingEnergies;
  vector<float> ParticlePotentialEnergies;

  SubIdList_t NestedSubhalos; // list of sub-in-subs.

  /* Methods relating to new merging approach */
  bool MergeRecursively(SubhaloList_t &Subhalos, const Snapshot_t &snap, Subhalo_t &ReferenceSubhalo);
  bool AreOverlappingInPhaseSpace(const Subhalo_t &ReferenceSubhalo);
  float PhaseSpaceDistance(const Subhalo_t &ReferenceSubhalo);
  void GetCorePhaseSpaceProperties();
  void SetMergerInformation(const HBTInt &ReferenceTrackId, const int &SnapshotIndex);

  /* Properties relating to the new merging approach */
  HBTxyz CoreComovingPosition;
  HBTxyz CorePhysicalVelocity;
  float CoreComovingSigmaR;
  float CorePhysicalSigmaV;

  Subhalo_t()
    : Nbound(0), Rank(0), Mbound(0), Depth(0)
#ifndef DM_ONLY
      ,
      NboundType{0}, MboundType{0.}
#endif
  {
    TrackId = SpecialConst::NullTrackId;
    SnapshotIndexOfLastIsolation = SpecialConst::NullSnapshotId;
    SnapshotIndexOfLastMaxMass = SpecialConst::NullSnapshotId;
    LastMaxMass = 0.;
    RmaxComovingOfLastMaxVmax = 0.;
    LastMaxVmaxPhysical = 0.;
    SnapshotIndexOfLastMaxVmax = SpecialConst::NullSnapshotId;
    SnapshotIndexOfBirth = SpecialConst::NullSnapshotId;
    SnapshotIndexOfDeath = SpecialConst::NullSnapshotId;
    SnapshotIndexOfSink = SpecialConst::NullSnapshotId;
    SinkTrackId = SpecialConst::NullTrackId;
    DescendantTrackId = SpecialConst::NullTrackId;
    MostBoundParticleId = SpecialConst::NullParticleId;
    TracerIndex = -1;
  }
  void Unbind(const Snapshot_t &epoch);
  void RecursiveUnbind(SubhaloList_t &Subhalos, const Snapshot_t &snap);
  HBTReal KineticDistance(const Halo_t &halo, const Snapshot_t &epoch);
  void TruncateSource();
  void RemoveOtherHostParticles(const HBTInt &GlobalHostHaloId);
  float GetMass() const
  {
    return Mbound; // accumulate(begin(MboundType), end(MboundType), (HBTReal)0.);
  }
  void UpdateTrack(const Snapshot_t &epoch);
  bool IsCentral()
  {
    return 0 == Rank;
  }
  void CalculateProfileProperties(const Snapshot_t &epoch);
  void CalculateShape();
  void AverageCoordinates();
  void CountParticleTypes();
  HBTInt KickNullParticles();
  void CountParticles();
  void LevelUpDetachedMembers(vector<Subhalo_t> &Subhalos);
  // for merger
  void MergeTo(Subhalo_t &host);
  bool IsTrapped()
  {
    return SinkTrackId != SpecialConst::NullTrackId;
  }
  bool IsAlive()
  {
    return SnapshotIndexOfDeath == SpecialConst::NullSnapshotId;
  }
  bool JustTrapped(int currentsnapshotindex)
  {
    return SnapshotIndexOfSink == currentsnapshotindex;
  }
  HBTInt GetTracerIndex()
  {
    assert(TracerIndex >= 0);
    assert(TracerIndex < Particles.size());
#ifdef CHECK_TRACER_INDEX
    if (TracerId != Particles[TracerIndex].Id)
    {
      cerr << "Tracer particle ID is incorrect!" << endl;
      abort();
    }
#endif
    return TracerIndex;
  }
  void SetTracerIndex(const HBTInt index)
  {
    assert((index >= 0));
    assert((index < Particles.size()) ||
           ((Particles.size() == 0) && (index == 0))); // allow index=0 for halos with no particles
    TracerIndex = index;
#ifdef CHECK_TRACER_INDEX
    if (Particles.size() > 0)
    {
      TracerId = Particles[index].Id;
    }
#endif
  }
  vector<HBTInt> GetMostBoundTracerIds(HBTInt n);
};

class MemberShipTable_t
/* list the subhaloes inside each host, rather than ordering the subhaloes
 *
 * the principle is to not move the objects, but construct a table of them, since moving objects will change their id
 * (or index at least), introducing the trouble to re-index them and update the indexes in any existence references.
 */
{
public:
  typedef VectorView_t<HBTInt> MemberList_t; // list of members in a group
private:
  void BindMemberLists();
  void FillMemberLists(const SubhaloList_t &Subhalos, bool include_orphans);
  void CountMembers(const SubhaloList_t &Subhalos, bool include_orphans);
  void SortSatellites(const SubhaloList_t &Subhalos);
  void CountEmptyGroups();
  /*avoid operating on the Mem_* below; use the public VectorViews whenever possible; only operate the Mem_* variables
   * when adjusting memory*/
  vector<MemberList_t> Mem_SubGroups; // list of subhaloes inside each host halo, with the storage of each subgroup
                                      // mapped to a location in AllMembers
public:
  vector<HBTInt> AllMembers; // the complete list of all the subhaloes in SubGroups. (contains local subhaloid, i.e.,
                             // the index of subhaloes in the local SubhaloSnapshot_t). do not resize this vector
                             // manually. resize only with ResizeAllMembers() function.
  VectorView_t<MemberList_t>
    SubGroups;   // list of subhaloes inside each host halo. contain one more group than halo catalogue, to hold field
                 // subhaloes. It is properly offseted so that SubGroup[hostid=-1] gives field subhaloes, and hostid>=0
                 // for the normal groups.
  HBTInt NBirth; // newly born halos, excluding fake halos
  HBTInt NFake;  // Fake (unbound) halos with no progenitors
  vector<vector<HBTInt>> SubGroupsOfHeads; // list of top-level subhaloes in each halo

  MemberShipTable_t() : Mem_SubGroups(), AllMembers(), SubGroups(), SubGroupsOfHeads(), NBirth(0), NFake(0)
  {
  }
  HBTInt GetNumberOfFieldSubs()
  {
    return SubGroups[-1].size();
  }
  void Init(const HBTInt nhalos, const HBTInt nsubhalos, const float alloc_factor = 1.2);
  void ResizeAllMembers(size_t n);
  void Build(const HBTInt nhalos, const SubhaloList_t &Subhalos, bool include_orphans);
  void SortMemberLists(const SubhaloList_t &Subhalos);
  void AssignRanks(SubhaloList_t &Subhalos);
  void SubIdToTrackId(const SubhaloList_t &Subhalos);
  void TrackIdToSubId(SubhaloList_t &Subhalos);
};
class SubhaloSnapshot_t : public Snapshot_t
{
private:
  bool ParallelizeHaloes;
  hid_t H5T_SubhaloInMem, H5T_SubhaloInDisk;
  MPI_Datatype MPI_HBT_SubhaloShell_t; // MPI datatype ignoring the particle list

  void RegisterNewTracks(MpiWorker_t &world);
  void DecideCentrals(const HaloSnapshot_t &halo_snap);
  void FeedCentrals(HaloSnapshot_t &halo_snap);
  void BuildHDFDataType();
  void BuildMPIDataType();
  void PurgeMostBoundParticles();

  /* I/O methods */
  void ReadFile(int iFile, const SubReaderDepth_t depth);
  void WriteBoundFiles(MpiWorker_t &world, const int &number_ranks_writing);
  void WriteSourceFiles(MpiWorker_t &world, const int &number_ranks_writing);
  void WriteBoundSubfile(int iFile, int nfiles, HBTInt NumSubsAll);
  void WriteSourceSubfile(int iFile, int nfiles);

  void LevelUpDetachedSubhalos();
  void ExtendCentralNest();
  void LocalizeNestedIds(MpiWorker_t &world);
  void GlobalizeTrackReferences();
  void NestSubhalos(MpiWorker_t &world);
  void MaskSubhalos();
  vector<int> RootNestSize; // buffer variable for temporary use.
  void FillDepthRecursive(HBTInt subid, int depth);
  void FillDepth();
  void SetNestedParentIds();

  void HandleTracerlessSubhalos(MpiWorker_t &world, vector<Subhalo_t> &LocalSubhalos);

public:
  SubhaloList_t Subhalos;
  MemberShipTable_t MemberTable;

  SubhaloSnapshot_t() : Snapshot_t(), Subhalos(), MemberTable(), ParallelizeHaloes(true)
  {
    BuildHDFDataType();
    BuildMPIDataType();
  }
  SubhaloSnapshot_t(MpiWorker_t &world, int snapshot_index,
                    const SubReaderDepth_t depth = SubReaderDepth_t::SubParticles)
    : SubhaloSnapshot_t()
  {
    Load(world, snapshot_index, depth);
  }
  ~SubhaloSnapshot_t()
  {
    // This destructor may be called after HDF5 and MPI have shut down
    if (H5Iis_valid(H5T_SubhaloInDisk))
      H5Tclose(H5T_SubhaloInDisk);
    if (H5Iis_valid(H5T_SubhaloInMem))
      H5Tclose(H5T_SubhaloInMem);
    My_Type_free(&MPI_HBT_SubhaloShell_t);
  }
  string GetSubDir();
  void GetSubFileName(string &filename, int iFile, const string &ftype = "Sub");
  void Load(MpiWorker_t &world, int snapshot_index, const SubReaderDepth_t depth = SubReaderDepth_t::SubParticles);
  void Save(MpiWorker_t &world);
  void Clear()
  {
    // TODO
    cout << "Clean() not implemented yet\n";
  }
  void UpdateParticles(MpiWorker_t &world, const ParticleSnapshot_t &snapshot);
  void UpdateSplitParticles(const ParticleSnapshot_t &snapshot);
  //   void ParticleIndexToId();
  void UpdateMostBoundPosition(MpiWorker_t &world, const ParticleSnapshot_t &part_snap);
  void AssignHosts(MpiWorker_t &world, HaloSnapshot_t &halo_snap, const ParticleSnapshot_t &part_snap);
  void ConstrainToSingleHost(const HaloSnapshot_t &halo_snap);
  void PrepareCentrals(MpiWorker_t &world, HaloSnapshot_t &halo_snap);
  void RefineParticles();
  void ReassignParticles(MpiWorker_t &world, HaloSnapshot_t &halo_snap);
  void UpdateTracks(MpiWorker_t &world, const HaloSnapshot_t &halo_snap);

  /* To remove duplicate particles from the source subgroup. */
  void CleanTracks();

  HBTInt size() const
  {
    return Subhalos.size();
  }
  HBTInt GetId(HBTInt index) const
  {
    return Subhalos[index].TrackId;
  }
  const HBTxyz &GetComovingPosition(HBTInt index) const
  {
    return Subhalos[index].ComovingMostBoundPosition;
  }
  const HBTxyz GetPhysicalVelocity(HBTInt index) const
  {
    return Subhalos[index].PhysicalAverageVelocity;
  }
  HBTReal GetMass(HBTInt index) const
  {
    return Subhalos[index].GetMass();
  }
};

inline HBTInt GetCoreSize(HBTInt nbound)
/* get the size of the core that determines the position of the subhalo.
 * coresize controlled by SubCoreSizeFactor and SubCoreSizeMin.
 * if you do not want a cored center, then
 * set SubCoreSizeFactor=0 and SubCoreSizeMin=1 to use most-bound particle;
 * set SubCoreSizeFactor=1 to use all the particles*/
{
  int coresize = nbound * HBTConfig.SubCoreSizeFactor;
  if (coresize < HBTConfig.SubCoreSizeMin)
    coresize = HBTConfig.SubCoreSizeMin;
  if (coresize > nbound)
    coresize = nbound;
  return coresize;
}

class TrackKeyList_t : public KeyList_t<HBTInt, HBTInt>
{
  typedef HBTInt Index_t;
  typedef HBTInt Key_t;
  const SubhaloSnapshot_t &Snap;

public:
  TrackKeyList_t(SubhaloSnapshot_t &snap) : Snap(snap){};
  Index_t size() const
  {
    return Snap.size();
  }
  Key_t GetKey(Index_t i) const
  {
    return Snap.GetId(i);
  }
  Index_t GetIndex(Index_t i) const
  {
    return i;
  }
};

#endif
