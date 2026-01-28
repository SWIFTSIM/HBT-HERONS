#include <assert.h>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <list>
#include <numeric>
#include <sstream>
#include <string>
#include <typeinfo>

#include "../../halo.h"
#include "../../mymath.h"
#include "../../hash_integers.h"
#include "../../particle_exchanger.h"
#include "../../halo_particle_iterator.h"

struct HaloFragment_t
{
  HBTInt HaloId;
  HBTInt NumberParticles;
  int OriginalRank;
  int OriginalOrder;
  int TargetRank;
};

static void create_MPI_HaloInfo_t(MPI_Datatype &dtype)
{
  HaloFragment_t p;
#define NumAttr 13
  MPI_Datatype oldtypes[NumAttr];
  int blockcounts[NumAttr];
  MPI_Aint offsets[NumAttr], origin, extent;

  MPI_Get_address(&p, &origin);
  MPI_Get_address((&p) + 1, &extent); // to get the extent of s
  extent -= origin;

  int i = 0;
#define RegisterAttr(x, type, count)                                                                                   \
  {                                                                                                                    \
    MPI_Get_address(&(p.x), offsets + i);                                                                              \
    offsets[i] -= origin;                                                                                              \
    oldtypes[i] = type;                                                                                                \
    blockcounts[i] = count;                                                                                            \
    i++;                                                                                                               \
  }
  RegisterAttr(HaloId, MPI_HBT_INT, 1);
  RegisterAttr(NumberParticles, MPI_HBT_INT, 1);
  RegisterAttr(OriginalRank, MPI_INT, 1);
  RegisterAttr(OriginalOrder, MPI_INT, 1);
  RegisterAttr(TargetRank, MPI_INT, 1);
#undef RegisterAttr
  assert(i <= NumAttr);

  MPI_Type_create_struct(i, blockcounts, offsets, oldtypes, &dtype);
  MPI_Type_create_resized(dtype, (MPI_Aint)0, extent, &dtype);
  MPI_Type_commit(&dtype);
#undef NumAttr
}

/* Populates receive and send count and displacement vectors with the values to
 * use within MPI.*/
static void ComputeMPICounts(MpiWorker_t &world, const std::vector<IdRank_t> &TargetRank,
  std::vector<int> &send_counts, std::vector<int> &send_displacement,
  std::vector<int> &recv_counts, std::vector<int> &recv_displacement)
{
  for (size_t i = 0; i < TargetRank.size(); i++)
    send_counts[TargetRank[i].Rank]++;

  MPI_Alltoall(send_counts.data(), 1, MPI_INT, recv_counts.data(), 1, MPI_INT, world.Communicator);
  CompileOffsets(send_counts, send_displacement);
  CompileOffsets(recv_counts, recv_displacement);
}

/* Populates receive and send count and displacement vectors with the values to
 * use within MPI.*/
 static void ComputeMPICounts(MpiWorker_t &world, const std::vector<HaloFragment_t> &HaloFragments,
  std::vector<int> &send_counts, std::vector<int> &send_displacement,
  std::vector<int> &recv_counts, std::vector<int> &recv_displacement)
{
  for (size_t i = 0; i < HaloFragments.size(); i++)
    send_counts[HaloFragments[i].TargetRank]++;

  MPI_Alltoall(send_counts.data(), 1, MPI_INT, recv_counts.data(), 1, MPI_INT, world.Communicator);
  CompileOffsets(send_counts, send_displacement);
  CompileOffsets(recv_counts, recv_displacement);
}

inline bool CompHaloFragment_HaloId(const HaloFragment_t &a, const HaloFragment_t &b)
{
  return a.HaloId < b.HaloId;
}

inline bool CompHaloFragment_Size(const HaloFragment_t &a, const HaloFragment_t &b)
{
  return a.NumberParticles > b.NumberParticles;
}

inline bool CompHaloFragment_OriginalOrder(const HaloFragment_t &a, const HaloFragment_t &b)
{
  return a.OriginalOrder < b.OriginalOrder;
}

inline bool CompHaloInfo_TargetRank(const HaloFragment_t &a, const HaloFragment_t &b)
{
  return a.TargetRank < b.TargetRank;
}

inline bool CompHaloId(const Halo_t &a, const Halo_t &b)
{
  return a.HaloId < b.HaloId;
}

/* Gets the total size of haloes. This function returns a vector in which each
 * halo is only present once across all ranks. */
std::vector<HaloFragment_t> GetTotalHaloSizes(MpiWorker_t &world, const std::vector<Halo_t> &Halos)
{
  /* We store the properties of local halo fragments. */
  std::vector<HaloFragment_t> DisjointHaloFragments(Halos.size());
  for (size_t fragment_i = 0; fragment_i < DisjointHaloFragments.size(); fragment_i++)
  {
    /* To get global size of each FoF. */
    DisjointHaloFragments[fragment_i].HaloId = Halos[fragment_i].HaloId;
    DisjointHaloFragments[fragment_i].NumberParticles = Halos[fragment_i].Particles.size();

    /* For MPI communications. */
    DisjointHaloFragments[fragment_i].OriginalRank = world.rank();
    DisjointHaloFragments[fragment_i].TargetRank = RankFromIdHash(Halos[fragment_i].HaloId, world.size());
  }

  /* Sort vectors according to rank. */
  std::sort(DisjointHaloFragments.begin(), DisjointHaloFragments.end(), CompHaloInfo_TargetRank);

  /* We obtain counts and offsets for MPI communication. */
  std::vector<int> SendHaloCounts(world.size()),
                   RecvHaloCounts(world.size()),
                   SendHaloDisps (world.size()),
                   RecvHaloDisps (world.size());
  ComputeMPICounts(world, DisjointHaloFragments, SendHaloCounts, SendHaloDisps, RecvHaloCounts, RecvHaloDisps);

  /* Collect fragments within their assigned ranks. */
  size_t TotalRecvCount = std::accumulate(RecvHaloCounts.begin(), RecvHaloCounts.end(), 0);
  std::vector<HaloFragment_t> MergedHaloFragments(TotalRecvCount);
  MPI_Datatype MPI_HaloFragment_t;
  create_MPI_HaloInfo_t(MPI_HaloFragment_t);
  MPI_Alltoallv(DisjointHaloFragments.data()   , SendHaloCounts.data(), SendHaloDisps.data(), MPI_HaloFragment_t,
                MergedHaloFragments.data(), RecvHaloCounts.data(), RecvHaloDisps.data(), MPI_HaloFragment_t, world.Communicator);
  MPI_Type_free(&MPI_HaloFragment_t);

  /* Accumulate the size of unique haloes in this rank. */
  std::unordered_map<HBTInt, HBTInt> HaloSizes;
  for (size_t fragment_i = 0; fragment_i < MergedHaloFragments.size(); fragment_i++)
    HaloSizes[MergedHaloFragments[fragment_i].HaloId] += MergedHaloFragments[fragment_i].NumberParticles;

  /* Assign total sizes to the halo fragments */
  MergedHaloFragments.resize(HaloSizes.size());
  size_t fragment_i = 0;
  for (auto &unique_halo_it: HaloSizes)
  {
    MergedHaloFragments[fragment_i].HaloId = unique_halo_it.first;
    MergedHaloFragments[fragment_i].NumberParticles = unique_halo_it.second;
    fragment_i++;
  }

  return MergedHaloFragments;
}

/* Assigns each unique halo to a MPI rank based on a modified rounb-robin scheme.
 * The assignment tries to keep the total number of FoF particles in a rank to
 * be close to the average expected number of FoF particles per rank. Nonetheless,
 * the largest haloes will always be placed in a rank, regardless of whether
 * it violates the above criterion. */
std::vector<IdRank_t> RoundRobinAssignment(MpiWorker_t &world, std::vector<HaloFragment_t> &LocalHaloSizes)
{
  /* NOTE: this function might become inefficient for large numbers of haloes,
   * because one rank is tasked with holding all information instead of being
   * spread across all ranks that we have. We could use batch sort to prevent
   * this (not implemented in HBT-HERONS yet). */

  /* How many unique FoF groups do we have across all ranks. */
  std::vector<HBTInt> GlobalNumHalos(world.size(), 0);
  GlobalNumHalos[world.rank()] = LocalHaloSizes.size();
  MPI_Allreduce(MPI_IN_PLACE, GlobalNumHalos.data(), world.size(), MPI_HBT_INT, MPI_SUM, world.Communicator);

  /* Create a list of offsets for receiving data in the root rank */
  HBTInt NHalosTotal = GlobalNumHalos[0];
  std::vector<int> vector_offset(world.size(), 0);
  for(int i = 1; i < world.size(); i++)
  {
    vector_offset[i] = vector_offset[i-1] + GlobalNumHalos[i-1];
    NHalosTotal += GlobalNumHalos[i];
  }

  /* Add warnings because we can only pass INT as count and offset dtype to
   * MPI_Gatherv. If we have more new subhaloes than INT_MAX per rank or across
   * all ranks, an overflow will happen. */
  if (GlobalNumHalos[world.rank()] > INT_MAX)
   throw runtime_error("Error: in RoundRobinAssignment(). The number of new halos in a single rank is larger than INT_MAX. "
                       "It will cause required MPI communications to overflow. Please try more MPI threads. Aborting.\n");
  if (NHalosTotal > INT_MAX)
    throw runtime_error("Error: in RoundRobinAssignment(). The number of halos across all ranks is larger than INT_MAX. "
                        "It will cause required MPI communications to overflow. Please try more MPI threads. Aborting.\n");

  /* Convert count and displacement vectors into int. */
  std::vector<int> recv_counts(GlobalNumHalos.begin(), GlobalNumHalos.end());
  std::vector<int> recv_disps(vector_offset.begin(), vector_offset.end());

  /* To hold all the sizes of all haloes in the simulation. We only resize in
   * root rank because that is where the ranking will take place. */
  std::vector<HaloFragment_t> GlobalHaloSizes;
  if(world.rank() == 0)
    GlobalHaloSizes.resize(NHalosTotal);

  /* Gather all into root rank*/
  MPI_Datatype MPI_HaloFragment_t;
  create_MPI_HaloInfo_t(MPI_HaloFragment_t);
  MPI_Gatherv(LocalHaloSizes.data() , LocalHaloSizes.size(), MPI_HaloFragment_t,
              GlobalHaloSizes.data(), recv_counts.data(), recv_disps.data(), MPI_HaloFragment_t,
              0, world.Communicator);

  /* We can now do the assignment */
  if (world.rank() == 0)
  {
    int NumProc = world.size();

    /* Store original order and then sort in descending halo size (we start
     * assigning ranks to the largest haloes first) */
    HBTInt TotalParticles = 0;
    for(size_t i = 0; i < GlobalHaloSizes.size(); i++)
    {
      GlobalHaloSizes[i].OriginalOrder = i;
      TotalParticles += GlobalHaloSizes[i].NumberParticles;
    }

    std::sort(GlobalHaloSizes.begin(), GlobalHaloSizes.end(), CompHaloFragment_Size);

    /* We will try to assign haloes such that the number of particles per MPI
     * rank remains below or close to this value (ignoring the largest FoFs). */
    HBTInt MaxPartPerRank = TotalParticles / NumProc;
    int ThisRank = 0;
    std::vector<HBTInt> ParticlesOnRank(NumProc, 0);
    for (size_t i = 0; i < GlobalHaloSizes.size(); i++)
    {
      /* We do a round-robin over the ranks. However, if the next rank in the round-robin is
       * full (MaxPartPerRank) we skip to the next one. If a rank is empty MaxPartPerRank is
       * ignored so that we can always assign at least one halo per rank so that the biggest
       * halos can be assigned. */
      for (int TriedRanks = 0; TriedRanks < NumProc; TriedRanks++)
      {
        if((ParticlesOnRank[ThisRank] + GlobalHaloSizes[i].NumberParticles < MaxPartPerRank) \
        || (ParticlesOnRank[ThisRank] == 0))
        {
          GlobalHaloSizes[i].TargetRank = ThisRank;
          ParticlesOnRank[ThisRank] += GlobalHaloSizes[i].NumberParticles;
          ThisRank++; // Next group will start with next rank.
          ThisRank %= NumProc;
          break;
        }
        ThisRank++;
        ThisRank %= NumProc;

        if(TriedRanks == NumProc)
        {
          // Tried all ranks for this halo and couldn't find a place - abort.
          // This should be impossible.
          stringstream error_message;
          error_message << "Failed to distribute group with ID " << GlobalHaloSizes[i].HaloId << "." << endl;
          error_message << "It has " << GlobalHaloSizes[i].NumberParticles << " particles." << endl;
          error_message << "Max particles per rank (summed over assigned groups) is set to " << MaxPartPerRank << "." << endl;
          error_message << "When trying to assign this group, particle load on ranks so far is:" << endl;
          for(int rank = 0; rank < NumProc; rank++)
            error_message << "  Rank " << rank << ": " << ParticlesOnRank[rank] << endl;

          throw runtime_error(error_message.str());
        }
      }
    }

    /* Return to original order to scatter back to original tasks. */
    std::sort(GlobalHaloSizes.begin(), GlobalHaloSizes.end(), CompHaloFragment_OriginalOrder);
  }

  MPI_Scatterv(GlobalHaloSizes.data(), recv_counts.data()   , recv_disps.data() , MPI_HaloFragment_t,
               LocalHaloSizes.data() , LocalHaloSizes.size(), MPI_HaloFragment_t, 0, world.Communicator);
  MPI_Type_free(&MPI_HaloFragment_t);

  std::vector<IdRank_t> TargetTask(LocalHaloSizes.size());
  for(size_t i = 0; i < LocalHaloSizes.size(); i++)
  {
    TargetTask[i].Id = LocalHaloSizes[i].HaloId;
    TargetTask[i].Rank = LocalHaloSizes[i].TargetRank;
  }
  return TargetTask;
}

/* This function makes the round-robin-based task assignment, which is done on a
 * halo-level after fragments have been merged, into the same size and ordering
 * as haloes. */
 std::vector<IdRank_t> CommunicateTaskAssignment(MpiWorker_t &world, const std::vector<IdRank_t> HaloTaskAssignment, const std::vector<Halo_t> &Halos)
{
  /* We store the properties of local halo fragments. */
  std::vector<HaloFragment_t> DisjointHaloFragments(Halos.size());
  for (size_t fragment_i = 0; fragment_i < DisjointHaloFragments.size(); fragment_i++)
  {
    DisjointHaloFragments[fragment_i].HaloId = Halos[fragment_i].HaloId;
    DisjointHaloFragments[fragment_i].TargetRank = RankFromIdHash(Halos[fragment_i].HaloId, world.size());
  }

  /* Sort vectors according to rank. */
  std::sort(DisjointHaloFragments.begin(), DisjointHaloFragments.end(), CompHaloInfo_TargetRank);

  /* We obtain counts and offsets for MPI communication. */
  std::vector<int> SendHaloCounts(world.size()),
                   RecvHaloCounts(world.size()),
                   SendHaloDisps (world.size()),
                   RecvHaloDisps (world.size());
  ComputeMPICounts(world, DisjointHaloFragments, SendHaloCounts, SendHaloDisps, RecvHaloCounts, RecvHaloDisps);

  /* Collect fragments within their assigned ranks. */
  size_t TotalRecvCount = std::accumulate(RecvHaloCounts.begin(), RecvHaloCounts.end(), 0);
  std::vector<HaloFragment_t> ProbingHaloFragments(TotalRecvCount);
  MPI_Datatype MPI_HaloFragment_t;
  create_MPI_HaloInfo_t(MPI_HaloFragment_t);
  MPI_Alltoallv(DisjointHaloFragments.data()   , SendHaloCounts.data(), SendHaloDisps.data(), MPI_HaloFragment_t,
                ProbingHaloFragments.data(), RecvHaloCounts.data(), RecvHaloDisps.data(), MPI_HaloFragment_t, world.Communicator);

  /* Store original order and sort by host halo ID */
  for(size_t i = 0; i < ProbingHaloFragments.size(); i++)
    ProbingHaloFragments[i].OriginalOrder = i;

  std::sort(ProbingHaloFragments.begin(), ProbingHaloFragments.end(), CompHaloFragment_HaloId);

  /* We iterate across all halo fragments and copy the assignment */
  std::vector<IdRank_t>::const_iterator CurrentHaloRankAssignment = HaloTaskAssignment.begin();
  for(size_t i = 0; i < ProbingHaloFragments.size(); i++)
  {
    while ((CurrentHaloRankAssignment->Id != ProbingHaloFragments[i].HaloId) \
          &(CurrentHaloRankAssignment != HaloTaskAssignment.end()))
      CurrentHaloRankAssignment++;

    if(ProbingHaloFragments[i].HaloId == CurrentHaloRankAssignment->Id)
      ProbingHaloFragments[i].TargetRank = CurrentHaloRankAssignment->Rank;
  }

  /* Return to the original order and then to original task. */
  std::sort(ProbingHaloFragments.begin(), ProbingHaloFragments.end(), CompHaloFragment_OriginalOrder);
  MPI_Alltoallv(ProbingHaloFragments.data() , RecvHaloCounts.data(), RecvHaloDisps.data(), MPI_HaloFragment_t,
                DisjointHaloFragments.data(), SendHaloCounts.data(), SendHaloDisps.data(), MPI_HaloFragment_t, world.Communicator);
  MPI_Type_free(&MPI_HaloFragment_t);

  std::vector<IdRank_t> RankAssignment(DisjointHaloFragments.size());
  for(size_t i = 0; i < DisjointHaloFragments.size(); i++)
  {
    RankAssignment[i].Id = DisjointHaloFragments[i].HaloId;
    RankAssignment[i].Rank = DisjointHaloFragments[i].TargetRank;
  }
  return RankAssignment;
}

/* Assigns MPI tasks to haloes so that the number of FoF particles is approximately
 * equal across ranks. */
static std::vector<IdRank_t> DecideTargetProcessor(MpiWorker_t &world, std::vector<Halo_t> &Halos)
{
  /* First we compute total FoF group sizes. Note that the ordering and vector
   * size of HaloSizes and Halos is not the same, as Halos contains disjoint
   * FoF fragments but HaloSizes has merged all of them (to get total size). */
  std::vector<HaloFragment_t> HaloSizes = GetTotalHaloSizes(world, Halos);

  /* Assign a target MPI rank to haloes based on their total size. The resulting
   * vector has the same global ordering as HaloSizes. */
  std::vector<IdRank_t> HaloTaskAssignment = RoundRobinAssignment(world, HaloSizes);

  /* We now need to propagate the rank assignment into a vector with the same
   * ordering as the disjoint halo fragments (Halos). */
  std::vector<IdRank_t> HaloFragmentTaskAssignment = CommunicateTaskAssignment(world, HaloTaskAssignment, Halos);

  return HaloFragmentTaskAssignment;
}

/* After communicating disjoint pieces of FoF groups from different MPI ranks,
 * merge them into a single FoF group. */
static void MergeHaloFragments(vector<Halo_t> &Halos)
{
  if (Halos.empty())
    return;

  std::sort(Halos.begin(), Halos.end(), CompHaloId);

  auto it1 = Halos.begin();
  for (auto it2 = it1 + 1; it2 != Halos.end(); ++it2)
  {
    if (it2->HaloId == it1->HaloId) // Piece of the same FoF, we can merge them.
    {
      it1->Particles.insert(it1->Particles.end(), it2->Particles.begin(), it2->Particles.end());
    }
    else // This piece is the start of a different FoF .
    {
      ++it1;
      if (it2 != it1)
        *it1 = move(*it2);
    }
  }
  Halos.resize(it1 - Halos.begin() + 1);

  for (auto &h : Halos)
  {
    h.AverageCoordinates();

    /* Sort particles by ID to get reproducible Particle vectors regardless of
     * the number of MPI ranks that we use. Important because Halos.Particles is
     * swapped with the Subhalo_t.Particles vector, which will lead to unbinding
     * differences if we do subsampling and change number of MPI ranks. */
     std::sort(h.Particles.begin(), h.Particles.end(), ParticleExchangeComp::CompParticleId);
  }
}

/* Communicates the basic properties of each FoF group segment to its assigned rank. */
void ExchangeHaloProperties(MpiWorker_t &world, const std::vector<IdRank_t> &TargetRank,
                            std::vector<Halo_t> &InHalosSorted, std::vector<Halo_t> &OutHalos,
                            const std::vector<int> &SendHaloDisps, const std::vector<int> &SendHaloCounts,
                            const std::vector<int> &RecvHaloDisps, const std::vector<int> &RecvHaloCounts)
{
  MPI_Datatype MPI_Halo_Shell_t;
  create_MPI_Halo_Id_type(MPI_Halo_Shell_t);
  MPI_Alltoallv(InHalosSorted.data(), SendHaloCounts.data(), SendHaloDisps.data(), MPI_Halo_Shell_t,
                OutHalos.data()     , RecvHaloCounts.data(), RecvHaloDisps.data(), MPI_Halo_Shell_t, world.Communicator);
  MPI_Type_free(&MPI_Halo_Shell_t);
}

/* Communicates the particles of each FoF group segment to its assigned rank. */
void ExchangeHaloParticles(MpiWorker_t &world, std::vector<IdRank_t> &TargetRank,
                           std::vector<Halo_t> &InHalosSorted, std::vector<Halo_t> &OutHalos,
                           std::vector<int> &SendHaloDisps, std::vector<int> &SendHaloCounts,
                           std::vector<int> &RecvHaloDisps, std::vector<int> &RecvHaloCounts)
{
  /* Measure required space for particle vector of each FoF fragment.  */
  std::vector<HBTInt> InHalosSortedSizes(InHalosSorted.size());
  for (size_t halo_i = 0; halo_i < InHalosSorted.size(); halo_i++)
    InHalosSortedSizes[halo_i] = InHalosSorted[halo_i].Particles.size();

  /* Communicate the required particle vector sizes.  */
  vector<HBTInt> OutHaloSizes(OutHalos.size());
  MPI_Alltoallv(InHalosSortedSizes.data(), SendHaloCounts.data(), SendHaloDisps.data(), MPI_HBT_INT,
                OutHaloSizes.data()      , RecvHaloCounts.data(), RecvHaloDisps.data(), MPI_HBT_INT, world.Communicator);

  /* Resize the particle vector of each FoF fragment so it can hold all required
   * particles. */
  for (size_t halo_i = 0; halo_i < OutHalos.size(); halo_i++)
    OutHalos[halo_i].Particles.resize(OutHaloSizes[halo_i]);

  /* Combined iterator for groups of haloes */
  typedef typename vector<Halo_t>::iterator HaloIterator_t;
  typedef HaloParticleIterator_t<HaloIterator_t> ParticleIterator_t;

  vector<ParticleIterator_t> InParticleIterator(world.size());
  vector<ParticleIterator_t> OutParticleIterator(world.size());
  for (int rank = 0; rank < world.size(); rank++)
  {
    InParticleIterator[rank] .init(InHalosSorted.begin() + SendHaloDisps[rank],
                                   InHalosSorted.begin() + SendHaloDisps[rank] + SendHaloCounts[rank]);
    OutParticleIterator[rank].init(OutHalos.begin() + RecvHaloDisps[rank],
                                   OutHalos.begin() + RecvHaloDisps[rank] + RecvHaloCounts[rank]);
  }

  /* Total particle count in each rank iterator. */
  vector<HBTInt> InParticleCount(world.size(), 0);
  for (HBTInt i = 0; i < InHalosSorted.size(); i++)
    InParticleCount[TargetRank[i].Rank] += InHalosSorted[i].Particles.size();

  /* Distribute halo particles */
  MPI_Datatype MPI_HBT_Particle;
  Particle_t().create_MPI_type(MPI_HBT_Particle);
  MyAllToAll<Particle_t, ParticleIterator_t, ParticleIterator_t>(world, InParticleIterator, InParticleCount,
                                                                 OutParticleIterator, MPI_HBT_Particle);
  MPI_Type_free(&MPI_HBT_Particle);
}

static void ExchangeHaloFragments(MpiWorker_t &world, std::vector<Halo_t> &InHalos)
{
  /* Decide how haloes will be distributed across MPI ranks. This vector is
   * sorted in ascending TargetRank order. */
  std::vector<IdRank_t> TargetRank = DecideTargetProcessor(world, InHalos);

  /* Compute counts and offsets for MPI communication */
  std::vector<int> SendHaloCounts(world.size()),
                   RecvHaloCounts(world.size()),
                   SendHaloDisps (world.size()),
                   RecvHaloDisps (world.size());
  ComputeMPICounts(world, TargetRank, SendHaloCounts, SendHaloDisps, RecvHaloCounts, RecvHaloDisps);

  /* Group local FoF fragments based on which MPI rank they will be sent to. */
  std::vector<Halo_t> InHalosSorted(InHalos.size());
  for (size_t halo_i = 0; halo_i < InHalos.size(); halo_i++)
    InHalosSorted[halo_i] = std::move(InHalos[TargetRank[halo_i].Id]);

  /* To hold incoming halo fragments. */
  size_t NumberIncomingHaloFragments = std::accumulate(RecvHaloCounts.begin(), RecvHaloCounts.end(), 0);
  std::vector<Halo_t> OutHalos(NumberIncomingHaloFragments);

  /* First we send the basic information of each FoF group, then the particle
   * information. These function calls handle vector re-sizing. */
  ExchangeHaloProperties(world, TargetRank, InHalosSorted, OutHalos, SendHaloCounts, SendHaloDisps, RecvHaloCounts, RecvHaloDisps);
  ExchangeHaloParticles (world, TargetRank, InHalosSorted, OutHalos, SendHaloCounts, SendHaloDisps, RecvHaloCounts, RecvHaloDisps);

  /* Done. */
  InHalos.swap(OutHalos);
}

/* Gathers fragments of FoF groups that have been read by different MPI ranks into
 * the same MPI rank. Once collected, each FoF fragment is merged together to get the
 * complete FoF group. */
void CollectHaloFragments(MpiWorker_t &world, vector<Halo_t> &Halos)
{
  ExchangeHaloFragments(world, Halos);
  MergeHaloFragments(Halos);
}