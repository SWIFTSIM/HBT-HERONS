#include <assert.h>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <list>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <typeinfo>
#include <unordered_map>

#include "../../halo.h"
#include "../../mymath.h"
#include "../../hash_integers.h"
#include "../../particle_exchanger.h"
#include "../../halo_particle_iterator.h"

struct HaloFragment_t
{
  HBTInt HaloId;
  HBTInt NumberParticles;
  double TotalMass;
  HBTxyz ComovingAveragePosition;
  int OriginalRank;
  int OriginalOrder;
  int TargetRank;
};

static void create_MPI_HaloInfo_t(MPI_Datatype &dtype)
{
  HaloFragment_t p;
#define NumAttr 7
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
  RegisterAttr(TotalMass, MPI_DOUBLE, 1);
  RegisterAttr(ComovingAveragePosition[0], MPI_HBT_REAL, 3);
  RegisterAttr(OriginalRank, MPI_INT, 1);
  RegisterAttr(OriginalOrder, MPI_INT, 1);
  RegisterAttr(TargetRank, MPI_INT, 1);
#undef RegisterAttr
  assert(i == NumAttr);

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

/* Gets the total particle count and centre of mass of each FoF group across all
 * ranks. Each FoF group may be split across multiple ranks as disjoint fragments;
 * this function collects and merges those fragments. Returns a vector where each
 * FoF group appears exactly once, with its total particle count and COM. */
static std::vector<HaloFragment_t> GetTotalHaloSizes(MpiWorker_t &world, const std::vector<Halo_t> &Halos)
{
  /* Compute the COM and particle count of each local FoF fragment. */
  std::vector<HaloFragment_t> DisjointHaloFragments(Halos.size());
  for (size_t fragment_i = 0; fragment_i < DisjointHaloFragments.size(); fragment_i++)
  {
    DisjointHaloFragments[fragment_i].HaloId = Halos[fragment_i].HaloId;
    DisjointHaloFragments[fragment_i].NumberParticles = Halos[fragment_i].Particles.size();
    DisjointHaloFragments[fragment_i].TotalMass =
      AveragePosition(DisjointHaloFragments[fragment_i].ComovingAveragePosition,
                      Halos[fragment_i].Particles.data(),
                      Halos[fragment_i].Particles.size());

    /* Route each fragment to a home rank via HaloId hash for accumulation. */
    DisjointHaloFragments[fragment_i].OriginalRank = world.rank();
    DisjointHaloFragments[fragment_i].TargetRank = RankFromIdHash(Halos[fragment_i].HaloId, world.size());
  }

  /* Sort by target rank for MPI communication. */
  std::sort(DisjointHaloFragments.begin(), DisjointHaloFragments.end(), CompHaloInfo_TargetRank);

  std::vector<int> SendHaloCounts(world.size()),
                   RecvHaloCounts(world.size()),
                   SendHaloDisps (world.size()),
                   RecvHaloDisps (world.size());
  ComputeMPICounts(world, DisjointHaloFragments, SendHaloCounts, SendHaloDisps, RecvHaloCounts, RecvHaloDisps);

  size_t TotalRecvCount = std::accumulate(RecvHaloCounts.begin(), RecvHaloCounts.end(), 0);
  std::vector<HaloFragment_t> MergedHaloFragments(TotalRecvCount);
  MPI_Datatype MPI_HaloFragment_t;
  create_MPI_HaloInfo_t(MPI_HaloFragment_t);
  MPI_Alltoallv(DisjointHaloFragments.data(), SendHaloCounts.data(), SendHaloDisps.data(), MPI_HaloFragment_t,
                MergedHaloFragments.data(), RecvHaloCounts.data(), RecvHaloDisps.data(), MPI_HaloFragment_t,
                world.Communicator);
  MPI_Type_free(&MPI_HaloFragment_t);

  /* Accumulate particle counts and particle-count-weighted COMs for unique FoFs.
   * For the COM accumulation, offsets relative to the first fragment's position
   * are used to handle periodic boundary conditions. */
  struct HaloAccumulator
  {
    HBTInt TotalParticles = 0;
    double TotalMass = 0.;
    double COM[3] = {0., 0., 0.};
    double Origin[3] = {0., 0., 0.};
    bool IsFirst = true;
  };
  std::unordered_map<HBTInt, HaloAccumulator> Accumulators;

  for (auto &frag : MergedHaloFragments)
  {
    auto &acc = Accumulators[frag.HaloId];
    if (acc.IsFirst)
    {
      for (int j = 0; j < 3; j++)
        acc.Origin[j] = frag.ComovingAveragePosition[j];
      acc.IsFirst = false;
    }
    double weight = frag.TotalMass;
    acc.TotalParticles += frag.NumberParticles;
    acc.TotalMass += frag.TotalMass;
    for (int j = 0; j < 3; j++)
    {
      double offset = frag.ComovingAveragePosition[j] - acc.Origin[j];
      if (HBTConfig.PeriodicBoundaryOn)
        offset = NEAREST(offset);
      acc.COM[j] += weight * offset;
    }
  }

  /* Build output: one entry per unique FoF with total size and final COM. */
  MergedHaloFragments.resize(Accumulators.size());
  size_t out_i = 0;
  for (auto &kv : Accumulators)
  {
    HaloAccumulator &acc = kv.second;
    MergedHaloFragments[out_i].HaloId = kv.first;
    MergedHaloFragments[out_i].NumberParticles = acc.TotalParticles;
    MergedHaloFragments[out_i].TotalMass = acc.TotalMass;
    for (int j = 0; j < 3; j++)
    {
      double com = acc.COM[j] / acc.TotalMass + acc.Origin[j];
      if (HBTConfig.PeriodicBoundaryOn)
        com = position_modulus(com, HBTConfig.BoxSize);
      MergedHaloFragments[out_i].ComovingAveragePosition[j] = com;
    }
    out_i++;
  }

  return MergedHaloFragments;
}

/* Assigns each FoF group to an MPI rank. Starts a the spatial cell-based
 * assignment (preserving locality and inter-snapshot stability),
 * then greedily moves FoFs from the most overloaded rank to the most underloaded
 * rank until all ranks are within a tolerance of the mean particle count. */
static std::vector<IdRank_t> SpatialAssignmentWithRebalancing(MpiWorker_t &world,
                                                               std::vector<HaloFragment_t> &LocalHaloSizes)
{
  /* How many unique FoF groups do we have across all ranks. */
  std::vector<HBTInt> GlobalNumHalos(world.size(), 0);
  GlobalNumHalos[world.rank()] = LocalHaloSizes.size();
  MPI_Allreduce(MPI_IN_PLACE, GlobalNumHalos.data(), world.size(), MPI_HBT_INT, MPI_SUM, world.Communicator);

  for (int r = 0; r < world.size(); r++)
    if (GlobalNumHalos[r] > INT_MAX)
      throw std::runtime_error("Error: in SpatialAssignmentWithRebalancing(). The number of halos in a single rank "
                               "is larger than INT_MAX, causing MPI overflow. Please use more MPI ranks. Aborting.\n");

  HBTInt NHalosTotal = GlobalNumHalos[0];
  std::vector<int> vector_offset(world.size(), 0);
  for(int i = 1; i < world.size(); i++)
  {
    vector_offset[i] = vector_offset[i-1] + GlobalNumHalos[i-1];
    NHalosTotal += GlobalNumHalos[i];
  }

  if (NHalosTotal > INT_MAX)
    throw std::runtime_error("Error: in SpatialAssignmentWithRebalancing(). The total number of halos "
                             "is larger than INT_MAX, causing MPI overflow. Please use more MPI ranks. Aborting.\n");

  std::vector<int> recv_counts(GlobalNumHalos.begin(), GlobalNumHalos.end());
  std::vector<int> recv_disps(vector_offset.begin(), vector_offset.end());

  std::vector<HaloFragment_t> GlobalHaloSizes;
  if(world.rank() == 0)
    GlobalHaloSizes.resize(NHalosTotal);

  MPI_Datatype MPI_HaloFragment_t;
  create_MPI_HaloInfo_t(MPI_HaloFragment_t);
  MPI_Gatherv(LocalHaloSizes.data() , LocalHaloSizes.size(), MPI_HaloFragment_t,
              GlobalHaloSizes.data(), recv_counts.data(), recv_disps.data(), MPI_HaloFragment_t,
              0, world.Communicator);

  if (world.rank() == 0)
  {
    int NumProc = world.size();

    /* Record original order before any sorting, so we can scatter back correctly. */
    for(size_t i = 0; i < GlobalHaloSizes.size(); i++)
      GlobalHaloSizes[i].OriginalOrder = i;

    /* Initial assignment: map each FoF's COM to a spatial cell. */
    auto dims = ClosestFactors(NumProc, 3);
    HBTxyz step;
    for (int j = 0; j < 3; j++)
      step[j] = HBTConfig.BoxSize / dims[j];

    std::vector<HBTInt> ParticlesOnRank(NumProc, 0);
    for (size_t i = 0; i < GlobalHaloSizes.size(); i++)
    {
      GlobalHaloSizes[i].TargetRank = AssignCell(GlobalHaloSizes[i].ComovingAveragePosition, step, dims);
      ParticlesOnRank[GlobalHaloSizes[i].TargetRank] += GlobalHaloSizes[i].NumberParticles;
    }

    /* Greedy rebalancing: repeatedly move the largest FoF that fits from the most
     * overloaded rank to the most underloaded rank. A FoF of size S "fits" when
     * 2*S < gap (= load[R_max] - load[R_min]), which ensures the move reduces the
     * max load without making R_min the new maximum. Stop when all ranks are within
     * the tolerance of the mean, or when no suitable FoF exists. */
    HBTInt TotalParticles = 0;
    for (int r = 0; r < NumProc; r++)
      TotalParticles += ParticlesOnRank[r];
    HBTInt MeanParticles = TotalParticles / NumProc;
    const double Tolerance = 0.05;

    /* Build per-rank lists of (NumberParticles, global_index), sorted descending. */
    std::vector<std::vector<std::pair<HBTInt, int>>> RankHalos(NumProc);
    for (int i = 0; i < (int)GlobalHaloSizes.size(); i++)
      RankHalos[GlobalHaloSizes[i].TargetRank].emplace_back(GlobalHaloSizes[i].NumberParticles, i);
    for (int r = 0; r < NumProc; r++)
      std::sort(RankHalos[r].begin(), RankHalos[r].end(), std::greater<std::pair<HBTInt, int>>());

    while (true)
    {
      int R_max = std::max_element(ParticlesOnRank.begin(), ParticlesOnRank.end()) - ParticlesOnRank.begin();
      int R_min = std::min_element(ParticlesOnRank.begin(), ParticlesOnRank.end()) - ParticlesOnRank.begin();

      if (ParticlesOnRank[R_max] <= (HBTInt)((1.0 + Tolerance) * MeanParticles))
        break;

      HBTInt gap = ParticlesOnRank[R_max] - ParticlesOnRank[R_min];

      /* Find the largest FoF on R_max satisfying 2*size < gap. */
      auto &list = RankHalos[R_max];
      auto it = list.begin();
      while (it != list.end() && 2 * it->first >= gap)
        ++it;

      if (it == list.end())
        break;

      int halo_idx = it->second;
      HBTInt halo_size = it->first;

      GlobalHaloSizes[halo_idx].TargetRank = R_min;
      ParticlesOnRank[R_max] -= halo_size;
      ParticlesOnRank[R_min] += halo_size;

      list.erase(it);
      auto insert_pos = std::lower_bound(
        RankHalos[R_min].begin(), RankHalos[R_min].end(),
        std::make_pair(halo_size, halo_idx),
        std::greater<std::pair<HBTInt, int>>()
      );
      RankHalos[R_min].insert(insert_pos, {halo_size, halo_idx});
    }

    /* Restore original order before scattering back. */
    std::sort(GlobalHaloSizes.begin(), GlobalHaloSizes.end(), CompHaloFragment_OriginalOrder);
  }

  MPI_Scatterv(GlobalHaloSizes.data(), recv_counts.data(), recv_disps.data() , MPI_HaloFragment_t,
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

/* Propagates the per-FoF rank assignment (computed on merged unique FoFs) back to
 * the disjoint fragment level, returning a vector indexed in the same order as
 * Halos where .Id is the array index and .Rank is the destination rank. */
static std::vector<IdRank_t> CommunicateTaskAssignment(MpiWorker_t &world,
                                                        const std::vector<IdRank_t> &HaloTaskAssignment,
                                                        const std::vector<Halo_t> &Halos)
{
  /* Build fragment descriptors in the same order as Halos. Record the original
   * index in OriginalOrder before sorting so we can map back after the exchange. */
  std::vector<HaloFragment_t> DisjointHaloFragments(Halos.size());
  for (size_t fragment_i = 0; fragment_i < DisjointHaloFragments.size(); fragment_i++)
  {
    DisjointHaloFragments[fragment_i].HaloId = Halos[fragment_i].HaloId;
    DisjointHaloFragments[fragment_i].OriginalOrder = fragment_i;
    DisjointHaloFragments[fragment_i].TargetRank = RankFromIdHash(Halos[fragment_i].HaloId, world.size());
  }

  /* Sort by target rank for the forward alltoallv. */
  std::sort(DisjointHaloFragments.begin(), DisjointHaloFragments.end(), CompHaloInfo_TargetRank);

  std::vector<int> SendHaloCounts(world.size()),
                   RecvHaloCounts(world.size()),
                   SendHaloDisps (world.size()),
                   RecvHaloDisps (world.size());
  ComputeMPICounts(world, DisjointHaloFragments, SendHaloCounts, SendHaloDisps, RecvHaloCounts, RecvHaloDisps);

  size_t TotalRecvCount = std::accumulate(RecvHaloCounts.begin(), RecvHaloCounts.end(), 0);
  std::vector<HaloFragment_t> ProbingHaloFragments(TotalRecvCount);
  MPI_Datatype MPI_HaloFragment_t;
  create_MPI_HaloInfo_t(MPI_HaloFragment_t);
  MPI_Alltoallv(DisjointHaloFragments.data(), SendHaloCounts.data(), SendHaloDisps.data(), MPI_HaloFragment_t,
                ProbingHaloFragments.data(), RecvHaloCounts.data(), RecvHaloDisps.data(), MPI_HaloFragment_t,
                world.Communicator);

  /* Build a map from HaloId to target rank for O(1) lookup. */
  std::unordered_map<HBTInt, int> HaloIdToRank;
  HaloIdToRank.reserve(HaloTaskAssignment.size());
  for (const auto &assignment : HaloTaskAssignment)
    HaloIdToRank[assignment.Id] = assignment.Rank;

  for (auto &frag : ProbingHaloFragments)
    frag.TargetRank = HaloIdToRank.at(frag.HaloId);

  /* ProbingHaloFragments is still in received order so send it back directly. */
  MPI_Alltoallv(ProbingHaloFragments.data(), RecvHaloCounts.data(), RecvHaloDisps.data(), MPI_HaloFragment_t,
                DisjointHaloFragments.data(), SendHaloCounts.data(), SendHaloDisps.data(), MPI_HaloFragment_t,
                world.Communicator);
  MPI_Type_free(&MPI_HaloFragment_t);

  /* DisjointHaloFragments is now in the TargetRank-sorted order with updated
   * TargetRanks. Use OriginalOrder (the original Halos index) to map back. */
  std::vector<IdRank_t> RankAssignment(Halos.size());
  for(size_t i = 0; i < DisjointHaloFragments.size(); i++)
  {
    int orig_idx = DisjointHaloFragments[i].OriginalOrder;
    RankAssignment[orig_idx].Id = orig_idx;
    RankAssignment[orig_idx].Rank = DisjointHaloFragments[i].TargetRank;
  }
  return RankAssignment;
}

/* Assigns MPI ranks to FoF fragments, balancing particle count while preserving
 * spatial locality. */
static std::vector<IdRank_t> DecideTargetProcessor(MpiWorker_t &world, std::vector<Halo_t> &Halos)
{
  /* Collect total sizes and COMs for each unique FoF across all ranks. */
  std::vector<HaloFragment_t> HaloSizes = GetTotalHaloSizes(world, Halos);

  /* Assign a target rank to each unique FoF. */
  std::vector<IdRank_t> HaloTaskAssignment = SpatialAssignmentWithRebalancing(world, HaloSizes);

  /* Propagate the assignment back to the fragment level. */
  std::vector<IdRank_t> HaloFragmentTaskAssignment = CommunicateTaskAssignment(world, HaloTaskAssignment, Halos);

  return HaloFragmentTaskAssignment;
}

/* After communicating disjoint pieces of FoF groups from different MPI ranks,
 * merge them into a single FoF group. */
static void MergeHaloFragments(std::vector<Halo_t> &Halos)
{
  if (Halos.empty())
    return;

  std::sort(Halos.begin(), Halos.end(), CompHaloId);

  auto it1 = Halos.begin();
  for (auto it2 = it1 + 1; it2 != Halos.end(); ++it2)
  {
    if (it2->HaloId == it1->HaloId) // Piece of the same FoF, we can merge them
    {
      it1->Particles.insert(it1->Particles.end(), it2->Particles.begin(), it2->Particles.end());
    }
    else // This piece is the start of a different FoF
    {
      ++it1;
      if (it2 != it1)
        *it1 = std::move(*it2);
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

/* Communicates the basic properties of each FoF fragment to its assigned rank. */
static void ExchangeHaloProperties(MpiWorker_t &world,
                                   std::vector<Halo_t> &InHalosSorted, std::vector<Halo_t> &OutHalos,
                                   const std::vector<int> &SendHaloCounts, const std::vector<int> &SendHaloDisps,
                                   const std::vector<int> &RecvHaloCounts, const std::vector<int> &RecvHaloDisps)
{
  MPI_Datatype MPI_Halo_Shell_t;
  create_MPI_Halo_Id_type(MPI_Halo_Shell_t);
  MPI_Alltoallv(InHalosSorted.data(), SendHaloCounts.data(), SendHaloDisps.data(), MPI_Halo_Shell_t,
                OutHalos.data()     , RecvHaloCounts.data(), RecvHaloDisps.data(), MPI_Halo_Shell_t,
                world.Communicator);
  MPI_Type_free(&MPI_Halo_Shell_t);
}

/* Communicates the particles of each FoF fragment to its assigned rank. */
static void ExchangeHaloParticles(MpiWorker_t &world, const std::vector<IdRank_t> &TargetRank,
                                  std::vector<Halo_t> &InHalosSorted, std::vector<Halo_t> &OutHalos,
                                  const std::vector<int> &SendHaloCounts, const std::vector<int> &SendHaloDisps,
                                  const std::vector<int> &RecvHaloCounts, const std::vector<int> &RecvHaloDisps)
{
  /* Communicate the particle vector size of each fragment so receivers can resize. */
  std::vector<HBTInt> InHalosSortedSizes(InHalosSorted.size());
  for (size_t halo_i = 0; halo_i < InHalosSorted.size(); halo_i++)
    InHalosSortedSizes[halo_i] = InHalosSorted[halo_i].Particles.size();

  std::vector<HBTInt> OutHaloSizes(OutHalos.size());
  MPI_Alltoallv(InHalosSortedSizes.data(), SendHaloCounts.data(), SendHaloDisps.data(), MPI_HBT_INT,
                OutHaloSizes.data()      , RecvHaloCounts.data(), RecvHaloDisps.data(), MPI_HBT_INT,
                world.Communicator);

  for (size_t halo_i = 0; halo_i < OutHalos.size(); halo_i++)
    OutHalos[halo_i].Particles.resize(OutHaloSizes[halo_i]);

  typedef typename std::vector<Halo_t>::iterator HaloIterator_t;
  typedef HaloParticleIterator_t<HaloIterator_t> ParticleIterator_t;

  std::vector<ParticleIterator_t> InParticleIterator(world.size());
  std::vector<ParticleIterator_t> OutParticleIterator(world.size());
  for (int rank = 0; rank < world.size(); rank++)
  {
    InParticleIterator[rank] .init(InHalosSorted.begin() + SendHaloDisps[rank],
                                   InHalosSorted.begin() + SendHaloDisps[rank] + SendHaloCounts[rank]);
    OutParticleIterator[rank].init(OutHalos.begin() + RecvHaloDisps[rank],
                                   OutHalos.begin() + RecvHaloDisps[rank] + RecvHaloCounts[rank]);
  }

  /* Total particle count destined for each rank. */
  std::vector<HBTInt> InParticleCount(world.size(), 0);
  for (size_t i = 0; i < InHalosSorted.size(); i++)
    InParticleCount[TargetRank[i].Rank] += InHalosSorted[i].Particles.size();

  MPI_Datatype MPI_HBT_Particle;
  Particle_t().create_MPI_type(MPI_HBT_Particle);
  MyAllToAll<Particle_t, ParticleIterator_t, ParticleIterator_t>(world, InParticleIterator, InParticleCount,
                                                                 OutParticleIterator, MPI_HBT_Particle);
  MPI_Type_free(&MPI_HBT_Particle);
}

static void ExchangeHaloFragments(MpiWorker_t &world, std::vector<Halo_t> &InHalos)
{
  /* Decide how fragments will be distributed across MPI ranks. Returns a vector
   * indexed in the same order as InHalos, with .Id = array index and
   * .Rank = destination rank. */
  std::vector<IdRank_t> TargetRank = DecideTargetProcessor(world, InHalos);

  /* Sort by destination rank so InHalosSorted is grouped correctly for Alltoallv. */
  std::sort(TargetRank.begin(), TargetRank.end(), CompareRank);

  /* Reorder InHalos to match the sorted TargetRank. */
  std::vector<Halo_t> InHalosSorted(InHalos.size());
  for (size_t halo_i = 0; halo_i < InHalos.size(); halo_i++)
    InHalosSorted[halo_i] = std::move(InHalos[TargetRank[halo_i].Id]);

  std::vector<int> SendHaloCounts(world.size()),
                   RecvHaloCounts(world.size()),
                   SendHaloDisps (world.size()),
                   RecvHaloDisps (world.size());
  ComputeMPICounts(world, TargetRank, SendHaloCounts, SendHaloDisps, RecvHaloCounts, RecvHaloDisps);

  size_t NumberIncomingHaloFragments = std::accumulate(RecvHaloCounts.begin(), RecvHaloCounts.end(), 0);
  std::vector<Halo_t> OutHalos(NumberIncomingHaloFragments);

  ExchangeHaloProperties(world, InHalosSorted, OutHalos, SendHaloCounts, SendHaloDisps, RecvHaloCounts, RecvHaloDisps);
  ExchangeHaloParticles (world, TargetRank, InHalosSorted, OutHalos, SendHaloCounts, SendHaloDisps, RecvHaloCounts, RecvHaloDisps);

  InHalos.swap(OutHalos);
}

/* Gathers fragments of FoF groups that have been read by different MPI ranks into
 * the same MPI rank. Once collected, each FoF fragment is merged together to get the
 * complete FoF group. */
void CollectHaloFragments(MpiWorker_t &world, std::vector<Halo_t> &Halos)
{
  ExchangeHaloFragments(world, Halos);
  MergeHaloFragments(Halos);
}
