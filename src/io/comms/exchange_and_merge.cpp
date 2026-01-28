using namespace std;
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
#include "../../halo_particle_iterator.h"
#include "../../particle_exchanger.h"

struct HaloFragment_t
{
  HBTInt HaloId;
  HBTInt NumberParticles;
  int OriginalRank;
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

/* Assigns MPI tasks to haloes so that the number of FoF particles is approximately
 * equal across ranks. */
static std::vector<IdRank_t> DecideTargetProcessor(MpiWorker_t &world, std::vector<Halo_t> &Halos)
{
  /* First we compute total FoF group sizes. Note that the ordering and vector
   * size of HaloSizes and Halos is not the same! */
  std::vector<HaloFragment_t> HaloSizes = GetTotalHaloSizes(world, Halos);
  std::vector<IdRank_t> TargetTask(Halos.size());
  return TargetTask;
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