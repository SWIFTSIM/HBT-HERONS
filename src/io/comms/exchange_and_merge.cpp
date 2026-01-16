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
#include "../../halo_particle_iterator.h"
#include "../../particle_exchanger.h"

struct HaloInfo_t
{
  HBTInt id;
  HBTReal m;
  HBTxyz x;
  int order;
};

static void create_MPI_HaloInfo_t(MPI_Datatype &dtype)
{
  HaloInfo_t p;
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
  RegisterAttr(id, MPI_HBT_INT, 1);
  RegisterAttr(m, MPI_HBT_REAL, 1);
  RegisterAttr(x[0], MPI_HBT_REAL, 3);
  RegisterAttr(order, MPI_INT, 1);
#undef RegisterAttr
  assert(i <= NumAttr);

  MPI_Type_create_struct(i, blockcounts, offsets, oldtypes, &dtype);
  MPI_Type_create_resized(dtype, (MPI_Aint)0, extent, &dtype);
  MPI_Type_commit(&dtype);
#undef NumAttr
}

inline bool CompHaloInfo_Id(const HaloInfo_t &a, const HaloInfo_t &b)
{
  return a.id < b.id;
}

inline bool CompHaloInfo_Order(const HaloInfo_t &a, const HaloInfo_t &b)
{
  return a.order < b.order;
}

inline bool CompHaloId(const Halo_t &a, const Halo_t &b)
{
  return a.HaloId < b.HaloId;
}

static double ReduceHaloPosition(vector<HaloInfo_t>::iterator it_begin, vector<HaloInfo_t>::iterator it_end, HBTxyz &x)
{
  HBTInt i, j;
  double sx[3], origin[3], msum;

  if (it_begin == it_end)
    return 0.;
  if (it_begin + 1 == it_end)
  {
    copyHBTxyz(x, it_begin->x);
    return it_begin->m;
  }

  sx[0] = sx[1] = sx[2] = 0.;
  msum = 0.;
  if (HBTConfig.PeriodicBoundaryOn)
    for (j = 0; j < 3; j++)
      origin[j] = it_begin->x[j];

  for (auto it = it_begin; it != it_end; ++it)
  {
    HBTReal m = it->m;
    msum += m;
    for (j = 0; j < 3; j++)
      if (HBTConfig.PeriodicBoundaryOn)
        sx[j] += NEAREST(it->x[j] - origin[j]) * m;
      else
        sx[j] += it->x[j] * m;
  }

  for (j = 0; j < 3; j++)
  {
    sx[j] /= msum;
    if (HBTConfig.PeriodicBoundaryOn)
    {
      sx[j] += origin[j];
      x[j] = position_modulus(sx[j], HBTConfig.BoxSize);
    }
    else
      x[j] = sx[j];
  }
  return msum;
}

static void ReduceHaloRank(vector<HaloInfo_t>::iterator it_begin, vector<HaloInfo_t>::iterator it_end, HBTxyz &step,
                           vector<int> &dims)
{
  HBTxyz x;
  ReduceHaloPosition(it_begin, it_end, x);
  int rank = AssignCell(x, step, dims);
  for (auto it = it_begin; it != it_end; ++it)
    it->id = rank; // store destination rank in id.
}

static vector<IdRank_t> DecideTargetProcessor(MpiWorker_t &world, vector<Halo_t> &Halos)
{
  int this_rank = world.rank();
  for (auto &&h : Halos)
    h.Mass = AveragePosition(h.ComovingAveragePosition, h.Particles.data(), h.Particles.size());

  vector<HaloInfo_t> HaloInfoSend(Halos.size()), HaloInfoRecv;
  for (HBTInt i = 0; i < Halos.size(); i++)
  {
    HaloInfoSend[i].id = Halos[i].HaloId;
    HaloInfoSend[i].m = Halos[i].Mass;
    HaloInfoSend[i].x = Halos[i].ComovingAveragePosition;
  }
  HBTInt MaxHaloId = 0;
  if (Halos.size())
    MaxHaloId = Halos.back().HaloId;
  MPI_Allreduce(MPI_IN_PLACE, &MaxHaloId, 1, MPI_HBT_INT, MPI_MAX, world.Communicator);
  HBTInt ndiv = (++MaxHaloId) / world.size();
  if (MaxHaloId % world.size())
    ndiv++;
  vector<int> SendSizes(world.size(), 0), SendOffsets(world.size()), RecvSizes(world.size()), RecvOffsets(world.size());
  for (HBTInt i = 0; i < Halos.size(); i++)
  {
    int idiv = Halos[i].HaloId / ndiv;
    SendSizes[idiv]++;
  }
  CompileOffsets(SendSizes, SendOffsets);
  MPI_Alltoall(SendSizes.data(), 1, MPI_INT, RecvSizes.data(), 1, MPI_INT, world.Communicator);
  int nhalo_recv = CompileOffsets(RecvSizes, RecvOffsets);
  HaloInfoRecv.resize(nhalo_recv);
  MPI_Datatype MPI_HaloInfo_t;
  create_MPI_HaloInfo_t(MPI_HaloInfo_t);
  MPI_Alltoallv(HaloInfoSend.data(), SendSizes.data(), SendOffsets.data(), MPI_HaloInfo_t, HaloInfoRecv.data(),
                RecvSizes.data(), RecvOffsets.data(), MPI_HaloInfo_t, world.Communicator);
  for (int i = 0; i < nhalo_recv; i++)
    HaloInfoRecv[i].order = i;
  sort(HaloInfoRecv.begin(), HaloInfoRecv.end(), CompHaloInfo_Id);
  list<int> haloid_offsets;
  HBTInt curr_id = -1;
  for (int i = 0; i < nhalo_recv; i++)
  {
    if (curr_id != HaloInfoRecv[i].id)
    {
      haloid_offsets.push_back(i);
      curr_id = HaloInfoRecv[i].id;
    }
  }
  haloid_offsets.push_back(nhalo_recv);
  // combine coordinates and determine target
  auto dims = ClosestFactors(world.size(), 3);
  HBTxyz step;
  for (int i = 0; i < 3; i++)
    step[i] = HBTConfig.BoxSize / dims[i];
  auto it_end = haloid_offsets.end();
  --it_end;
  for (auto it = haloid_offsets.begin(); it != it_end; it++)
  {
    auto it_next = it;
    ++it_next;
    ReduceHaloRank(HaloInfoRecv.begin() + *it, HaloInfoRecv.begin() + *it_next, step, dims);
  }
  sort(HaloInfoRecv.begin(), HaloInfoRecv.end(), CompHaloInfo_Order);
  // send back
  MPI_Alltoallv(HaloInfoRecv.data(), RecvSizes.data(), RecvOffsets.data(), MPI_HaloInfo_t, HaloInfoSend.data(),
                SendSizes.data(), SendOffsets.data(), MPI_HaloInfo_t, world.Communicator);
  MPI_Type_free(&MPI_HaloInfo_t);

  vector<IdRank_t> TargetRank(Halos.size());
  for (HBTInt i = 0; i < TargetRank.size(); i++)
  {
    TargetRank[i].Id = i;
    TargetRank[i].Rank = HaloInfoSend[i].id;
  }
  std::sort(TargetRank.begin(), TargetRank.end(), CompareRank);
  return TargetRank;
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
                            //  std::vector<Halo_t>::iterator NewHalos,
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
