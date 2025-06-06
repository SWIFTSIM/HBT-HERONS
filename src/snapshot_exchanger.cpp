using namespace std;
#include <iostream>
// #include <iomanip>
#include <assert.h>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <sstream>
#include <string>
#include <typeinfo>

#include "mpi_wrapper.h"
#include "mymath.h"
#include "pairwise_alltoallv.h"
#include "snapshot.h"
#include "sort_by_hash.h"

inline int GetGrid(HBTReal x, HBTReal step, int dim)
{
  int i = floor(x / step);
  if (i < 0)
    i = 0;
  if (i >= dim)
    i = dim - 1;
  return i;
}
inline int AssignCell(HBTxyz &Pos, const HBTReal step[3], const vector<int> &dims)
{
#define GRIDtoRank(g0, g1, g2) (((g0) * dims[1] + (g1)) * dims[2] + (g2))
#define GID(i) GetGrid(Pos[i], step[i], dims[i])
  return GRIDtoRank(GID(0), GID(1), GID(2));
}

vector<HBTInt> ParticleSnapshot_t::PartitionParticles(MpiWorker_t &world)
{
  return sort_by_hash(Particles, world.size());
}

inline bool CompParticleId(const Particle_t &a, const Particle_t &b)
{
  return a.Id < b.Id;
}

/*
  This sends each particle to an MPI rank based on the hash of its ID.
  This is so that when we need to find a particle by ID we can easily
  compute which rank it is stored on by hashing the ID.

  Within a rank particles are sorted by ID.
*/
void ParticleSnapshot_t::ExchangeParticles(MpiWorker_t &world)
{

  {
    HBTInt np = Particles.size();
    MPI_Allreduce(&np, &NumberOfParticlesOnAllNodes, 1, MPI_HBT_INT, MPI_SUM, world.Communicator);
  }

  /* We should not communicate particles for GADGET-4 just yet, as the groups
   * are defined based on particle offsets. Communicating will therefore 
   * mix up the order of particles and give bogus results. */
  if(HBTConfig.GroupFileFormat == "gadget4_hdf")
    return;

  vector<HBTInt> SendOffsets = PartitionParticles(world);
  SendOffsets.push_back(Particles.size());
  vector<HBTInt> SendSizes(world.size(), 0);
  for (int i = 0; i < world.size(); i++)
    SendSizes[i] = SendOffsets[i + 1] - SendOffsets[i];

  vector<HBTInt> ReceiveSizes(world.size(), 0), ReceiveOffsets(world.size());
  MPI_Alltoall(SendSizes.data(), 1, MPI_HBT_INT, ReceiveSizes.data(), 1, MPI_HBT_INT, world.Communicator);
  vector<Particle_t> ReceivedParticles;
  ReceivedParticles.resize(CompileOffsets(ReceiveSizes, ReceiveOffsets));

#ifndef NDEBUG
  // Bounds check offsets and counts in debug mode
  for (int rank = 0; rank < world.size(); rank += 1)
  {
    if (SendSizes[rank] > 0)
    {
      assert(SendOffsets[rank] >= 0);
      assert(SendOffsets[rank] + SendSizes[rank] <= Particles.size());
    }
    if (ReceiveSizes[rank] > 0)
    {
      assert(ReceiveOffsets[rank] >= 0);
      assert(ReceiveOffsets[rank] + ReceiveSizes[rank] <= ReceivedParticles.size());
    }
  }
#endif

  MPI_Datatype MPI_HBT_Particle;
  Particle_t().create_MPI_type(MPI_HBT_Particle);

  Pairwise_Alltoallv(Particles, SendSizes, SendOffsets, MPI_HBT_Particle, ReceivedParticles, ReceiveSizes,
                     ReceiveOffsets, MPI_HBT_Particle, world.Communicator);

  MPI_Barrier(world.Communicator);

  MPI_Type_free(&MPI_HBT_Particle);

  Particles.swap(ReceivedParticles);

  sort(Particles.begin(), Particles.end(), CompParticleId);
}
