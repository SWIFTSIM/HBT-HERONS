using namespace std;
#include <iostream>
// #include <iomanip>
#include <assert.h>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <numeric>
#include <sstream>
#include <string>
#include <typeinfo>

#include "../mpi_wrapper.h"
#include "../mymath.h"
#include "../snapshot.h"
#include "apostle_io.h"
#include "gadget_io.h"
#include "swiftsim_io.h"

void ParticleSnapshot_t::Load(MpiWorker_t &world, int snapshot_index)
{
  Clear();
  SetSnapshotIndex(snapshot_index);

  if (HBTConfig.SnapshotFormat == "gadget")
  {
    GadgetReader_t(world, SnapshotId, Particles, Cosmology);
  }
  else if (HBTConfig.SnapshotFormat == "apostle")
  {
    ApostleReader_t().LoadSnapshot(world, SnapshotId, Particles, Cosmology);
  }
  else if (HBTConfig.SnapshotFormat == "swiftsim")
  {
    SwiftSimReader_t().LoadSnapshot(world, SnapshotId, Particles, Cosmology);

    // Load particle splitting information if we need it.
#ifndef DM_ONLY
    if (HBTConfig.ParticlesSplit)
    {
      SwiftSimReader_t().ReadParticleSplits(ParticleSplitMap, SnapshotId);
    }
#endif
  }
  else if (HBTConfig.SnapshotFormat == "mysnapshot")
  {
    /* Insert your snapshot reader here, and include relevant information in the header
     * if necessary. Essential things are to set cosmology and fill up Particle vector.
     * LoadMySnapshot(SnapshotId, Particles, Cosmology); */
  }
  else
    throw(runtime_error("unknown SnapshotFormat " + HBTConfig.SnapshotFormat));

  global_timer.Tick("snap_io", world.Communicator);

  ExchangeParticles(world);
  global_timer.Tick("snap_exchange", world.Communicator);

  FillParticleHash();
  global_timer.Tick("snap_hash", world.Communicator);

  if (world.rank() == 0)
    std::cout << NumberOfParticlesOnAllNodes << " particles loaded from snapshot " << SnapshotId << " (SnapshotIndex = " \
              << snapshot_index << ")" << std::endl;
}

#ifdef TEST_snapshot_io
#include "../config_parser.h"

int main(int argc, char **argv)
{
  HBTConfig.ParseConfigFile(argv[1]);
  ParticleSnapshot_t snapshot;
  snapshot.Load(HBTConfig.MaxSnapshotIndex, true);
  cout << snapshot.GetNumberOfParticles() << endl;
  cout << snapshot.GetParticleId(10) << endl;
  cout << snapshot.GetComovingPosition(10) << endl;
  cout << snapshot.GetParticleMass(10) << ',' << snapshot.GetParticleMass(100) << endl;
  cout << snapshot.GetParticleIndex(snapshot.GetParticleId(10)) << endl;
  return 0;
}
#endif
