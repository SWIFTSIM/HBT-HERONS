#include <mpi.h>

#include "verify.h"
#include "mpi_wrapper.h"
#include "subhalo.h"
#include "merger_tree.h"

//
// Basic test of merger tree building code
//
// We create subhalos with TrackId and bound particle IDs set, and store their
// tracer particle IDs as if they were resolved in the previous output. We then
// change the subhalos to represent the current output, where:
//
//   i % 3 == 0 : the subhalo is still resolved, and accretes the particles of
//                subhalo i+1
//   i % 3 == 1 : the subhalo disrupted, and its particles are now bound to
//                subhalo i-1, which should be its descendant
//   i % 3 == 2 : the subhalo disrupted, and its particles are not bound to any
//                subhalo, so it should have no descendant
//
// Only disrupted subhalos should be assigned a descendant, so resolved
// subhalos keep DescendantTrackId = -1.
//
int main(int argc, char *argv[])
{

  MPI_Init(&argc, &argv);
  int comm_size;
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
  int comm_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
  MpiWorker_t world(MPI_COMM_WORLD);

  MergerTreeInfo merger_tree;

  const int subhalos_per_rank = 100;

  // Make a vector of subhalos
  SubhaloList_t Subhalo(subhalos_per_rank);

  // Assign sequential TrackIds across MPI ranks
  for (int i = 0; i < subhalos_per_rank; i += 1)
  {
    Subhalo[i].TrackId = (comm_rank * subhalos_per_rank) + i;
  }

  // Set bound particle IDs in the subhalos, as in the previous output
  const int particles_per_sub = 100;
  for (int i = 0; i < subhalos_per_rank; i += 1)
  {
    Subhalo[i].Nbound = particles_per_sub;
    Subhalo[i].Particles.resize(particles_per_sub);
    for (int j = 0; j < particles_per_sub; j += 1)
    {
      Subhalo[i].Particles[j].Id = (Subhalo[i].TrackId * particles_per_sub) + j;
#ifndef DM_ONLY
      Subhalo[i].Particles[j].Type = TypeDM;
#endif
    }
  }

  // Identify tracer particle IDs for each subhalo (just the first nr_tracers particles)
  const int nr_tracers = 5;
  merger_tree.StoreTracerIds(Subhalo, nr_tracers);

  // Change the subhalos to represent the current output
  const int disrupted_snapshot = 1;
  for (int i = 0; i < subhalos_per_rank; i += 1)
  {
    if (i % 3 == 0)
      continue;

    // Disrupted subhalos have no bound particles
    if (i % 3 == 1)
    {
      Subhalo[i - 1].Particles.insert(Subhalo[i - 1].Particles.end(), Subhalo[i].Particles.begin(),
                                      Subhalo[i].Particles.end());
      Subhalo[i - 1].Nbound = Subhalo[i - 1].Particles.size();
    }
    Subhalo[i].Particles.clear();
    Subhalo[i].Nbound = 0;
    Subhalo[i].SnapshotOfDeath = disrupted_snapshot;
  }

  // Use tracers to identify subhalo descendants
  merger_tree.FindDescendants(Subhalo, world);

  // Check results
  for (int i = 0; i < subhalos_per_rank; i += 1)
  {
    if (i % 3 == 0)
      verify(Subhalo[i].DescendantTrackId == SpecialConst::NullTrackId);
    else if (i % 3 == 1)
      verify(Subhalo[i].DescendantTrackId == Subhalo[i - 1].TrackId);
    else
      verify(Subhalo[i].DescendantTrackId == SpecialConst::NullTrackId);
  }

  MPI_Finalize();

  return 0;
}
