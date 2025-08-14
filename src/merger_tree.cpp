#include <vector>
#include <cassert>

#include "mpi_wrapper.h"
#include "merger_tree.h"
#include "datatypes.h"
#include "subhalo.h"
#include "locate_ids.h"
#include "argsort.h"

void MergerTreeInfo::Clear()
{
  // Discard the vectors of tracer IDs for the current snapshot
  DescendantTracerIds.clear();
}

void MergerTreeInfo::StoreTracerIds(SubhaloList_t &subhalos, int nr_tracers)
{
  // Store tracers for all resolved subhalos
  for (auto &sub : subhalos)
  {
    if (sub.Nbound > 1)
    {
      DescendantTracerIds[sub.TrackId] = sub.GetMostBoundTracerIds(nr_tracers);
    }
  }
}

void MergerTreeInfo::FindDescendants(SubhaloList_t &Subhalos, MpiWorker_t world)
{

  // Compute descendants for all subhalos in DescendantTracerIds and update
  // DescendantTrackId in the affected Subhalo_t class instances.
  //
  // For subhalos which were resolved in the prevous snapshot,
  // DescendantTrackId stores the TrackId of the subhalo in this snapshot
  // which received the largest number of the tracer particles. Subhalos
  // which were not resolved in the previous snapshot or didn't exist will
  // have DescendantTrackId=-1.
  //
  // This means that if a resolved subhalo in snapshot i becomes unresolved
  // in snapshot i+1 then you can find the TrackId that it merged with by
  // looking at its DescendantTrackId at snapshot i+1. This provides a way to
  // identify a descendant when subhalos become unresolved without sinking.
  //
  // In most cases subhalos which remain resolved will have
  // DescendantTrackId==TrackId and subhalos which sink will have
  // DescendantTrackId==SinkTrackId, although this is not guaranteed.
  if (world.rank() == 0)
  {
    std::cout << "Finding descendants of merged and disrupted subhaloes..." << std::endl;
  }

  std::vector<HBTInt> count_found(0);
  std::vector<HBTInt> trackids_found(0);
  {
    // Count particles in resolved halos
    HBTInt nr_particles = 0;
    for (auto &subhalo : Subhalos)
    {
      if (subhalo.Nbound > 1)
        nr_particles += subhalo.Nbound;
    }

    // Populate vectors of IDs of all particles in resolved subhalos
    std::vector<HBTInt> all_trackids;
    all_trackids.reserve(nr_particles);
    std::vector<HBTInt> all_particle_ids;
    all_particle_ids.reserve(nr_particles);
    for (const auto &subhalo : Subhalos)
    {
      if (subhalo.Nbound > 1)
      {
        for (HBTInt i = 0; i < subhalo.Nbound; i += 1)
        {
          all_trackids.push_back(subhalo.TrackId);
          all_particle_ids.push_back(subhalo.Particles[i].Id);
        }
      }
    }
    assert(all_trackids.size() == nr_particles);
    assert(all_particle_ids.size() == nr_particles);

    // Find the total number of particles we need to search for.
    // This includes tracers from all resolved subhalos
    HBTInt nr_to_find = 0;
    for (const auto &p : DescendantTracerIds)
    {
      // Here p is (TrackId, vector of particle IDs) pair
      nr_to_find += p.second.size();
    }

    // Populate vectors of particles to find
    std::vector<HBTInt> find_particle_ids;
    find_particle_ids.reserve(nr_to_find);
    for (const auto &p : DescendantTracerIds)
    {
      const HBTInt &LostTrackId = p.first;
      const std::vector<HBTInt> &LostTracerIds = p.second;
      for (const auto part_id : LostTracerIds)
      {
        find_particle_ids.push_back(part_id);
      }
    }
    assert(find_particle_ids.size() == nr_to_find);

    // Determine subhalo membership of the particle IDs to find:
    //
    // For each particle ID in find_particle_ids we identify all matching IDs in
    // all_particle_ids (IDs of all bound subhalo particles) and return the TrackId
    // of the subhalo that particle belongs to. Since HBT halos are not fully
    // exclusive some particles may be found multiple times, possibly on other
    // MPI ranks.
    //
    // count_found contains the number of times each entry in find_particle_ids was
    // found.
    //
    // trackids_found contains the trackids of the matching particle IDs and has
    // sum(count_found) entries.
    LocateValuesById(all_particle_ids, all_trackids, MPI_HBT_INT, find_particle_ids, count_found, trackids_found,
                     world.Communicator);
  }

  // Now identify a descendant for each subhalo in the DescendantTracerIds map
  std::vector<std::pair<HBTInt, HBTInt>> DescendantTrackIds;
  {
    HBTInt count_found_offset = 0;
    HBTInt trackids_found_offset = 0;
    DescendantTrackIds.reserve(DescendantTracerIds.size());
    for (const auto &p : DescendantTracerIds)
    {

      // TrackId of the subhalo which we want to find a descendant for
      const HBTInt &LostTrackId = p.first;

      // Vector of tracer particle IDs from the lost subhalo
      const std::vector<HBTInt> &LostTracerIds = p.second;

      // Make a map where the keys are TrackIds of the possible descendants
      // and the values are counts of how many particles went to each one.
      std::map<HBTInt, HBTInt> nr_desc_particles;
      // Loop over tracers to find for this subhalo
      for (HBTInt i = 0; i < LostTracerIds.size(); i += 1)
      {
        // Loop over number of times each tracer was found
        assert(count_found_offset < count_found.size());
        for (HBTInt j = 0; j < count_found[count_found_offset]; j += 1)
        {
          assert(trackids_found_offset < trackids_found.size());
          const HBTInt desc_trackid = trackids_found[trackids_found_offset];
          if (nr_desc_particles.count(desc_trackid) == 0)
          {
            // This is the first particle going to this descendant
            nr_desc_particles[desc_trackid] = 1;
          }
          else
          {
            // Increment the number of particles going to this descendant
            nr_desc_particles[desc_trackid] += 1;
          }
          trackids_found_offset += 1;
        }
        count_found_offset += 1;
      }
      // Now check which subhalo received the largest number of particles.
      // This is the TrackId key with the largest associated value in nr_desc_particles.
      HBTInt max_trackid = SpecialConst::NullTrackId;
      HBTInt max_count = 0;
      for (const auto &p : nr_desc_particles)
      {
        const HBTInt desc_trackid = p.first;
        const HBTInt desc_count = p.second;
        if (desc_count > max_count)
        {
          max_trackid = desc_trackid;
          max_count = desc_count;
        }
      }
      // Store the result in the same order as the DescendantTracerIds map keys
      std::pair<HBTInt, HBTInt> desc_pair(LostTrackId, max_trackid);
      DescendantTrackIds.push_back(desc_pair);
    }
    assert(count_found_offset == count_found.size());
    assert(trackids_found_offset == trackids_found.size());
  }

  // Now we need to update DescendantTrackId in the Subhalo_t structs.
  // Here we assume that subhalos have not migrated between MPI ranks since
  // Unbind() was called, so all entries in the DescendantTracerIds map
  // correspond to subhalos which are stored on this MPI rank.
  {
    // Make a sorting index for the subhalos
    std::vector<HBTInt> subhalo_trackid;
    subhalo_trackid.reserve(Subhalos.size());
    for (auto &sub : Subhalos)
      subhalo_trackid.push_back(sub.TrackId);
    std::vector<HBTInt> subhalo_order = argsort<HBTInt, HBTInt>(subhalo_trackid);

    // Iterate through (TrackId, DescendantTrackId) pairs finding matching subhalos
    HBTInt sub_nr = 0, LocalUpdatedSubhaloes = 0, LocalNoDescendants = 0;
    for (const auto &p : DescendantTrackIds)
    {
      const HBTInt TrackId1 = p.first;  // TrackId of the progenitor
      const HBTInt TrackId2 = p.second; // TrackId of the descendant
      while (Subhalos[subhalo_order[sub_nr]].TrackId < TrackId1)
      {
        sub_nr += 1;
        assert(sub_nr < Subhalos.size());
      }
      assert(Subhalos[subhalo_order[sub_nr]].TrackId == TrackId1);

      /* Only update descendant for subhaloes that are disrupted but which
       * do not have an assigned descendant already. */
      if(!Subhalos[subhalo_order[sub_nr]].IsAlive() && Subhalos[subhalo_order[sub_nr]].DescendantTrackId == SpecialConst::NullTrackId)
      {
        Subhalos[subhalo_order[sub_nr]].DescendantTrackId = TrackId2;
        LocalUpdatedSubhaloes++;
        LocalNoDescendants += (TrackId2 == -1);
      }
    }

    HBTInt TotalUpdatedSubhaloes = 0, TotalNoDescendants = 0;
    MPI_Allreduce(&LocalUpdatedSubhaloes, &TotalUpdatedSubhaloes, 1, MPI_HBT_INT, MPI_SUM, world.Communicator);
    MPI_Allreduce(&LocalNoDescendants   , &TotalNoDescendants   , 1, MPI_HBT_INT, MPI_SUM, world.Communicator);

    if (world.rank() == 0)
    {
      std::cout << "    Number of subhaloes assigned a descendant = " << TotalUpdatedSubhaloes - TotalNoDescendants << std::endl;
      std::cout << "    Number of subhaloes without identified descendants = " << TotalNoDescendants << std::endl;
    }
  }

  // We can now free the vectors of tracer IDs
  Clear();
}
