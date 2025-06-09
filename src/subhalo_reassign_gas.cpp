#include <algorithm>
#include <iostream>
#include <numeric>
#include <map>
#include <omp.h>

#include "datatypes.h"
#include "subhalo.h"
#include "config_parser.h"
#include "geometric_tree.h"

/*
 * Class providing the minimum snapshot functionality for the octree
 * neighbour search to work, given a vector of positions.
 */
class TracerSnapshot_t : public Snapshot_t
{

private:

  // A zero vector we can return a reference to
  HBTxyz zero = {0., 0., 0.};

  // Vector of particle coordinates
  std::vector<HBTxyz> all_positions;

public:

  TracerSnapshot_t(std::vector<HBTxyz> pos)
  {
    all_positions = std::move(pos); // Invalidates input vector
  }

  HBTInt size() const
  {
    return all_positions.size();
  }

  HBTInt GetId(const HBTInt index) const
  {
    return index;
  }

  const HBTxyz &GetComovingPosition(const HBTInt index) const
  {
    return all_positions[index];
  }

  const HBTxyz GetPhysicalVelocity(const HBTInt index) const
  {
    return zero;
  }

  HBTReal GetMass(const HBTInt index) const
  {
    return 1.0;
  }

  HBTReal GetInternalEnergy(HBTInt index) const
  {
    return 0.0;
  }
};

void SubhaloSnapshot_t::ReassignParticles(MpiWorker_t &world, HaloSnapshot_t &halo_snap)
{
  HBTInt nr_reassigned = 0;

  // Check if particle reassignment is enabled
  if(!HBTConfig.ReassignParticles) {
    if (world.rank() == 0)
      cout << "  Not reassigning particles\n";
    return;
  }

  // Loop over FoF groups
  HBTInt NumHalos = MemberTable.SubGroups.size();
#pragma omp parallel for schedule(dynamic, 1) reduction(+:nr_reassigned)
  for (HBTInt haloid = 0; haloid < NumHalos; haloid++)
    {
      // Get indexes of subhalos in this FoF group
      auto &subgroup = MemberTable.SubGroups[haloid];

      // Count tracer and non-tracer type particles
      HBTInt nr_tracers = 0;
      HBTInt nr_non_tracers = 0;
      for(auto subid : subgroup)
        {
          for(auto &part : Subhalos[subid].Particles)
            {
              assert(part.Id != SpecialConst::NullParticleId); // Should not have null particles at this point
              if(part.IsTracer())
                {
                  nr_tracers += 1;
                }
              else
                {
                  nr_non_tracers += 1;
                }
            }
        }

      // If there are no non-tracers, we can skip this halo
      if(nr_non_tracers > 0)
        {
          // Store the number of particles in each subhalo.
          // This is so that moved particles can be appended to the particle
          // arrays but we can still iterate over the original particles only.
          std::map<HBTInt,HBTInt> sublen;
          for(auto subid : subgroup)
            sublen[subid] = Subhalos[subid].Particles.size();

          // Store position and subhalo index of the tracers
          std::vector<HBTxyz> tracer_pos(nr_tracers);
          std::vector<HBTInt> tracer_subid(nr_tracers);
          nr_tracers = 0;
          for(auto subid : subgroup)
            {
              for(HBTInt i=0; i<sublen[subid]; i+=1)
                {
                  auto &part = Subhalos[subid].Particles[i];
                  if(part.IsTracer())
                    {
                      tracer_pos[nr_tracers] = part.ComovingPosition;
                      tracer_subid[nr_tracers] = subid;
                      nr_tracers += 1;
                    }
                }
            }

          // Build a tree for the neighbour search
          TracerSnapshot_t tracer_snap(tracer_pos);
          GeoTree_t tree;
          tree.Build(tracer_snap);

          // For each non-tracer particle (bound or not), identify the nearest neighbour tracer type particles
          for(auto subid : subgroup)
            {
              for(HBTInt i=0; i<sublen[subid]; i+=1)
                {
                  auto &part = Subhalos[subid].Particles[i];
                  assert(part.Id != SpecialConst::NullParticleId); // Should not have null particles at this point
                  if(!part.IsTracer())
                    {
                      const HBTInt nr_ngbs = HBTConfig.NumNeighboursForReassignment;
                      const HBTxyz centre = part.ComovingPosition;
                      // TODO: better initial search radius guess
                      //       could use size of the tree node containing the point?
                      //       not sure why we have to provide a guess at all!
                      std::vector<HBTInt> neighbour_list = tree.NearestNeighbours(centre, nr_ngbs, HBTConfig.SofteningHalo*0.1);
                      // Check if any of the neighbours are in the current subhalo
                      bool found = false;
                      for(auto ngb_idx : neighbour_list)
                        {
                          if(tracer_subid[ngb_idx] == subid)
                            {
                              found=true;
                              break;
                            }
                        }
                      if(!found)
                        {
                          // Current subhalo is not on the neighbour list, so we
                          // may want to move this particle. Find nearest neighbour.
                          const HBTInt ngb_idx = neighbour_list.back();
                          const HBTInt ngb_subid = tracer_subid[ngb_idx];
                          if((ngb_subid >= 0) && (ngb_subid != subid))
                            {
                              // Append the particle to the other subhalo's source list
                              Subhalos[ngb_subid].Particles.push_back(part);
                              // Flag the particle for removal from this subhalo
                              part.Id = SpecialConst::NullParticleId;
                              nr_reassigned += 1;
                            }
                        }
                    }
                }
            }

          // Now tidy up any particles we flagged for removal
          for(auto subid : subgroup) {
            // Ensure Nbound is not greater than the total number of particles
            if(Subhalos[subid].Nbound > Subhalos[subid].Particles.size())
              Subhalos[subid].Nbound = Subhalos[subid].Particles.size();
            // Remove null particles and update tracer index
            Subhalos[subid].KickNullParticles();
          }
        }
      // Next FoF halo
    }

  // Compute total number of particles which were reassigned
  HBTInt nr_reassigned_tot;
  MPI_Allreduce(&nr_reassigned, &nr_reassigned_tot, 1, MPI_HBT_INT, MPI_SUM, world.Communicator);
  if (world.rank() == 0)
    cout << "  Total number of particles reassigned = " << nr_reassigned_tot << "\n";
}
