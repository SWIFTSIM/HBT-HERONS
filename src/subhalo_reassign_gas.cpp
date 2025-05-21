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

void SubhaloSnapshot_t::ReassignParticles()
{
  switch(HBTConfig.ReassignParticles)
    {
    case 0:
      // We're not reassigning particles
      return;
    case 1:
      // Non-tracer-type particles are moved to tracer type neighbours
      ReassignNonTracerParticles();
      return;
    case 2:
      // Particles of any type are moved if all neighbours are in another halo
      ReassignParticlesOfAnyType();
      return;
    default:
      // Invalid parameter value
      std::cerr << "Invalid value for ReassignParticles" << std::endl;
      std::abort();
    }
}

/*
  For each non-tracer type particle, identify the nearest tracer type particle
  in any subhalo in the same FoF group and assign it to the subhalo the tracer
  is bound to, if any.

  All non-tracer type particles, bound or not, are considered for reassignment.
*/
void SubhaloSnapshot_t::ReassignNonTracerParticles()
{
  // Loop over FoF groups
  HBTInt NumHalos = MemberTable.SubGroups.size();
#pragma omp parallel for schedule(dynamic, 1)
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
                  if((part.Id != SpecialConst::NullParticleId) && (!part.IsTracer()))
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
  // Done.
}


/*
  For each particle, identify N nearest neighbours. If all neighbours
  are in the same halo (but not the one the particle is in) then move
  the particle to that halo.
*/
void SubhaloSnapshot_t::ReassignParticlesOfAnyType()
{
  // Loop over FoF groups
  HBTInt NumHalos = MemberTable.SubGroups.size();
#pragma omp parallel for schedule(dynamic, 1)
  for (HBTInt haloid = 0; haloid < NumHalos; haloid++)
    {
      // Get indexes of subhalos in this FoF group
      auto &subgroup = MemberTable.SubGroups[haloid];

      // Store the number of particles in each subhalo.
      // This is so that moved particles can be appended to the particle
      // arrays but we can still iterate over the original particles only.
      HBTInt nr_tot = 0;
      std::map<HBTInt,HBTInt> sublen;
      for(auto subid : subgroup)
        {
          sublen[subid] = Subhalos[subid].Particles.size();
          nr_tot += Subhalos[subid].Particles.size();
        }

      // Store position and subhalo index of all subhalo particles
      std::vector<HBTxyz> part_pos(nr_tot);
      std::vector<HBTInt> part_subid(nr_tot);
      nr_tot = 0;
      for(auto subid : subgroup)
        {
          for(HBTInt i=0; i<sublen[subid]; i+=1)
            {
              auto &part = Subhalos[subid].Particles[i];
              part_pos[nr_tot] = part.ComovingPosition;
              part_subid[nr_tot] = subid;
              nr_tot += 1;
            }
        }

      // If there are no particles (is that possible?), skip this FoF group
      if(nr_tot > 0)
        {
          // Build a tree for the neighbour search
          TracerSnapshot_t part_snap(part_pos);
          GeoTree_t tree;
          tree.Build(part_snap);

          // Loop over subhalos in the FoF
          for(auto subid : subgroup)
            {
              // Loop over particles in the subhalo
              for(HBTInt i=0; i<sublen[subid]; i+=1)
                {
                  auto &part = Subhalos[subid].Particles[i];
                  if(part.Id != SpecialConst::NullParticleId)
                    {
                      // Find neighbours. Note that neighbour_list is in descending
                      // order of distance and the last neighbour will usually be
                      // the particle we're considering reassigning (but might not be
                      // if several particles have identical coordinates).
                      const HBTInt nr_ngbs = HBTConfig.NumNeighboursForReassignment;
                      const HBTxyz centre = part.ComovingPosition;
                      std::vector<HBTInt> neighbour_list = tree.NearestNeighbours(centre, nr_ngbs, HBTConfig.SofteningHalo*0.1);
                      // Now check the what subhalos the neighbours are in. Here
                      // we're going to store the subid of the first neighbour which
                      // is in a different subhalo and count how many neighbours are in
                      // that subhalo.
                      HBTInt new_subid = -1;
                      HBTInt new_subid_count = 0;
                      for(auto ngb_idx : neighbour_list)
                        {
                          HBTInt ngb_subid = part_subid[ngb_idx];
                          if(ngb_subid != subid)
                            {
                              // This particle is in a different subhalo
                              if(new_subid == -1)new_subid = ngb_subid;
                              assert(ngb_subid >= 0); // We only search for neighbours in subhalos
                              assert(new_subid >= 0); // Should have been set now
                              if(ngb_subid == new_subid)new_subid_count += 1;
                            }
                        }
                      // The "neighbours" include the particle itself so if nr_ngbs-1
                      // neighbours are in another subhalo, we move the particle.
                      if((new_subid >= 0) && (new_subid_count==nr_ngbs-1))
                        {
                          // Append the particle to the other subhalo's source list
                          Subhalos[new_subid].Particles.push_back(part);
                          // Flag the particle for removal from this subhalo
                          part.Id = SpecialConst::NullParticleId;
                        }
                    }
                }
            }

          // Now tidy up any particles we flagged for removal
          for(auto subid : subgroup)
            {
              // Ensure Nbound is not greater than the total number of particles
              if(Subhalos[subid].Nbound > Subhalos[subid].Particles.size())
                Subhalos[subid].Nbound = Subhalos[subid].Particles.size();
              // Remove null particles and update tracer index
              Subhalos[subid].KickNullParticles();
            }
        }
      // Next FoF halo
    }
  // Done.
}
