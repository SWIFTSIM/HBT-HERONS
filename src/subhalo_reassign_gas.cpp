#include <algorithm>
#include <iostream>
#include <numeric>
#include <map>
#include <omp.h>

#include "datatypes.h"
#include "subhalo.h"
#include "config_parser.h"
#include "geometric_tree.h"

#ifndef DM_ONLY

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


/*
  For each gas particle, identify the nearest tracer type particle in
  any subhalo in the same FoF group and assign the gas particle to the
  subhalo the tracer is bound to, if any.

  All gas particles, bound or not, are considered for reassignment.
*/
void SubhaloSnapshot_t::ReassignGasParticles()
{

  // Loop over FoF groups
  HBTInt NumHalos = MemberTable.SubGroups.size();
  for (HBTInt haloid = 0; haloid < NumHalos; haloid++)
    {
      // Get indexes of subhalos in this FoF group
      auto &subgroup = MemberTable.SubGroups[haloid];

      // Count gas and tracer type particles
      HBTInt nr_tracers = 0;
      HBTInt nr_gas = 0;
      for(auto subid : subgroup)
        {
          for(auto &part : Subhalos[subid].Particles)
            {
              if(part.IsTracer())nr_tracers += 1;
              if(part.Type==TypeGas)nr_gas += 1;
            }
        }

      // If there's no gas, we can skip this halo
      if(nr_gas > 0)
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
              const HBTInt nbound = Subhalos[subid].Nbound;
              for(HBTInt i=0; i<sublen[subid]; i+=1)
                {
                  auto &part = Subhalos[subid].Particles[i];
                  if(part.IsTracer())
                    {
                      tracer_pos[nr_tracers] = part.ComovingPosition;
                      if(i < nbound)
                        tracer_subid[nr_tracers] = subid;
                      else
                        tracer_subid[nr_tracers] = -1; // Tracer is not bound to any subhalo
                      nr_tracers += 1;
                    }
                }
            }

          // Build a tree for the neighbour search
          TracerSnapshot_t tracer_snap(tracer_pos);
          GeoTree_t tree;
          tree.Build(tracer_snap);

          // For each gas particle (bound or not), identify the nearest neighbour tracer type particle
          for(auto subid : subgroup)
            {
              for(HBTInt i=0; i<sublen[subid]; i+=1)
                {
                  auto &part = Subhalos[subid].Particles[i];
                  const HBTxyz centre = part.ComovingPosition;
                  const HBTInt ngb_idx = tree.NearestNeighbour(centre, 1.0e-4);
                  const HBTInt ngb_subid = tracer_subid[ngb_idx];

                  // Check if we need to move this gas particle
                  if((ngb_subid >= 0) && (ngb_subid != subid))
                    {
                      // Append the particle to the other subhalo's source list
                      Subhalos[ngb_subid].Particles.push_back(part);
                      // Flag the particle for removal from this subhalo
                      part.Id = SpecialConst::NullParticleId;
                    }
                }
            }

          // Now tidy up any particles we flagged for removal and update Nbound
          for(auto subid : subgroup)
            Subhalos[subid].KickNullParticles();
        }
      // Next FoF halo
    }
  // Done.
}
#endif
