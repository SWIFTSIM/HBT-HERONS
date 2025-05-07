#include <vector>
#include <set>

#include "datatypes.h"
#include "subhalo.h"
#include "geometric_tree.h"

/*
 * Class providing the minimum snapshot functionality for the octree
 * neighbour search to work, given a vector of positions.
 */
class SubhaloParticleSnapshot_t : public Snapshot_t
{

private:

  // A zero vector we can return a reference to
  HBTxyz zero = {0., 0., 0.};

  // Vector of particle coordinates
  std::vector<HBTxyz> all_positions;

public:

  SubhaloParticleSnapshot_t(std::vector<HBTxyz> pos)
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
  Assign newly accreted FoF particles to resolved subhalos

  For each FoF particle not already in a subhalo's source list, we
  find the nearest neighbour bound subhalo particle and add the particle
  to that subhalo.

  Only called if we have at least one pre-existing subhalo in the halo.
  Note that this is called inside a multithreaded loop over FoF halos.
 */
void SubhaloSnapshot_t::FeedSubhalos(Halo_t &Host, MemberShipTable_t::MemberList_t &Members)
{
  // Count resolved subhalos in this FoF
  HBTInt nr_resolved = 0;
  for(auto subid : Members)
    if(Subhalos[subid].IsAlive() && (Subhalos[subid].Nbound > 0))nr_resolved += 1;

  // If we have no resolved subhalos, fall back to HBTs original behaviour
  // of assigning all FoF particles to the central.
  if(nr_resolved == 0)
    {
      auto &central = Subhalos[Members[0]];
      central.Particles.swap(Host.Particles);
      central.Nbound = central.Particles.size();
      return;
    }

  // Otherwise, build an octree containing all bound particles in this FoF
  // group's subhalos. First count how many bound particles we have.
  HBTInt nr_bound = 0;
  for(auto subid : Members)
    if(Subhalos[subid].IsAlive())nr_bound += Subhalos[subid].Nbound;

  // Then store the particle positions and subids
  std::vector<HBTxyz> bound_pos(nr_bound);
  std::vector<HBTInt> bound_subid(nr_bound);
  nr_bound = 0;
  for(auto subid : Members)
    {
      if(Subhalos[subid].IsAlive())
        {
          // Loop over bound particles only
          for(HBTInt i=0; i<Subhalos[subid].Nbound; i+=1)
            {
              bound_pos[nr_bound] = Subhalos[subid].Particles[i].ComovingPosition;
              bound_subid[nr_bound] = subid;
              nr_bound += 1;
            }
        }
    }

  // And build the tree
  SubhaloParticleSnapshot_t bound_snap(bound_pos);
  GeoTree_t tree;
  tree.Build(bound_snap);

  // Now we need to determine which FoF particles are new accretion, i.e. are
  // not on the source list of any subhalo in this halo. Make a set containing
  // the IDs of all of the bound particles.
  std::set<HBTInt> all_bound_ids;
  for(auto subid : Members)
    {
      if(Subhalos[subid].IsAlive())
        {
          // Loop over bound and unbound source particles
          for(auto p : Subhalos[subid].Particles)
            all_bound_ids.insert(p.Id);
        }
    }

  // Now determine which subhalo gets each newly accreted FoF particle
  for(auto p : Host.Particles) {
    if(all_bound_ids.count(p.Id) == 0) {
      // This is a newly accreted particle, so find its nearest bound subhalo
      // neighbour particle.
      HBTInt neighbour_index = tree.NearestNeighbour(p.ComovingPosition, 1.0e-4);
      assert(neighbour_index >= 0);
      assert(neighbour_index < nr_bound);
      // Add this particle to the identified subhalo
      HBTInt subid = bound_subid[neighbour_index];
      Subhalos[subid].Particles.push_back(p);
    }
  }

  return;
}
