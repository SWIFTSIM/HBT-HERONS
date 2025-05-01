#include <iostream>
#include <cmath>
#include <algorithm>
#include <random>

#include "config_parser.h"
#include "snapshot.h"
#include "geometric_tree.h"
#include "verify.h"

//
// Make a snapshot subclass with random particle coordinates in a box
//
class RandomSnapshot_t : public Snapshot_t
{

private:

  // A zero vector we can return a reference to
  HBTxyz zero = {0., 0., 0.};

  // Vector of particle coordinates
  std::vector<HBTxyz> all_positions;

public:

  RandomSnapshot_t(HBTInt n, HBTReal boxsize, std::mt19937 &rng)
  {
    // Allocate storage
    all_positions = std::vector<HBTxyz>(n);

    // Assign random coordinates
    std::uniform_real_distribution<HBTReal> dist(0, boxsize);
    for(auto &pos : all_positions) {
      pos[0] = dist(rng);
      pos[1] = dist(rng);
      pos[2] = dist(rng);
    }
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

HBTReal distance_squared(HBTxyz pos1, HBTxyz pos2) {

  if(HBTConfig.PeriodicBoundaryOn) {
    return PeriodicDistance(pos1, pos2);
  } else {
    return (pos1[0]-pos2[0])*(pos1[0]-pos2[0]) +
      (pos1[1]-pos2[1])*(pos1[1]-pos2[1]) +
      (pos1[2]-pos2[2])*(pos1[2]-pos2[2]);
  }
}

int main(int argc, char *argv[])
{
  // Set up repeatable RNG
  std::mt19937 rng;
  rng.seed(0);

  // Make a random snapshot
  const HBTInt N = 1000;
  const HBTReal boxsize = 1.0;
  RandomSnapshot_t snap(N, boxsize, rng);

  // Disable periodicity for now
  HBTConfig.BoxSize = boxsize;
  HBTConfig.BoxHalf = boxsize/2;
  HBTConfig.PeriodicBoundaryOn = false;

  // Build the octree
  GeoTree_t tree;
  tree.Build(snap, N);

  // Make an array of centres for the neighbour search
  const HBTInt Nsearch = 100;
  std::vector<HBTxyz> centre(Nsearch);
  std::uniform_real_distribution<HBTReal> dist(0, boxsize);
  for(auto &pos : centre) {
    pos[0] = dist(rng);
    pos[1] = dist(rng);
    pos[2] = dist(rng);
  }

  // Allocate storage for the results
  std::vector<HBTInt> ngb_idx(Nsearch);

  // Carry out the search
  cout << "Finding neighbours using octree" << std::endl;
  for(HBTInt i=0; i<Nsearch; i+=1) {
    ngb_idx[i] = tree.NearestNeighbour(centre[i], 0.01);
    verify(ngb_idx[i] >= 0);
    verify(ngb_idx[i] < N);
  }

  // Now check the results by brute force:
  // There should be no points closer than the neighbour we found.
  cout << "Checking neighbour distances" << std::endl;
  for(HBTInt i=0; i<Nsearch; i+=1) {

    // Compute distance (squared) to the identified neighbour
    const HBTxyz &cen = centre[i];
    const HBTxyz &ngb_pos = snap.GetComovingPosition(ngb_idx[i]);
    const HBTReal ngb_r2 = distance_squared(cen, ngb_pos);

    // Check distance to all particles
    for(HBTInt j=0; j<N; j+=1) {
      const HBTxyz &part_pos = snap.GetComovingPosition(j);
      if(j != ngb_idx[i]) {
        // Other particles should not be closer than the nearest neighbour
        verify(distance_squared(cen, part_pos) >= ngb_r2);
      }
    }
  }

  cout << "Done." << std::endl;
  return 0;
}
