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

void test_neighbour_search(HBTInt N, HBTInt Nsearch, HBTReal boxsize, bool periodic, std::mt19937 &rng)
{
  const HBTInt nr_ngb = 10;

  std::cout << "Start test with N=" << N << ", Nsearch = " << Nsearch << ", boxsize = " << boxsize;
  if(periodic)std::cout << " (periodic)";
  std::cout << endl;

  // Make a random snapshot
  RandomSnapshot_t snap(N, boxsize, rng);

  // Box size and periodicity used by the tree are set in the HBTConfig object
  HBTConfig.BoxSize = boxsize;
  HBTConfig.BoxHalf = boxsize/2;
  HBTConfig.PeriodicBoundaryOn = false;

  // Build the octree
  GeoTree_t tree;
  tree.Build(snap, N);

  // Make an array of centres for the neighbour search
  std::vector<HBTxyz> centre(Nsearch);
  std::uniform_real_distribution<HBTReal> dist(0, boxsize);
  for(auto &pos : centre) {
    pos[0] = dist(rng);
    pos[1] = dist(rng);
    pos[2] = dist(rng);
  }

  // Allocate storage for the results
  std::vector<std::vector<HBTInt>> ngb_idx(Nsearch);

  // Carry out the search
  cout << "  Finding neighbours using octree" << std::endl;
  for(HBTInt i=0; i<Nsearch; i+=1)
    ngb_idx[i] = tree.NearestNeighbours(centre[i], nr_ngb);

  // Now check the results by brute force:
  // There should be no points closer than the neighbours we found.
  cout << "  Checking neighbour distances" << std::endl;
  for(HBTInt i=0; i<Nsearch; i+=1) {

    // Compute distance to the most distant identified neighbour of this centre
    HBTReal max_dist2 = 0;
    for(auto j : ngb_idx[i]) {
      const HBTxyz &cen = centre[i];
      const HBTxyz &ngb_pos = snap.GetComovingPosition(j);
      const HBTReal ngb_r2 = distance_squared(cen, ngb_pos);
      if(ngb_r2 > max_dist2)max_dist2 = ngb_r2;
    }
    verify(max_dist2 > 0);

    // Check distance to all particles:
    // All neighbours should have distances <= max_dist and non-neighbours
    // should have distances >= max_dist.
    for(HBTInt j=0; j<N; j+=1) {

      // Check if this particle is an identified neighbour
      bool found = false;
      for(auto k : ngb_idx[i]) {
        if(k == j)found=true;
      }

      // Check distance to this particle
      const HBTxyz &cen = centre[i];
      const HBTxyz &part_pos = snap.GetComovingPosition(j);
      HBTReal dist2 = distance_squared(cen, part_pos);
      if(found) {
        verify(dist2 <= max_dist2);
      } else {
        verify(dist2 >= max_dist2);
      }
    }
  }
  cout << "  Test done." << std::endl;
}


int main(int argc, char *argv[]) {

  // Set up repeatable RNG
  std::mt19937 rng;
  rng.seed(0);

  // Test neighbour search with and without periodic boundary
  const int nr_reps = 50;
  for(int rep_nr=0; rep_nr<nr_reps; rep_nr+=1)
  {
    const HBTInt N = 1000;
    const HBTInt Nsearch = 100;
    const HBTReal boxsize = 1.0;
    test_neighbour_search(N, Nsearch, boxsize, /* periodic = */ false, rng);
    test_neighbour_search(N, Nsearch, boxsize, /* periodic = */ true, rng);
  }

  // Test neighbour search with few particles
  for(int rep_nr=0; rep_nr<nr_reps; rep_nr+=1)
  {
    const HBTInt N = 10;
    const HBTInt Nsearch = 100;
    const HBTReal boxsize = 1.0;
    test_neighbour_search(N, Nsearch, boxsize, /* periodic = */ false, rng);
    test_neighbour_search(N, Nsearch, boxsize, /* periodic = */ true, rng);
  }

  // Test neighbour search with one particle
  for(int rep_nr=0; rep_nr<nr_reps; rep_nr+=1)
  {
    const HBTInt N = 1;
    const HBTInt Nsearch = 100;
    const HBTReal boxsize = 1.0;
    test_neighbour_search(N, Nsearch, boxsize, /* periodic = */ false, rng);
    test_neighbour_search(N, Nsearch, boxsize, /* periodic = */ true, rng);
  }
  std::cout << "All tests done." << endl;
}
