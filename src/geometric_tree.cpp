#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "config_parser.h"
#include "geometric_tree.h"
#include "mymath.h"

inline void shift_center(const double oldcenter[3], int son, double delta, double newcenter[3])
{
  for (int dim = 0; dim < 3; dim++)
  {
    int bit = get_bit(son, dim);
    if (bit)
      newcenter[dim] = oldcenter[dim] + delta;
    else
      newcenter[dim] = oldcenter[dim] - delta;
  }
}

void GeoTree_t::ProcessNode(HBTInt nodeid, HBTInt nextid, int sonid, HBTInt &mass, double len, const double center[3])
{
  if (nodeid < NumberOfParticles)
  {
    mass++;
    NextnodeFromParticle[nodeid] = nextid;
  }
  else
  {
    if (len >= HBTConfig.TreeNodeResolution) // only divide if above resolution;
    {
      double newcenter[3];
      shift_center(center, sonid, len / 4., newcenter);
      UpdateInternalNodes(nodeid, nextid, len / 2., newcenter);
    }
    else
      UpdateInternalNodes(nodeid, nextid, len,
                          center); // otherwise we don't divide the node seriouly so we don't have finer node length

    mass += Nodes[nodeid].way.mass; // update mass after updating internal nodes
  }
}

inline void GeoTree_t::FillNodeCenter(HBTInt nodeid, const double center[3])
{
  copyXYZ(Nodes[nodeid].way.s, center);
}

void GeoTree_t::UpdateInternalNodes(HBTInt no, HBTInt sib, double len, const double center[3])
{
  HBTInt p, pp, sons[8];
  int j, jj, i;
  HBTInt mass = 0;

  for (j = 0; j < 8; j++)
    sons[j] = Nodes[no].sons[j]; // backup sons
  Nodes[no].way.len = len;
  Nodes[no].way.sibling = sib;
  for (i = 0; sons[i] < 0; i++)
    ; // find first son
  jj = i;
  pp = sons[jj];
  Nodes[no].way.nextnode = pp;
  for (i++; i < 8; i++) // find sons in pairs,ie. find sibling
  {
    if (sons[i] >= 0) // ok, found a sibling
    {
      j = jj;
      p = pp;
      jj = i;
      pp = sons[jj];
      ProcessNode(p, pp, j, mass, len, center);
    }
  }
  ProcessNode(pp, sib, jj, mass, len, center);
  Nodes[no].way.mass = mass;
  FillNodeCenter(no, center);
}

HBTInt GeoTree_t::NearestNeighbour(const HBTxyz &cen)
// return the particle_index of the nearest neighbour
{
  if(NumberOfParticles == 0) {

    // We have zero particles, so there is no neighbour
    return -1;

  } else if(NumberOfParticles == 1) {

    // We have one particle, so it must be the neighbour
    return 0;

  } else {

    // We have multiple particles. Get initial guess for the search radius.
    // Need >= 2 particles for this to work.
    HBTReal number_density_guess = EstimateNumberDensity(cen);
    HBTReal rguess = GuessNeighbourRange(1, number_density_guess);

    // Search for the neighbours
    NearestNeighbourCollector_t collector;
    Search(cen, rguess, collector);
    while (collector.IsEmpty()) // WARNING: dead loop if tree is empty.
      {
        rguess *= 1.26; // double the guess volume
        Search(cen, rguess, collector);
      }
    return collector.Index;
  }
}

std::vector<HBTInt> GeoTree_t::NearestNeighbours(const HBTxyz &cen, HBTInt max_neighbours)
// return the particle_index of (up to) max_neighbours nearest neighbours
{
  // Determine how many neighbours we should find
  HBTInt nr_to_find = (NumberOfParticles < max_neighbours) ? NumberOfParticles : max_neighbours;

  // Class to store indexes of identified neighbours
  NumNearestNeighboursCollector_t collector(max_neighbours);

  if(nr_to_find == 0)
    {
      // Nothing to do in this case
    }
  else if(NumberOfParticles == nr_to_find)
    {
      // Handle the case where we find all of the particles as neighbours
      for(HBTInt pid=0; pid<nr_to_find; pid+=1)
        {
          bool IsPeriodic = HBTConfig.PeriodicBoundaryOn;
          auto &pos = Snapshot->GetComovingPosition(pid);
          double x0 = cen[0], y0 = cen[1], z0 = cen[2];
          double dx = pos[0] - x0;
          if (IsPeriodic)
            dx = NEAREST(dx);
          double dy = pos[1] - y0;
          if (IsPeriodic)
            dy = NEAREST(dy);
          double dz = pos[2] - z0;
          if (IsPeriodic)
            dz = NEAREST(dz);
          double r2 = dx * dx + dy * dy + dz * dz;
          collector.Collect(pid, r2);
        }
    }
  else
    {
      // This case should have been handled above
      assert(NumberOfParticles > 1);

      // Get an initial guess for the search radius from the octree
      HBTReal number_density_guess = EstimateNumberDensity(cen);
      HBTReal rguess = GuessNeighbourRange(max_neighbours, number_density_guess);

      // Search for the neighbours
      Search(cen, rguess, collector);
      while (collector.NumberFound() < nr_to_find)
        {
          rguess *= 1.26; // double the guess volume
          Search(cen, rguess, collector);
        }
    }

  // Make a vector of neighbour indexes
  std::vector<HBTInt> result;
  while(collector.neighbours.size() > 0)
    {
      LocatedParticle_t lp = collector.neighbours.top();
      result.push_back(lp.index);
      collector.neighbours.pop();
    }
  return result;
}

void GeoTree_t::Search(const HBTxyz &searchcenter, HBTReal radius, ParticleCollector_t &collector)
{ /*find a list of particles from the tree, located within radius around searchcenter,
   * and process the particles with collector */
  bool IsPeriodic = HBTConfig.PeriodicBoundaryOn;
  double x0 = searchcenter[0], y0 = searchcenter[1], z0 = searchcenter[2];
  double h2 = radius * radius;

  HBTInt numngb = 0;
  HBTInt node_id = RootNodeId;

  while (node_id >= 0)
  {
    if (node_id < RootNodeId) /* single particle */
    {
      HBTInt pid = node_id;
      node_id = NextnodeFromParticle[node_id];

      auto &pos = Snapshot->GetComovingPosition(pid);
      double dx = pos[0] - x0;
      if (IsPeriodic)
        dx = NEAREST(dx);
      if (dx > radius || dx < -radius)
        continue;

      double dy = pos[1] - y0;
      if (IsPeriodic)
        dy = NEAREST(dy);
      if (dy > radius || dy < -radius)
        continue;

      double dz = pos[2] - z0;
      if (IsPeriodic)
        dz = NEAREST(dz);
      if (dz > radius || dz < -radius)
        continue;

      double r2 = dx * dx + dy * dy + dz * dz;

      if (r2 < h2)
        collector.Collect(pid, r2);
    }
    else
    {
      auto &node = Nodes[node_id];

      node_id = node.way.sibling; /* in case the node can be discarded */
      double rmax = node.way.len / 2.;
      rmax += radius;

      auto &pos = node.way.s;
      double dx = pos[0] - x0;
      if (IsPeriodic)
        dx = NEAREST(dx);
      if (dx > rmax || dx < -rmax)
        continue;

      double dy = pos[1] - y0;
      if (IsPeriodic)
        dy = NEAREST(dy);
      if (dy > rmax || dy < -rmax)
        continue;

      double dz = pos[2] - z0;
      if (IsPeriodic)
        dz = NEAREST(dz);
      if (dz > rmax || dz < -rmax)
        continue;

      node_id = node.way.nextnode; /* ok, we need to open the node */
    }
  }
}

HBTReal GeoTree_t::EstimateNumberDensity(const HBTxyz &searchcenter)
{
  // Find the node containing searchcentre and use it to estimate the local number density
  bool IsPeriodic = HBTConfig.PeriodicBoundaryOn;
  double x0 = searchcenter[0], y0 = searchcenter[1], z0 = searchcenter[2];

  HBTInt numngb = 0;
  HBTInt node_id = RootNodeId;

  // Side length of the root node
  HBTReal length = Nodes[node_id].way.len;
  assert(length > 0.0); // Will fail if the tree has only one particle

  // Mean number density
  HBTReal mean_density = NumberOfParticles / (length*length*length);

  while (node_id >= 0)
  {
    if (node_id < RootNodeId)
    {
      /* single particle, so estimate number density and return */
      return pow(length, -3.0);
    }
    else
    {
      auto &node = Nodes[node_id];

      node_id = node.way.sibling; /* in case the node can be discarded */
      double rmax = node.way.len / 2.;

      auto &pos = node.way.s;
      double dx = pos[0] - x0;
      if (IsPeriodic)
        dx = NEAREST(dx);
      if (dx > rmax || dx < -rmax)
        continue;

      double dy = pos[1] - y0;
      if (IsPeriodic)
        dy = NEAREST(dy);
      if (dy > rmax || dy < -rmax)
        continue;

      double dz = pos[2] - z0;
      if (IsPeriodic)
        dz = NEAREST(dz);
      if (dz > rmax || dz < -rmax)
        continue;

      length = node.way.len;
      node_id = node.way.nextnode; /* ok, we need to open the node */
    }
  }

  // We get here if the search centre is not in any tree node
  return mean_density;
}


double GeoTree_t::SphDensity(const HBTxyz &cen, HBTReal &hguess)
{
  LocatedParticleCollector_t collector(NumNeighbourSPH * 2);
  vector<LocatedParticle_t> &founds = collector.Founds;
  Search(cen, hguess, collector);
  int numngb = founds.size();
  while (numngb < NumNeighbourSPH)
  {
    if (numngb)
      hguess *= pow(1. * NumNeighbourSPH / numngb, 1.0 / 3.0) *
                1.1; // update hguess adaptively, and conservatively to keep it slightly larger
    else             // zero ngb, double hguess
      hguess *= 2.;

    founds.clear();
    Search(cen, hguess, collector);
    numngb = founds.size();
  }

  auto pivot_particle = founds.begin() + NumNeighbourSPH - 1;
  nth_element(founds.begin(), pivot_particle, founds.end(), CompLocatedDistance);
  double h = sqrt(pivot_particle->d2);
  // 	h=sqrtf(h);
  hguess = h * 1.01;
  double hinv3 = 1.0 / (h * h * h);

  double rho = 0.;
  for (auto it = founds.begin(); it <= pivot_particle; ++it)
  {
    double r = sqrt(it->d2);
    double u = r / h, wk;

    if (u < 0.5)
      wk = hinv3 * (2.546479089470 + 15.278874536822 * (u - 1) * u * u);
    else
      wk = hinv3 * 5.092958178941 * (1.0 - u) * (1.0 - u) * (1.0 - u);

    rho += wk;
  }
  return rho;
}

template class OctTree_t<GeoTreeCell_t>; // to wake up the functions for this type; trick!
