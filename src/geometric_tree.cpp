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

HBTInt GeoTree_t::NearestNeighbour(const HBTxyz &cen, HBTReal rguess)
// return the particle_index of the nearest neighbour
{
  NearestNeighbourCollector_t collector;
  Search(cen, rguess, collector);
  while (collector.IsEmpty()) // WARNING: dead loop if tree is empty.
  {
    rguess *= 1.26; // double the guess volume
    Search(cen, rguess, collector);
  }
  return collector.Index;
}

std::vector<HBTInt> GeoTree_t::NearestNeighbours(const HBTxyz &cen, HBTReal rguess, HBTInt max_neighbours)
// return the particle_index of (up to) max_neighbours nearest neighbours
{
  NumNearestNeighboursCollector_t collector(max_neighbours);
  HBTInt nr_to_find = (NumberOfParticles < max_neighbours) ? NumberOfParticles : max_neighbours;
  Search(cen, rguess, collector);
  while (collector.NumberFound() < nr_to_find)
  {
    rguess *= 1.26; // double the guess volume
    Search(cen, rguess, collector);
  }
  // Make a vector of neighbour indexes
  std::vector<HBTInt> result;
  while(collector.neighbours.size() > 0) {
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
