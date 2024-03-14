#include "linkedlist_base.h"
#include "mymath.h"

HBTReal LinkedlistBase_t::Distance2(const HBTxyz &x, const HBTxyz &y)
{
  HBTxyz dx;
  dx[0] = x[0] - y[0];
  dx[1] = x[1] - y[1];
  dx[2] = x[2] - y[2];
  if (PeriodicBoundary)
  {
#define _NEAREST(x) (((x) > BoxHalf) ? ((x) - BoxSize) : (((x) < -BoxHalf) ? ((x) + BoxSize) : (x)))
    dx[0] = _NEAREST(dx[0]);
    dx[1] = _NEAREST(dx[1]);
    dx[2] = _NEAREST(dx[2]);
#undef _NEAREST
  }
  return (dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);
}
void LinkedlistBase_t::init(int ndiv, PositionData_t *data, HBTReal boxsize, bool periodic)
{ // initialize the variables, some book-keeping, and allocate memory
  NDiv = ndiv;
  NDiv2 = NDiv * NDiv;
  Particles = data;
  PositionData_t &particles = *Particles;
  HBTInt np = particles.size();
  BoxSize = boxsize;
  PeriodicBoundary = periodic;
  if (BoxSize == 0.)
    PeriodicBoundary = false; // only effective when boxsize is specified
  BoxHalf = BoxSize / 2.;

  HBTInt i, j;
  // cout<<"creating linked list for "<<particles.size()<<" particles..."<<endl;
  /*determining enclosing cube*/
  if (BoxSize)
  {
    for (i = 0; i < 3; i++)
    {
      Range[i][0] = 0.;
      Range[i][1] = BoxSize;
    }
    for (j = 0; j < 3; j++)
      Step[j] = BoxSize / NDiv;
  }
  else
  {
    for (i = 0; i < 3; i++)
      for (j = 0; j < 2; j++)
        Range[i][j] = particles[0][j];
    for (i = 1; i < np; i++)
      for (j = 0; j < 3; j++)
      {
        auto x = particles[i][j];
        if (x < Range[j][0])
          Range[j][0] = x;
        else if (x > Range[j][1])
          Range[j][1] = x;
      }
    for (j = 0; j < 3; j++)
      Step[j] = (Range[j][1] - Range[j][0]) / NDiv;
  }

  HOC.assign(NDiv2 * NDiv, -1);
  List.resize(np);
}
void LinkedlistBase_t::build(int ndiv, PositionData_t *data, HBTReal boxsize, bool periodic)
{
  init(ndiv, data, boxsize, periodic);

  HBTInt i, j, grid[3];
  HBTInt ind;

  PositionData_t &particles = *Particles;
  HBTInt np = particles.size();
  for (i = 0; i < np; i++)
  {
    for (j = 0; j < 3; j++)
    {
      grid[j] = floor((particles[i][j] - Range[j][0]) / Step[j]);
      grid[j] = FixGridId(grid[j]);
    }
    ind = Sub2Ind(grid[0], grid[1], grid[2]);
    List[i] = HOC[ind];
    HOC[ind] = i; /*use hoc[ind] as swap varible to temporarily
                                  store last ll index, and finally the head*/
  }
}
void LinkedlistBase_t::SearchShell(HBTReal rmin, HBTReal rmax, const HBTxyz &searchcenter,
                                   ParticleCollector_t &collector)
{ /*search in [rmin, rmax] radial range (inclusive), and append results to founds
   */
  PositionData_t &particles = *Particles;
  HBTReal rmin2 = rmin * rmin, rmax2 = rmax * rmax;
  int i, j, k, subbox_grid[3][2];

  if (rmin <= 0)
  {
    SearchSphere(rmax, searchcenter, collector);
    return;
  }

  for (i = 0; i < 3; i++)
  {
    subbox_grid[i][0] = floor((searchcenter[i] - rmax - Range[i][0]) / Step[i]);
    subbox_grid[i][1] = floor((searchcenter[i] + rmax - Range[i][0]) / Step[i]);
    if (!PeriodicBoundary)
    { // do not fix if periodic, since the search sphere is allowed to overflow the box in periodic case.
      subbox_grid[i][0] = FixGridId(subbox_grid[i][0]);
      subbox_grid[i][1] = FixGridId(subbox_grid[i][1]);
    }
  }
  for (i = subbox_grid[0][0]; i < subbox_grid[0][1] + 1; i++)
    for (j = subbox_grid[1][0]; j < subbox_grid[1][1] + 1; j++)
      for (k = subbox_grid[2][0]; k < subbox_grid[2][1] + 1; k++)
      {
        HBTInt pid = GetHOCSafe(i, j, k); // in case the grid-id is out of box, in the periodic case
        while (pid >= 0)
        {
          HBTReal dr2 = Distance2(particles[pid], searchcenter);
          if (dr2 < rmax2 && dr2 > rmin2)
            collector.Collect(pid, dr2);
          pid = List[pid];
        }
      }
}
void LinkedlistBase_t::SearchSphere(HBTReal radius, const HBTxyz &searchcenter, ParticleCollector_t &collector)
{
  PositionData_t &particles = *Particles;
  HBTReal x0 = searchcenter[0], y0 = searchcenter[1], z0 = searchcenter[2];
  HBTReal radius2 = radius * radius;
  int i, j, k, subbox_grid[3][2];

  for (i = 0; i < 3; i++)
  {
    subbox_grid[i][0] = floor((searchcenter[i] - radius - Range[i][0]) / Step[i]);
    subbox_grid[i][1] = floor((searchcenter[i] + radius - Range[i][0]) / Step[i]);
    if (!PeriodicBoundary)
    { // do not fix if periodic, since the search sphere is allowed to overflow the box in periodic case.
      subbox_grid[i][0] = FixGridId(subbox_grid[i][0]);
      subbox_grid[i][1] = FixGridId(subbox_grid[i][1]);
    }
  }
  for (i = subbox_grid[0][0]; i < subbox_grid[0][1] + 1; i++)
    for (j = subbox_grid[1][0]; j < subbox_grid[1][1] + 1; j++)
      for (k = subbox_grid[2][0]; k < subbox_grid[2][1] + 1; k++)
      {
        HBTInt pid = GetHOCSafe(i, j, k); // in case the grid-id is out of box, in the periodic case
        while (pid >= 0)
        {
          HBTInt p = pid;
          pid = List[p];

          auto &pos = particles[p];
          HBTReal dx = pos[0] - x0;
          if (PeriodicBoundary)
            dx = NEAREST(dx);
          if (dx > radius || dx < -radius)
            continue;

          HBTReal dy = pos[1] - y0;
          if (PeriodicBoundary)
            dy = NEAREST(dy);
          if (dy > radius || dy < -radius)
            continue;

          HBTReal dz = pos[2] - z0;
          if (PeriodicBoundary)
            dz = NEAREST(dz);
          if (dz > radius || dz < -radius)
            continue;

          HBTReal r2 = dx * dx + dy * dy + dz * dz;

          if (r2 < radius2)
            collector.Collect(p, r2);
        }
      }
}

void LinkedlistBase_t::SearchCylinder(HBTReal radius_z, HBTReal radius_p, const HBTxyz &searchcenter,
                                      ParticleCollector_t &collector)
{
  PositionData_t &particles = *Particles;
  HBTReal x0 = searchcenter[0], y0 = searchcenter[1], z0 = searchcenter[2];
  HBTReal radius2 = radius_p * radius_p;
  int i, j, k, subbox_grid[3][2];
  HBTReal radvec[3] = {radius_p, radius_p, radius_z};

  for (i = 0; i < 3; i++)
  {
    subbox_grid[i][0] = floor((searchcenter[i] - radvec[i] - Range[i][0]) / Step[i]);
    subbox_grid[i][1] = floor((searchcenter[i] + radvec[i] - Range[i][0]) / Step[i]);
    if (!PeriodicBoundary)
    { // do not fix if periodic, since the search sphere is allowed to overflow the box in periodic case.
      subbox_grid[i][0] = FixGridId(subbox_grid[i][0]);
      subbox_grid[i][1] = FixGridId(subbox_grid[i][1]);
    }
  }

  for (i = subbox_grid[0][0]; i < subbox_grid[0][1] + 1; i++)
    for (j = subbox_grid[1][0]; j < subbox_grid[1][1] + 1; j++)
      for (k = subbox_grid[2][0]; k < subbox_grid[2][1] + 1; k++)
      {
        HBTInt pid = GetHOCSafe(i, j, k); // in case the grid-id is out of box, in the periodic case
        while (pid >= 0)
        {
          HBTInt p = pid;
          pid = List[p];

          auto &pos = particles[p];
          HBTReal dx = pos[0] - x0;
          if (PeriodicBoundary)
            dx = NEAREST(dx);
          if (dx > radius_p || dx < -radius_p)
            continue;

          HBTReal dy = pos[1] - y0;
          if (PeriodicBoundary)
            dy = NEAREST(dy);
          if (dy > radius_p || dy < -radius_p)
            continue;

          HBTReal dz = pos[2] - z0;
          if (PeriodicBoundary)
            dz = NEAREST(dz);
          if (dz > radius_z || dz < -radius_z)
            continue;

          HBTReal r2 = dx * dx + dy * dy;

          if (r2 < radius2)
            collector.Collect(p, r2);
        }
      }
}

HBTInt LinkedlistBase_t::TagFriendsOfFriends(HBTInt seed, HBTInt grpid, vector<HBTInt> &group_tags, HBTReal LinkLength)
/*tag all the particles that are linked to seed with grpid
 * Note if system stack size is too small, this recursive routine may crash
 * in that case you should set:  ulimit -s unlimited  (bash) before running.
 **/
{
  PositionData_t &particles = *Particles;
  auto &searchcenter = particles[seed];
  HBTReal x0 = searchcenter[0], y0 = searchcenter[1], z0 = searchcenter[2];
  HBTReal LinkLength2 = LinkLength * LinkLength;
  int i, j, k, subbox_grid[3][2];

  vector<HBTInt> friends;

  for (i = 0; i < 3; i++)
  {
    subbox_grid[i][0] = floor((searchcenter[i] - LinkLength - Range[i][0]) / Step[i]);
    subbox_grid[i][1] = floor((searchcenter[i] + LinkLength - Range[i][0]) / Step[i]);
    if (!PeriodicBoundary)
    { // do not fix if periodic, since the search sphere is allowed to overflow the box in periodic case.
      subbox_grid[i][0] = FixGridId(subbox_grid[i][0]);
      subbox_grid[i][1] = FixGridId(subbox_grid[i][1]);
    }
  }
  for (i = subbox_grid[0][0]; i < subbox_grid[0][1] + 1; i++)
    for (j = subbox_grid[1][0]; j < subbox_grid[1][1] + 1; j++)
      for (k = subbox_grid[2][0]; k < subbox_grid[2][1] + 1; k++)
      {
        HBTInt pid = GetHOCSafe(i, j, k); // in case the grid-id is out of box, in the periodic case
        while (pid >= 0)
        {
          HBTInt p = pid;
          pid = List[p];
          if (group_tags[p] >= 0) // already tagged
            continue;

          auto &pos = particles[p];
          HBTReal dx = pos[0] - x0;
          if (PeriodicBoundary)
            dx = NEAREST(dx);
          if (dx > LinkLength || dx < -LinkLength)
            continue;

          HBTReal dy = pos[1] - y0;
          if (PeriodicBoundary)
            dy = NEAREST(dy);
          if (dy > LinkLength || dy < -LinkLength)
            continue;

          HBTReal dz = pos[2] - z0;
          if (PeriodicBoundary)
            dz = NEAREST(dz);
          if (dz > LinkLength || dz < -LinkLength)
            continue;

          HBTReal r2 = dx * dx + dy * dy + dz * dz;

          if (r2 < LinkLength2)
          {
            group_tags[p] = grpid;
            friends.emplace_back(p);
          }
        }
      }

  HBTInt nfriends = friends.size();
  for (auto &&p : friends)
    nfriends += TagFriendsOfFriends(p, grpid, group_tags, LinkLength);

  return nfriends; // excluding the seed particle
}
