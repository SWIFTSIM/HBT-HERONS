#include <vector>

#include "config_parser.h"
#include "datatypes.h"
#include "verify.h"

// Defined in src/subhalo_tracking.cpp
IdRank_t DecideLocalHostId(std::vector<IdRank_t>::const_iterator particle_hosts);

//
// Test the host decision used when the tracers of a subhalo are spread over
// several MPI ranks. Each entry is the (local FoF index, rank) of one tracer,
// ordered from most to least bound. FoF indices are only unique within a rank,
// so FoFs with the same index on different ranks must be scored separately.
//
// Tracer weights are 1 / (1 + sqrt(i)) for i = 0, ..., 9:
// 1.000 0.500 0.414 0.366 0.333 0.309 0.290 0.274 0.261 0.240
//
int main(int argc, char *argv[])
{
  HBTConfig.NumTracerHostFinding = 10;

  // FoF 3 on rank 0 has the highest score (1.50). FoF 5 on rank 1 (1.42) and
  // FoF 5 on rank 2 (1.07) are different groups and must not be combined.
  {
    std::vector<IdRank_t> hosts{{3, 0}, {3, 0}, {5, 1}, {5, 1}, {5, 1},
                                {5, 1}, {5, 2}, {5, 2}, {5, 2}, {5, 2}};
    IdRank_t host = DecideLocalHostId(hosts.cbegin());
    verify(host.Id == 3);
    verify(host.Rank == 0);
  }

  // FoF 5 on rank 1 has most of the tracers, so the subhalo must go to rank 1,
  // not to rank 2 where the least bound tracers are.
  {
    std::vector<IdRank_t> hosts{{5, 1}, {5, 1}, {5, 1}, {5, 1}, {5, 1},
                                {5, 1}, {5, 2}, {5, 2}, {5, 2}, {5, 2}};
    IdRank_t host = DecideLocalHostId(hosts.cbegin());
    verify(host.Id == 5);
    verify(host.Rank == 1);
  }

  // Hostless tracers (-1) on several ranks count as a single candidate, placed
  // on the rank of the most bound hostless tracer.
  {
    std::vector<IdRank_t> hosts{{2, 0},  {-1, 3}, {-1, 1}, {-1, 1}, {-1, 2},
                                {-1, 1}, {-1, 2}, {-1, 1}, {2, 0},  {2, 0}};
    IdRank_t host = DecideLocalHostId(hosts.cbegin());
    verify(host.Id == -1);
    verify(host.Rank == 3);
  }

  // An orphan has a single tracer; the other entries are not valid (-2).
  {
    std::vector<IdRank_t> hosts{{4, 2},  {-2, 0}, {-2, 0}, {-2, 0}, {-2, 0},
                                {-2, 0}, {-2, 0}, {-2, 0}, {-2, 0}, {-2, 0}};
    IdRank_t host = DecideLocalHostId(hosts.cbegin());
    verify(host.Id == 4);
    verify(host.Rank == 2);
  }

  return 0;
}
