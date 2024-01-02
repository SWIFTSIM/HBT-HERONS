#ifndef GADGET_GROUP_IO_INCLUDED
#define GADGET_GROUP_IO_INCLUDED

#include "../halo.h"
#include "../mpi_wrapper.h"

namespace GadgetGroup
{
struct GroupV4Header_t
{
  int Ngroups;
  int Nsubgroups;
  int Nids;
  int TotNgroups;
  int TotNsubgroups;
  int TotNids;
  int num_files; // long? no, but padding may exist.
  double time;
  double redshift;
  double HubbleParam;
  double BoxSize;
  double Omega0;
  double OmegaLambda;
  int flag_doubleprecision; // long? no, but padding may exist
};

extern void Load(MpiWorker_t &world, int SnapshotId, vector<Halo_t> &Halos);
extern bool IsGadgetGroup(const string &GroupFileFormat);

} // namespace GadgetGroup

#endif