using namespace std;
#include <iostream>
#include <numeric>
// #include <iomanip>
#include <sstream>
#include <string>
#include <typeinfo>
#include <assert.h>
#include <cstdlib>
#include <cstdio>
#include <chrono>
#include <list>
#include <stdexcept>

#include "../../snapshot.h"
#include "../../mymath.h"
#include "../../hdf_wrapper.h"
#include "gadget4_io.h"
#include "../comms/exchange_and_merge.h"

/* Checks whether a file is accessible. */
inline bool is_it_valid(const std::string& name) {
  struct stat buffer;
  return (stat (name.c_str(), &buffer) == 0);
}

namespace Gadget4Reader
{

/* Creates the Header struct data type for MPI communication */
void create_Gadget4Header_MPI_type(MPI_Datatype &dtype)
{
  Gadget4Header_t p;
#define NumAttr 13
  MPI_Datatype oldtypes[NumAttr];
  int blockcounts[NumAttr];
  MPI_Aint offsets[NumAttr], origin, extent;

  MPI_Get_address(&p, &origin);
  MPI_Get_address((&p) + 1, &extent); // to get the extent of s
  extent -= origin;

  int i = 0;
#define RegisterAttr(x, type, count)                                                                                   \
  {                                                                                                                    \
    MPI_Get_address(&(p.x), offsets + i);                                                                              \
    offsets[i] -= origin;                                                                                              \
    oldtypes[i] = type;                                                                                                \
    blockcounts[i] = count;                                                                                            \
    i++;                                                                                                               \
  }
  RegisterAttr(NumberOfFiles, MPI_INT, 1);
  RegisterAttr(BoxSize, MPI_DOUBLE, 1);
  RegisterAttr(ScaleFactor, MPI_DOUBLE, 1);
  RegisterAttr(OmegaM0, MPI_DOUBLE, 1);
  RegisterAttr(OmegaLambda0, MPI_DOUBLE, 1);
  RegisterAttr(mass, MPI_DOUBLE, TypeMax);
  RegisterAttr(npart[0], MPI_INT, TypeMax);
  RegisterAttr(npartTotal[0], MPI_HBT_INT, TypeMax);

  RegisterAttr(MassInMsunh, MPI_DOUBLE, 1);
  RegisterAttr(LengthInMpch, MPI_DOUBLE, 1);
  RegisterAttr(VelInKmS, MPI_DOUBLE, 1);

  RegisterAttr(DM_comoving_softening, MPI_DOUBLE, 1);
  RegisterAttr(DM_maximum_physical_softening, MPI_DOUBLE, 1);

#undef RegisterAttr
          assert(i <= NumAttr);

  MPI_Type_create_struct(i, blockcounts, offsets, oldtypes, &dtype);
  MPI_Type_create_resized(dtype, (MPI_Aint)0, extent, &dtype);
  MPI_Type_commit(&dtype);
#undef NumAttr
}

void Gadget4Reader_t::SetSnapshot(int snapshotId)
{
  SnapshotId = snapshotId;
  if (HBTConfig.SnapshotNameList.empty())
  {
    stringstream formatter;
    if (HBTConfig.SnapshotDirBase.length() > 0)
      formatter << HBTConfig.SnapshotDirBase << "_" << setw(3) << setfill('0') << snapshotId << "/";
    formatter << HBTConfig.SnapshotFileBase << "_" << setw(3) << setfill('0') << snapshotId;
    SnapshotName = formatter.str();
  }
  else
    SnapshotName = HBTConfig.SnapshotNameList[snapshotId];
}

/* Generates the path to the snapshot (sub)file */
void Gadget4Reader_t::GetFileName(int ifile, string &filename)
{
  stringstream formatter;

  if (HBTConfig.SnapshotDirBase.length() > 0)
   formatter << HBTConfig.SnapshotPath << "/" << SnapshotName << "." << ifile << ".hdf5";
  else
    formatter << HBTConfig.SnapshotPath << "/" << SnapshotName << ".hdf5";

  filename = formatter.str();
}

void Gadget4Reader_t::GetGroupFileName(int ifile, string &filename)
{
  string snap_idname = SnapshotName.substr(SnapshotName.size() - 3); // last 3 chars
  stringstream formatter;

  /* Halo catalogues have the same number of files as snapshots. Within each block
  /* we try two possible combinations, since the name of the FoF catalogue
   * depends on whether GADGET-4 only ran FoF or if it also did SUBFIND/SUBFIND-HBT. */
  if (HBTConfig.SnapshotDirBase.length() > 0)
  {
    /* If only FoF was run. */
    formatter << HBTConfig.HaloPath << "groups_" << setw(3) << setfill('0') << SnapshotId << "/fof_tab_" << snap_idname << "." << ifile << ".hdf5";
    filename = formatter.str();
    bool is_valid_file = is_it_valid(filename);

    if(is_valid_file)
      return;

    std::stringstream().swap(formatter); // Empty stringstream

    /* Subfind or Subfind-HBT was also run. */
    formatter << HBTConfig.HaloPath << "groups_" << setw(3) << setfill('0') << SnapshotId << "/fof_subhalo_tab_" << snap_idname << "." << ifile << ".hdf5";
    filename = formatter.str();
  }
  else
  {
    /* If only FoF was run. */
    formatter << HBTConfig.HaloPath << "/fof_tab_" << snap_idname << ".hdf5";
    filename = formatter.str();
    bool is_valid_file = is_it_valid(filename);

    if(is_valid_file)
      return;

    std::stringstream().swap(formatter); // Empty stringstream

    /* Subfind or Subfind-HBT was also run. */
    formatter << HBTConfig.HaloPath << "/fof_subhalo_tab_" << snap_idname << ".hdf5";
    filename = formatter.str();
  }
  return;
}

void Gadget4Reader_t::ReadHeader(int ifile, Gadget4Header_t &header)
{
  string filename;
  GetFileName(ifile, filename);
  hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

  ReadAttribute(file, "Header", "NumFilesPerSnapshot", H5T_NATIVE_INT, &header.NumberOfFiles);
  ReadAttribute(file, "Header", "BoxSize", H5T_NATIVE_DOUBLE, &header.BoxSize);
  ReadAttribute(file, "Header", "Time", H5T_NATIVE_DOUBLE, &header.ScaleFactor);
  ReadAttribute(file, "Parameters", "Omega0", H5T_NATIVE_DOUBLE, &header.OmegaM0);
  ReadAttribute(file, "Parameters", "OmegaLambda", H5T_NATIVE_DOUBLE, &header.OmegaLambda0);
  for (int i = 0; i < TypeMax; i++) // initialize
  {
    header.mass[i] = 0.;
    header.npart[i] = 0;
    header.npartTotal[i] = 0;
  }
  ReadAttribute(file, "Header", "MassTable", H5T_NATIVE_DOUBLE, header.mass);
  ReadAttribute(file, "Header", "NumPart_ThisFile", H5T_NATIVE_INT, header.npart);
  ReadAttribute(file, "Header", "NumPart_Total", H5T_HBTInt, header.npartTotal);

  /* Read unit system used by GADGET-4 */
  double length_cgs;
  ReadAttribute(file, "Parameters", "UnitLength_in_cm", H5T_NATIVE_DOUBLE, &length_cgs);
  double mass_cgs;
  ReadAttribute(file, "Parameters", "UnitMass_in_g", H5T_NATIVE_DOUBLE, &mass_cgs);
  double velocity_cgs;
  ReadAttribute(file, "Parameters", "UnitVelocity_in_cm_per_s", H5T_NATIVE_DOUBLE, &velocity_cgs);
  cout << length_cgs << " " << mass_cgs << " " << velocity_cgs << " " << endl;

  /* To later tell HBT-HERONS the units of the snapshot */
  Header.MassInMsunh  = mass_cgs / 1.98841e33;
  Header.LengthInMpch = length_cgs / 3.08567758e24;
  Header.VelInKmS = velocity_cgs / 1e5;

  /* For now we assume that all particles have the same softening as PartType1. */
  // TODO: generalise to allow for different softenings per particle type.
  int DM_softening_class;
  ReadAttribute(file, "Parameters", "SofteningClassOfPartType1", H5T_NATIVE_INT, &DM_softening_class);

  stringstream DM_comoving_softening_dataset;
  DM_comoving_softening_dataset << "SofteningComovingClass" << DM_softening_class;
  ReadAttribute(file, "Parameters", DM_comoving_softening_dataset.str().c_str(), H5T_NATIVE_DOUBLE, &Header.DM_comoving_softening);

  stringstream DM_max_softening_dataset;
  DM_max_softening_dataset << "SofteningComovingClass" << DM_softening_class;
  ReadAttribute(file, "Parameters", DM_max_softening_dataset.str().c_str(), H5T_NATIVE_DOUBLE, &Header.DM_maximum_physical_softening);

  H5Fclose(file);
}

int Gadget4Reader_t::ReadGroupFileCounts(int ifile)
{
  int nfiles;
  string filename;
  GetGroupFileName(ifile, filename);
  hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  ReadAttribute(file, "Header", "NumFiles", H5T_NATIVE_INT, &nfiles);
  H5Fclose(file);

  return nfiles;
}

HBTInt Gadget4Reader_t::ReadGroupFileTotParticles(int ifile)
{ // total np in groups from all files
  HBTInt np_tot;
  string filename;
  GetGroupFileName(ifile, filename);
  hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  ReadAttribute(file, "Header", "Nids_Total", H5T_HBTInt, &np_tot);
  H5Fclose(file);

  return np_tot;
}

void Gadget4Reader_t::GetParticleCountInFile(hid_t file, int np[])
{
  ReadAttribute(file, "Header", "NumPart_ThisFile", H5T_NATIVE_INT, np);
#ifdef DM_ONLY
  for (int i = 0; i < TypeMax; i++)
    if (i != TypeDM)
      np[i] = 0;
#endif
}
HBTInt Gadget4Reader_t::CompileFileOffsets(int nfiles)
{
  HBTInt offset = 0;
  np_file.reserve(nfiles);
  offset_file.reserve(nfiles);
  for (int ifile = 0; ifile < nfiles; ifile++)
  {
    offset_file.push_back(offset);

    int np_this[TypeMax];
    string filename;
    GetFileName(ifile, filename);
    hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    GetParticleCountInFile(file, np_this);
    H5Fclose(file);
    HBTInt np = accumulate(begin(np_this), end(np_this), (HBTInt)0);

    np_file.push_back(np);
    offset += np;
  }
  return offset;
}

HBTInt Gadget4Reader_t::CompileGroupFileOffsets(vector<HBTInt> &nhalo_per_groupfile,
                                                vector<HBTInt> &offsethalo_per_groupfile)
{
  HBTInt offset = 0, nhalo_total;
  int nfiles = nhalo_per_groupfile.size();
  for (int ifile = 0; ifile < nfiles; ifile++)
  {
    offsethalo_per_groupfile[ifile] = offset;

    HBTInt nhalo_this;
    string filename;
    GetGroupFileName(ifile, filename);
    hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    ReadAttribute(file, "Header", "Ngroups_ThisFile", H5T_HBTInt, &nhalo_this);
    if (ifile == 0)
      ReadAttribute(file, "Header", "Ngroups_Total", H5T_HBTInt, &nhalo_total);
    H5Fclose(file);

    nhalo_per_groupfile[ifile] = nhalo_this;
    offset += nhalo_this;
  }
  assert(nhalo_total == offset);
  return offset;
}

static void check_id_size(hid_t loc)
{
  hid_t dset = H5Dopen2(loc, "ParticleIDs", H5P_DEFAULT);
  hid_t dtype = H5Dget_type(dset);
  size_t ParticleIDStorageSize = H5Tget_size(dtype);
  assert(sizeof(HBTInt) >= ParticleIDStorageSize); // use HBTi8 or HBTdouble if you need long int for id
  H5Tclose(dtype);
  H5Dclose(dset);
}

/* Reads the snapshot particle information. Note that this method does not load
 * the HostId of particles, since that information needs to be obtained from the
 * halo lengths provided in the GADGET4 group catalogues. */
void Gadget4Reader_t::ReadSnapshot(int ifile, Particle_t *ParticlesInFile)
{
  string filename;
  GetFileName(ifile, filename);
  hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  vector<int> np_this(TypeMax);
  vector<HBTInt> offset_this(TypeMax);
  GetParticleCountInFile(file, np_this.data());
  CompileOffsets(np_this, offset_this);

  HBTReal vunit = sqrt(Header.ScaleFactor);
  HBTReal boxsize = Header.BoxSize;
  for (int itype = 0; itype < TypeMax; itype++)
  {
    int np = np_this[itype];
    if (np == 0)
      continue;
    auto ParticlesThisType = ParticlesInFile + offset_this[itype];
    stringstream grpname;
    grpname << "PartType" << itype;
    if (!H5Lexists(file, grpname.str().c_str(), H5P_DEFAULT))
      continue;
    hid_t particle_data = H5Gopen2(file, grpname.str().c_str(), H5P_DEFAULT);
    // 	if(particle_data<0) continue; //skip non-existing type

    check_id_size(particle_data);

    { // read position
      vector<HBTxyz> x(np);
      ReadDataset(particle_data, "Coordinates", H5T_HBTReal, x.data());
      if (HBTConfig.PeriodicBoundaryOn)
      {
#pragma omp parallel for
        for (int i = 0; i < np; i++)
          for (int j = 0; j < 3; j++)
            x[i][j] = position_modulus(x[i][j], boxsize);
      }
#pragma omp parallel for
      for (int i = 0; i < np; i++)
        copyHBTxyz(ParticlesThisType[i].ComovingPosition, x[i]);
    }

    { // velocity
      vector<HBTxyz> v(np);
      if (H5Lexists(particle_data, "Velocities", H5P_DEFAULT))
        ReadDataset(particle_data, "Velocities", H5T_HBTReal, v.data());
      else
        ReadDataset(particle_data, "Velocity", H5T_HBTReal, v.data());
#pragma omp parallel for
      for (int i = 0; i < np; i++)
        for (int j = 0; j < 3; j++)
          ParticlesThisType[i].PhysicalVelocity[j] = v[i][j] * vunit;
    }

    { // id
      vector<HBTInt> id(np);
      ReadDataset(particle_data, "ParticleIDs", H5T_HBTInt, id.data());
#pragma omp parallel for
      for (int i = 0; i < np; i++)
        ParticlesThisType[i].Id = id[i];
    }

    // mass
    if (Header.mass[itype] == 0)
    {
      vector<HBTReal> m(np);
      if (H5Lexists(particle_data, "Masses", H5P_DEFAULT))
        ReadDataset(particle_data, "Masses", H5T_HBTReal, m.data());
      else
        ReadDataset(particle_data, "Mass", H5T_HBTReal, m.data());
#pragma omp parallel for
      for (int i = 0; i < np; i++)
        ParticlesThisType[i].Mass = m[i];
    }
    else
    {
#pragma omp parallel for
      for (int i = 0; i < np; i++)
        ParticlesThisType[i].Mass = Header.mass[itype];
    }

#ifndef DM_ONLY
    // internal energy
#ifdef HAS_THERMAL_ENERGY
    if (itype == 0)
    {
      vector<HBTReal> u(np);
      ReadDataset(particle_data, "InternalEnergy", H5T_HBTReal, u.data());
#pragma omp parallel for
      for (int i = 0; i < np; i++)
        ParticlesThisType[i].InternalEnergy = u[i];
    }
/*	else
    {
      for(int i=0;i<np;i++)
        ParticlesThisType[i].InternalEnergy=0; //necessary? maybe default initialized?
    }*/
#endif

    { // type
      ParticleType_t t = static_cast<ParticleType_t>(itype);
#pragma omp parallel for
      for (int i = 0; i < np; i++)
        ParticlesThisType[i].Type = t;
    }
#endif

    H5Gclose(particle_data);
  }

  H5Fclose(file);
}

void Gadget4Reader_t::ReadGroupLen(int ifile, HBTInt *buf)
{
  string filename;
  GetGroupFileName(ifile, filename);
  hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  HBTInt ngroups;
  ReadAttribute(file, "Header", "Ngroups_ThisFile", H5T_HBTInt, &ngroups);
  if (ngroups > 0)
    ReadDataset(file, "Group/GroupLen", H5T_HBTInt, buf);
  H5Fclose(file);
}

void Gadget4Reader_t::LoadSnapshot(MpiWorker_t &world, int snapshotId, vector<Particle_t> &Particles,
                                   Cosmology_t &Cosmology)
{
  SetSnapshot(snapshotId);

  if (world.rank() == root_node)
  {
    ReadHeader(0, Header);
    CompileFileOffsets(Header.NumberOfFiles);
  }

  MPI_Bcast(&Header, 1, MPI_Gadget4Header_t, root_node, world.Communicator);
  world.SyncContainer(np_file, MPI_HBT_INT, root_node);
  world.SyncContainer(offset_file, MPI_HBT_INT, root_node);

  Cosmology.Set(Header.ScaleFactor, Header.OmegaM0, Header.OmegaLambda0);

  /* Assign the box size read in from the Header */
  HBTConfig.BoxSize = Header.BoxSize;
  HBTConfig.BoxHalf = HBTConfig.BoxSize / 2;

  /* Use the softening values we read in from the Header */
  HBTConfig.SofteningHalo = Header.DM_comoving_softening;
  HBTConfig.MaxPhysicalSofteningHalo = Header.DM_maximum_physical_softening;

  /* Update the tree parameters based on the softenings */
  HBTConfig.TreeNodeResolution = HBTConfig.SofteningHalo * 0.1;
  HBTConfig.TreeNodeResolutionHalf = HBTConfig.TreeNodeResolution / 2.;

  /* Update the units */
  HBTConfig.MassInMsunh = Header.MassInMsunh;
  HBTConfig.LengthInMpch = Header.LengthInMpch;
  HBTConfig.VelInKmS = Header.VelInKmS;

  /* Update the gravitational constant and H based on the units loaded from
   * GADGET-4 */
  PhysicalConst::G =
    43.0071 * (HBTConfig.MassInMsunh / 1e10) / HBTConfig.VelInKmS / HBTConfig.VelInKmS / HBTConfig.LengthInMpch;
  PhysicalConst::H0 = 100. * (1. / HBTConfig.VelInKmS) / (1. / HBTConfig.LengthInMpch);


  /* Read physical properties of particles and their particle type. */
  LoadParticleProperties(world, Particles);

  /* Read host halo of each particle. GADGET4 has particle and group information
   * separate, hence the need for this separate function call. */
  LoadParticleHosts(world, Particles);
}

typedef vector<HBTInt> CountBuffer_t;

/* Loads halo particle sizes from all files in parallel and gathers to root MPI
 * rank. */
void Gadget4Reader_t::LoadHaloSizes(MpiWorker_t &world)
{

  /* Read file metadata */
  int FileCounts;
  vector<HBTInt> nhalo_per_groupfile;
  vector<HBTInt> offsethalo_per_groupfile;
  HBTInt NhaloTotal;
  if (world.rank() == root_node)
  {
    FileCounts = ReadGroupFileCounts(0);
    nhalo_per_groupfile.resize(FileCounts);
    offsethalo_per_groupfile.resize(FileCounts);
    NhaloTotal = CompileGroupFileOffsets(nhalo_per_groupfile, offsethalo_per_groupfile);
    TotNumPartInGroups = ReadGroupFileTotParticles(0);
  }
  world.SyncAtom(FileCounts, MPI_INT, root_node);
  world.SyncAtom(NhaloTotal, MPI_HBT_INT, root_node);
  world.SyncContainer(nhalo_per_groupfile, MPI_HBT_INT, root_node);
  world.SyncContainer(offsethalo_per_groupfile, MPI_HBT_INT, root_node);
  world.SyncAtom(TotNumPartInGroups, MPI_HBT_INT, root_node);

  /* Read halo sizes in parallel. */
  HBTInt nfiles_skip, nfiles_end;
  AssignTasks(world.rank(), world.size(), FileCounts, nfiles_skip, nfiles_end);
  {
    int nhalo_local = 0;
    nhalo_local =
      accumulate(nhalo_per_groupfile.begin() + nfiles_skip, nhalo_per_groupfile.begin() + nfiles_end, nhalo_local);

    vector<HBTInt> HaloSizesLocal;
    HaloSizesLocal.resize(nhalo_local);

    for (int i = 0, ireader = 0; i < world.size(); i++, ireader++)
    {
      if (ireader == HBTConfig.MaxConcurrentIO)
      {
        ireader = 0;                     // reset reader count
        MPI_Barrier(world.Communicator); // wait for every thread to arrive.
      }
      if (i == world.rank()) // read
      {
        for (int iFile = nfiles_skip; iFile < nfiles_end; iFile++)
        {
          if (nhalo_per_groupfile[iFile]) // some files do not have groups
            ReadGroupLen(iFile, HaloSizesLocal.data() + offsethalo_per_groupfile[iFile] -
                                  offsethalo_per_groupfile[nfiles_skip]);
        }
      }
    }

    vector<int> buffersizes, bufferoffsets;
    if (world.rank() == root_node)
    {
      HaloSizesAll.resize(NhaloTotal);
      buffersizes.resize(world.size());
    }
    MPI_Gather(&nhalo_local, 1, MPI_INT, buffersizes.data(), 1, MPI_INT, root_node, world.Communicator);
    if (world.rank() == root_node)
      CompileOffsets(buffersizes, bufferoffsets);
    MPI_Gatherv(HaloSizesLocal.data(), HaloSizesLocal.size(), MPI_HBT_INT, HaloSizesAll.data(), buffersizes.data(),
                bufferoffsets.data(), MPI_HBT_INT, root_node, world.Communicator);

    /* Sanity check */
    HBTInt np_local = accumulate(HaloSizesLocal.begin(), HaloSizesLocal.end(), (HBTInt)0);
    HBTInt np_allproc = 0;
    MPI_Reduce(&np_local, &np_allproc, 1, MPI_HBT_INT, MPI_SUM, root_node, world.Communicator);
    if (world.rank() == root_node)
      assert(np_allproc == TotNumPartInGroups);
  }
}

/* Gathers in the root MPI rank how many particles each MPI rank has. */
void Gadget4Reader_t::CollectProcSizes(MpiWorker_t &world, const std::vector<Particle_t> &Particles)
{
  if (world.rank() == root_node)
    ProcLen.resize(world.size());

  HBTInt NumPartThisProc = Particles.size();
  MPI_Gather(&NumPartThisProc, 1, MPI_HBT_INT, ProcLen.data(), 1, MPI_HBT_INT, root_node, world.Communicator);
}

/* Gathers in the root MPI rank how many particles each MPI rank has. */
void Gadget4Reader_t::CollectProcSizes(MpiWorker_t &world, const ParticleSnapshot_t &partsnap)
{
  if (world.rank() == root_node)
    ProcLen.resize(world.size());

  HBTInt NumPartThisProc = partsnap.Particles.size();
  MPI_Gather(&NumPartThisProc, 1, MPI_HBT_INT, ProcLen.data(), 1, MPI_HBT_INT, root_node, world.Communicator);
}

/* Identifies how haloes are partitioned into segments across MPI ranks for in
 * GADGET4 group outputs. */
struct HaloPartitioner_t
{
  vector<HBTInt> ProcFirstHalo, ProcLastHalo;
  vector<CountBuffer_t> HaloSizesOnProc;

  void Fill(vector<HBTInt> HaloSizes, vector<HBTInt> ProcLen)
  {
    /* Number of ranks, what the first and last halo number is for each rank, and
     * corresponding halo particles. */
    int nproc = ProcLen.size();
    ProcFirstHalo.resize(nproc);
    ProcLastHalo.resize(nproc);

    /* Offset of halo particle numbers */
    vector<HBTInt> HaloOffsets;
    HBTInt NumPartInHalos = CompileOffsets(HaloSizes, HaloOffsets);
    HaloOffsets.push_back(NumPartInHalos);

    /* Offset of MPI rank particle numbers */
    vector<HBTInt> ProcOffsets;
    HBTInt NumPartInProcs = CompileOffsets(ProcLen, ProcOffsets);
    ProcOffsets.push_back(NumPartInProcs);

    assert(NumPartInHalos < NumPartInProcs);

    ProcFirstHalo[0] = 0; // First segment of first MPI rank is always halo 0.
    int iproc = 1;

    /* We iterate over MPI ranks to find what the HostId of its first and last
     * halo segment is. We make use of the fact that GADGET4 groups are sorted
     * in ascending group number. */
    if (nproc > 1)
    {
      for (HBTInt ihalo = 0; ihalo < HaloOffsets.size(); ihalo++)
      {
        while (ProcOffsets[iproc] <= HaloOffsets[ihalo])
        {
          /* Identify whether the first halo segment of an MPI rank is a
           * continuation of the last halo segment of the previous rank, or a
           * completely new halo. */
          if (HaloOffsets[ihalo] == ProcOffsets[iproc]) // New halo
          {
            ProcFirstHalo[iproc] = ihalo;
          }
          else // Continuation of the last halo segment of the previous rank.
          {
            ProcFirstHalo[iproc] = ihalo - 1;
          }

          /* The value of iproc will keep on increasing until we find the last
           * segment of a halo that is split across MPI ranks. */
          ProcLastHalo[iproc - 1] = ihalo - 1;
          iproc++;
        }
      }
    }

    /* Last segment of last MPI rank with at least one halo segment is always the
     * last halo. */
    ProcLastHalo[iproc - 1] = HaloSizes.size() - 1;

    /* The remainder of the processors have no haloes (i.e. HostId = -1) */
    for (; iproc < nproc; iproc++)
    {
      ProcFirstHalo[iproc] = -1;
      ProcLastHalo[iproc] = -1;
    }

    /* Get particle size of each local halo segment. */
    HaloSizesOnProc.resize(nproc);
#pragma omp parallel for
    for (iproc = 0; iproc < nproc; iproc++)
    {

      HBTInt first_halo = ProcFirstHalo[iproc];
      HBTInt last_halo = ProcLastHalo[iproc];
      HBTInt nhalos = last_halo - first_halo + 1;
      if (first_halo < 0 || last_halo < 0)
        nhalos = 0;
      auto &local_halo_sizes = HaloSizesOnProc[iproc];
      local_halo_sizes.resize(nhalos);

      if (nhalos > 0)
        local_halo_sizes[0] = min(HaloOffsets[first_halo + 1], ProcOffsets[iproc + 1]) - ProcOffsets[iproc];
      for (HBTInt i = 1; i < nhalos - 1; i++)
        local_halo_sizes[i] = HaloSizes[i + first_halo];
      if (nhalos > 1)
        local_halo_sizes.back() = min(HaloSizes[last_halo], ProcOffsets[iproc + 1] - HaloOffsets[last_halo]);
    }
  }
};

/* Loads particle properties in parallel */
void Gadget4Reader_t::LoadParticleProperties(MpiWorker_t &world, vector<Particle_t> &Particles)
{

  HBTInt nfiles_skip, nfiles_end;
  AssignTasks(world.rank(), world.size(), Header.NumberOfFiles, nfiles_skip, nfiles_end);
  {
    HBTInt np = 0;
    np = accumulate(np_file.begin() + nfiles_skip, np_file.begin() + nfiles_end, np);
    Particles.resize(np);
  }

  for (int i = 0, ireader = 0; i < world.size(); i++, ireader++)
  {
    if (ireader == HBTConfig.MaxConcurrentIO)
    {
      ireader = 0;                     // reset reader count
      MPI_Barrier(world.Communicator); // wait for every thread to arrive.
    }
    if (i == world.rank()) // read
    {
      for (int iFile = nfiles_skip; iFile < nfiles_end; iFile++)
      {
        ReadSnapshot(iFile, Particles.data() + offset_file[iFile] - offset_file[nfiles_skip]);
      }
    }
  }
}

/* Load particle host halo. */
void Gadget4Reader_t::LoadParticleHosts(MpiWorker_t &world, vector<Particle_t> &Particles)
{

  /* We obtain the total length of haloes and how many particles each MPI rank
   * contains in the root MPI rank. */
  LoadHaloSizes(world);
  CollectProcSizes(world, Particles);

  /* Identify how haloes are partitioned across MPI ranks. */
  HaloPartitioner_t HaloPartitioner;
  if (world.rank() == root_node)
  {
    HaloPartitioner.Fill(HaloSizesAll, ProcLen);
  }

  /* Tell each MPI rank what the minimum and maximum halo number they contain */
  HBTInt first_halo, last_halo;
  MPI_Scatter(HaloPartitioner.ProcFirstHalo.data(), 1, MPI_HBT_INT, &first_halo, 1, MPI_HBT_INT, root_node,
              world.Communicator);
  MPI_Scatter(HaloPartitioner.ProcLastHalo.data(), 1, MPI_HBT_INT, &last_halo, 1, MPI_HBT_INT, root_node,
              world.Communicator);

  if (world.rank() == root_node)
  {
    cout << "First halo for each rank is: " <<  HaloPartitioner.ProcFirstHalo <<  endl;
    cout << "Last halo for each rank is: " <<  HaloPartitioner.ProcLastHalo <<  endl;
  }

  /* The root rank will send the size of each halo segment contained in non-root
   * MPI ranks. */
  vector<HBTInt> local_halo_sizes;
  if (world.rank() == root_node)
  {
    for (int rank = 0; rank < world.size(); rank++)
    {
      if (rank == root_node)
        local_halo_sizes = HaloPartitioner.HaloSizesOnProc[rank];
      else
      {
        auto &sendarr = HaloPartitioner.HaloSizesOnProc[rank];
        MPI_Send(sendarr.data(), sendarr.size(), MPI_HBT_INT, rank, 0, world.Communicator);
      }
    }
  }
  else
  {
    MPI_Status stat;
    MPI_Probe(root_node, 0, world.Communicator, &stat);
    int nhalos;
    MPI_Get_count(&stat, MPI_HBT_INT, &nhalos);
    local_halo_sizes.resize(nhalos);
    MPI_Recv(local_halo_sizes.data(), nhalos, MPI_HBT_INT, root_node, 0, world.Communicator, MPI_STATUS_IGNORE);
  }

  /* Retrieve the particle offsets of each halo segment, measured relative to
   * the first particle in each MPI rank. */
  vector<HBTInt> local_halo_offsets;
  HBTInt np = CompileOffsets(local_halo_sizes, local_halo_offsets);
  local_halo_offsets.push_back(np);

#pragma omp parallel for default(shared)
  for (HBTInt i = 0; i < local_halo_sizes.size(); i++)
  {
    /* The range of particles that belong to the current halo segment. */
    auto start_particle = Particles.begin() + local_halo_offsets[i];
    auto end_particle = Particles.begin() + local_halo_offsets[i + 1];

    /* We assign a HostHalo to each particle based on its host. This needs to be done because
     * when we find host haloes we need up-to-date information for particles that are part of
     * subhaloes. */
    for (auto particle_it = start_particle; particle_it != end_particle; ++particle_it)
      particle_it->HostId = i + first_halo;
  }

  /* The rest of the particles have no host, so initialise them to -1 */
  for (auto particle_it = Particles.begin() + local_halo_offsets[local_halo_sizes.size()]; particle_it != Particles.end(); ++particle_it)
    particle_it->HostId = -1;

  /* Sanity check */
  HBTInt np_tot = 0;
  MPI_Reduce(&np, &np_tot, 1, MPI_HBT_INT, MPI_SUM, root_node, world.Communicator);
  if (world.rank() == root_node)
    assert(np_tot == TotNumPartInGroups);
}

/* Load group segments residing on the current proc */
HBTInt Gadget4Reader_t::LoadLocalGroups(MpiWorker_t &world, const vector<Particle_t> &Particles, vector<Halo_t> &Halos)
{
  /* Identify how haloes are partitioned across MPI ranks. */
  HaloPartitioner_t HaloPartitioner;
  if (world.rank() == root_node)
  {
    HaloPartitioner.Fill(HaloSizesAll, ProcLen);
  }

  /* Tell each MPI rank what the minimum and maximum halo number they contain */
  HBTInt first_halo, last_halo;
  MPI_Scatter(HaloPartitioner.ProcFirstHalo.data(), 1, MPI_HBT_INT, &first_halo, 1, MPI_HBT_INT, root_node,
              world.Communicator);
  MPI_Scatter(HaloPartitioner.ProcLastHalo.data(), 1, MPI_HBT_INT, &last_halo, 1, MPI_HBT_INT, root_node,
              world.Communicator);

  /* The root rank will send the size of each halo segment contained in non-root
   * MPI ranks. */
  vector<HBTInt> local_halo_sizes;
  if (world.rank() == root_node)
  {
    for (int rank = 0; rank < world.size(); rank++)
    {
      if (rank == root_node)
        local_halo_sizes = HaloPartitioner.HaloSizesOnProc[rank];
      else
      {
        auto &sendarr = HaloPartitioner.HaloSizesOnProc[rank];
        MPI_Send(sendarr.data(), sendarr.size(), MPI_HBT_INT, rank, 0, world.Communicator);
      }
    }
  }
  else
  {
    MPI_Status stat;
    MPI_Probe(root_node, 0, world.Communicator, &stat);
    int nhalos;
    MPI_Get_count(&stat, MPI_HBT_INT, &nhalos);
    local_halo_sizes.resize(nhalos);
    MPI_Recv(local_halo_sizes.data(), nhalos, MPI_HBT_INT, root_node, 0, world.Communicator, MPI_STATUS_IGNORE);
  }

  /* We prepare the halo segments to receive particles. */
  Halos.resize(local_halo_sizes.size());
  for (HBTInt i = 0; i < Halos.size(); i++)
  {
    Halos[i].HaloId = i + first_halo;
    Halos[i].Particles.resize(local_halo_sizes[i]);
  }

  if (Halos.size())
    assert(Halos.back().HaloId == last_halo);

  /* Assign local particles to halo */
  vector<HBTInt> local_halo_offsets;
  HBTInt np = CompileOffsets(local_halo_sizes, local_halo_offsets);
  local_halo_offsets.push_back(np);

#pragma omp parallel for default(shared)
  for (HBTInt i = 0; i < Halos.size(); i++)
  {
    /* The range of particles that belong to the current halo segment. */
    auto start_particle = Particles.begin() + local_halo_offsets[i];
    auto end_particle = Particles.begin() + local_halo_offsets[i + 1];

    /* We copy the particle information (e.g. coordinates, velocities, etc) to halo segment,
     * because this information will be used for central subhaloes. */
    Halos[i].Particles.assign(start_particle, end_particle);
  }

  /* Sanity check */
  HBTInt np_tot = 0;
  MPI_Reduce(&np, &np_tot, 1, MPI_HBT_INT, MPI_SUM, root_node, world.Communicator);
  if (world.rank() == root_node)
    assert(np_tot == TotNumPartInGroups);

  return np;
}

void Gadget4Reader_t::LoadGroups(MpiWorker_t &world, const ParticleSnapshot_t &partsnap, vector<Halo_t> &Halos)
{
  SetSnapshot(partsnap.GetSnapshotId());

  /* Particles have host information already, so we create halo segments. */
  struct HaloLen_t
  {
    HBTInt haloid;
    HBTInt np;
    HaloLen_t(){};
    HaloLen_t(HBTInt haloid, HBTInt np) : haloid(haloid), np(np)
    {
    }
  };
  vector<HaloLen_t> HaloLen;

  HBTInt curr_host_id = -1;
  for (auto &&p : partsnap.Particles)
  {
    if (p.HostId == -1)
      break; // NullGroupId comes last
    if (p.HostId != curr_host_id)
    {
      curr_host_id = p.HostId;
      HaloLen.emplace_back(curr_host_id, 1);
    }
    else
      HaloLen.back().np++;
  }
  Halos.resize(HaloLen.size());
  for (HBTInt i = 0; i < Halos.size(); i++)
  {
    Halos[i].HaloId = HaloLen[i].haloid;
    Halos[i].Particles.resize(HaloLen[i].np);
  }
  auto p_in = partsnap.Particles.begin();
  for (auto &&h : Halos)
  {
    for (auto &&p : h.Particles)
    {
      p = *p_in;
      ++p_in;
    }
  }

  ExchangeAndMerge(world, Halos);
  global_timer.Tick("halo_comms", world.Communicator);
  HBTConfig.GroupLoadedFullParticle = true;
}

bool IsGadget4Group(const string &GroupFileFormat)
{
  return GroupFileFormat.substr(0, 11) == "gadget4_hdf";
}
} // namespace Gadget4Reader
