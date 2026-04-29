using namespace std;
#include <iostream>
#include <numeric>
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

inline bool CompParticleHost(const Particle_t &a, const Particle_t &b)
{
  return a.HostId < b.HostId;
}

namespace Gadget4Reader
{

typedef vector<HBTInt> CountBuffer_t;

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

/* Gets the number of particles of each type in the current file. */
void Gadget4Reader_t::GetParticleCountInFile(hid_t file, int np[])
{
  ReadAttribute(file, "Header", "NumPart_ThisFile", H5T_NATIVE_INT, np);
#ifdef DM_ONLY
  for (int i = 0; i < TypeMax; i++)
    if (i != TypeDM)
      np[i] = 0;
#endif // DM_ONLY
}

HBTInt Gadget4Reader_t::CompileFileOffsets(int nfiles)
{
  HBTInt offset = 0;
  NumberParticlesPerFile.reserve(nfiles);
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

    NumberParticlesPerFile.push_back(np);
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
  /* Open the simulation output to read */
  string filename;
  GetFileName(ifile, filename);
  hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

  /* Get the number of particles per type in this file. */
  std::vector<int> NumberParticlesPerType(TypeMax);
  GetParticleCountInFile(file, NumberParticlesPerType.data());

  /* Create an offset vector, used to position particles */
  std::vector<HBTInt> offset_this(TypeMax);
  CompileOffsets(NumberParticlesPerType, offset_this);

  /* Load each particle type. */
  for (int itype = 0; itype < TypeMax; itype++)
  {
    int NumberParticlesThisType = NumberParticlesPerType[itype];

    /* Nothing to read. */
    if (NumberParticlesThisType == 0)
      continue;

    /* Location of the Particle array the upcoming particles will be placed. */
    auto ParticlesThisType = ParticlesInFile + offset_this[itype];

    stringstream grpname;
    grpname << "PartType" << itype;

    if (!H5Lexists(file, grpname.str().c_str(), H5P_DEFAULT))
      continue;

    hid_t particle_data = H5Gopen2(file, grpname.str().c_str(), H5P_DEFAULT);
    check_id_size(particle_data);

    /* Particle cordinates */
    {
      std::vector<HBTxyz> x(NumberParticlesThisType);
      ReadDataset(particle_data, "Coordinates", H5T_HBTReal, x.data());
      if (HBTConfig.PeriodicBoundaryOn)
      {
#pragma omp parallel for
        for (int i = 0; i < NumberParticlesThisType; i++)
          for (int j = 0; j < 3; j++)
            x[i][j] = position_modulus(x[i][j], Header.BoxSize);
      }
#pragma omp parallel for
      for (int i = 0; i < NumberParticlesThisType; i++)
        copyHBTxyz(ParticlesThisType[i].ComovingPosition, x[i]);
    }

    /* Particle velocities */
    {
      /* Get the a-scale factor dependence on velocity. */
      HBTReal aexp;
      ReadAttribute(particle_data, "Velocities", "a_scaling", H5T_HBTReal, &aexp);

      vector<HBTxyz> v(NumberParticlesThisType);
      ReadDataset(particle_data, "Velocities", H5T_HBTReal, v.data());

#pragma omp parallel for
      for (int i = 0; i < NumberParticlesThisType; i++)
        for (int j = 0; j < 3; j++)
          ParticlesThisType[i].PhysicalVelocity[j] = v[i][j] * pow(Header.ScaleFactor, aexp);
    }

    /* Particle IDs */
    {
      vector<HBTInt> id(NumberParticlesThisType);
      ReadDataset(particle_data, "ParticleIDs", H5T_HBTInt, id.data());
#pragma omp parallel for
      for (int i = 0; i < NumberParticlesThisType; i++)
        ParticlesThisType[i].Id = id[i];
    }

    /* Particle masses */
    if (Header.mass[itype] == 0)
    {
      vector<HBTReal> m(NumberParticlesThisType);
      ReadDataset(particle_data, "Masses", H5T_HBTReal, m.data());
#pragma omp parallel for
      for (int i = 0; i < NumberParticlesThisType; i++)
        ParticlesThisType[i].Mass = m[i];
    }
    else
    {
#pragma omp parallel for
      for (int i = 0; i < NumberParticlesThisType; i++)
        ParticlesThisType[i].Mass = Header.mass[itype];
    }

    /* Particle internal energies */
#ifndef DM_ONLY
#ifdef HAS_THERMAL_ENERGY
    if (itype == 0)
    {
      vector<HBTReal> u(NumberParticlesThisType);
      ReadDataset(particle_data, "InternalEnergy", H5T_HBTReal, u.data());
#pragma omp parallel for
      for (int i = 0; i < NumberParticlesThisType; i++)
        ParticlesThisType[i].InternalEnergy = u[i];
    }
    else // Zero out internal energy for non-gas particles
    {
#pragma omp parallel for
      for (int i = 0; i < NumberParticlesThisType; i++)
        ParticlesThisType[i].InternalEnergy = 0.0;
    }
#endif // HAS_THERMAL_ENERGY

    /* Particle types */
    {
      ParticleType_t t = static_cast<ParticleType_t>(itype);
#pragma omp parallel for
      for (int i = 0; i < NumberParticlesThisType; i++)
        ParticlesThisType[i].Type = t;
    }
#endif // DM_ONLY

    H5Gclose(particle_data);
  }
  H5Fclose(file);
}

/* Reads the total number of all particle types for the haloes in the current
 * group file. */
void Gadget4Reader_t::ReadGroupLen(int ifile, HBTInt *buf)
{
  /* Open the group file. */
  string filename;
  GetGroupFileName(ifile, filename);
  hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

  /* How many groups are there in this file. */
  HBTInt NumGroups;
  ReadAttribute(file, "Header", "Ngroups_ThisFile", H5T_HBTInt, &NumGroups);

  /* Read the dataset. */
  if (NumGroups > 0)
    ReadDataset(file, "Group/GroupLen", H5T_HBTInt, buf);

  H5Fclose(file);
}

/* Reads the number of different particle types for the haloes in the current
 * group file. */
void Gadget4Reader_t::ReadGroupLenPerType(int ifile, std::array<HBTInt, TypeMax> *buf)
{
  /* Open the group file. */
  string filename;
  GetGroupFileName(ifile, filename);
  hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

  /* How many groups are there in this file. */
  HBTInt NumGroups;
  ReadAttribute(file, "Header", "Ngroups_ThisFile", H5T_HBTInt, &NumGroups);

  /* Read the dataset. */
  if (NumGroups > 0)
    ReadDataset(file, "Group/GroupLenType", H5T_HBTInt, buf);

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
  world.SyncContainer(NumberParticlesPerFile, MPI_HBT_INT, root_node);
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

  /* This will be used to determine which particles are hostless when
   * constraining subhaloes to their assigned hosts. The value of NullGroupId
   * is fixed to be the maximum positive representation of HBTInt. */
  HBTConfig.ParticleNullGroupId = NullGroupId;

  /* Read physical properties of particles and their particle type. */
  LoadParticleProperties(world, Particles);

  /* Read host halo of each particle. GADGET4 has particle and group information
   * separate, hence the need for this separate function call. */
  LoadParticleHosts(world, Particles);
}

/* Loads halo particle sizes from all files in parallel and gathers to root MPI
 * rank. */
void Gadget4Reader_t::LoadHaloSizes(MpiWorker_t &world)
{

  /* File metadata */
  int NumberGroupFiles;
  HBTInt TotalNumberGroups;
  std::vector<HBTInt> NumberGroupsPerFile;
  std::vector<HBTInt> OffsetGroupsPerFile;

  if (world.rank() == root_node)
  {
    NumberGroupFiles = ReadGroupFileCounts(0);
    NumberGroupsPerFile.resize(NumberGroupFiles);
    OffsetGroupsPerFile.resize(NumberGroupFiles);
    TotalNumberGroups = CompileGroupFileOffsets(NumberGroupsPerFile, OffsetGroupsPerFile);
    TotalNumberGroupParticles = ReadGroupFileTotParticles(0);
  }

  world.SyncAtom(NumberGroupFiles, MPI_INT, root_node);
  world.SyncAtom(TotalNumberGroups, MPI_HBT_INT, root_node);
  world.SyncContainer(NumberGroupsPerFile, MPI_HBT_INT, root_node);
  world.SyncContainer(OffsetGroupsPerFile, MPI_HBT_INT, root_node);
  world.SyncAtom(TotalNumberGroupParticles, MPI_HBT_INT, root_node);

  /* Assign each task to a range of files to read in parallel. */
  HBTInt FirstFileIndex, LastFileIndex;
  AssignTasks(world.rank(), world.size(), NumberGroupFiles, FirstFileIndex, LastFileIndex);

  /* Allocate sufficient memory to hold halo lengths in this rank. */
  HBTInt NumberGroupsInRank = std::accumulate(NumberGroupsPerFile.begin() + FirstFileIndex, NumberGroupsPerFile.begin() + LastFileIndex, 0);
  std::vector<HBTInt> LocalHaloSizes(NumberGroupsInRank);
  std::vector<std::array<HBTInt, TypeMax>> LocalHaloSizesPerType(NumberGroupsInRank);

  for (int i = 0, ireader = 0; i < world.size(); i++, ireader++)
  {
    if (ireader == HBTConfig.MaxConcurrentIO)
    {
      ireader = 0;                     // reset reader count
      MPI_Barrier(world.Communicator); // wait for every thread to arrive.
    }
    if (world.rank() == i) // read
    {
      for (int FileIndex = FirstFileIndex; FileIndex < LastFileIndex; FileIndex++)
      {
        if (NumberGroupsPerFile[FileIndex]) // some files do not have groups
        {
          ReadGroupLen(FileIndex, LocalHaloSizes.data() + OffsetGroupsPerFile[FileIndex] -
                                OffsetGroupsPerFile[FirstFileIndex]);
          ReadGroupLenPerType(FileIndex, LocalHaloSizesPerType.data() + OffsetGroupsPerFile[FileIndex] -
                                OffsetGroupsPerFile[FirstFileIndex]);
        }
      }
    }
  }

  /* Sanity check: the sum of each particle type part of the halo should be equal to
   * the total number. */
#ifndef NDEBUG
  for (HBTInt halo_i = 0; halo_i < NumberGroupsInRank; halo_i++)
  {
    HBTInt SumParticleTypes = 0;
    for (int PartType = 0; PartType < TypeMax; PartType++)
      SumParticleTypes += LocalHaloSizesPerType[halo_i][PartType];

    assert(LocalHaloSizes[halo_i] == SumParticleTypes);
  }
#endif

  /* In this block we will gather the halo sizes loaded by each rank into the
   * root rank. The information will eventually be used to determine particle
   * offsets. NOTE: this is currently redudant, but it can be useful for
   * debugging */
  {
    std::vector<int> NumberHalosPerRank, OffsetHaloesPerRank;

    if (world.rank() == root_node)
    {
      AllHaloSizes.resize(TotalNumberGroups);
      NumberHalosPerRank.resize(world.size());
    }
    MPI_Gather(&NumberGroupsInRank, 1, MPI_INT, NumberHalosPerRank.data(), 1, MPI_INT, root_node, world.Communicator);

    /* Create an offset value for each halo value. */
    if (world.rank() == root_node)
      CompileOffsets(NumberHalosPerRank, OffsetHaloesPerRank);

    MPI_Gatherv(LocalHaloSizes.data(), LocalHaloSizes.size(), MPI_HBT_INT, AllHaloSizes.data(), NumberHalosPerRank.data(),
                OffsetHaloesPerRank.data(), MPI_HBT_INT, root_node, world.Communicator);
  }

  /* In this block we will gather the halo sizes (per particle type) loaded by
   * each rank into the root rank. The information will eventually be used to
   * determine particle offsets per type. */
  {
    std::vector<int> NumberHalosPerRank, OffsetHaloesPerRank;
    if (world.rank() == root_node)
    {
      AllHaloSizesPerType.resize(TotalNumberGroups);
      NumberHalosPerRank.resize(world.size());
    }
    MPI_Gather(&NumberGroupsInRank, 1, MPI_INT, NumberHalosPerRank.data(), 1, MPI_INT, root_node, world.Communicator);

    /* Create an offset value for each halo value. */
    if (world.rank() == root_node)
      CompileOffsets(NumberHalosPerRank, OffsetHaloesPerRank);

    /* Gather in the root rank. We create a custom MPI dtype because each vector
     * element is an array of size MaxType. */
    MPI_Datatype MPI_HBT_ARRAY;
    MPI_Type_contiguous(TypeMax, MPI_HBT_INT, &MPI_HBT_ARRAY);
    MPI_Type_commit(&MPI_HBT_ARRAY);
    MPI_Gatherv(LocalHaloSizesPerType.data(), LocalHaloSizesPerType.size(), MPI_HBT_ARRAY, AllHaloSizesPerType.data(), NumberHalosPerRank.data(),
                OffsetHaloesPerRank.data(), MPI_HBT_ARRAY, root_node, world.Communicator);
    MPI_Type_free(&MPI_HBT_ARRAY);
  }

  /* Sanity check: the sum of each particle type part of the halo should be equal to
   * the total number (after communicating). */
#ifndef NDEBUG
  if (world.rank() == 0)
  {
    for (HBTInt halo_i = 0; halo_i < TotalNumberGroups; halo_i++)
    {
      HBTInt SumParticleTypes = 0;
      for (int PartType = 0; PartType < TypeMax; PartType++)
        SumParticleTypes += AllHaloSizesPerType[halo_i][PartType];

      assert(AllHaloSizes[halo_i] == SumParticleTypes);
    }
  }
#endif

  /* Sanity check: the sum of all halo sizes should equal the total number of particles
   * in haloes. */
  {
    HBTInt np_local = std::accumulate(LocalHaloSizes.begin(), LocalHaloSizes.end(), 0);
    HBTInt np_allproc = 0;

    MPI_Reduce(&np_local, &np_allproc, 1, MPI_HBT_INT, MPI_SUM, root_node, world.Communicator);
    if (world.rank() == root_node)
      assert(np_allproc == TotalNumberGroupParticles);
  }
}

/* Gathers in the root MPI rank how many particles each MPI rank has. */
void Gadget4Reader_t::GetNumberParticlesPerRank(MpiWorker_t &world, const std::vector<Particle_t> &Particles)
{
  if (world.rank() == root_node)
    NumberParticlesPerRank.resize(world.size());

  HBTInt NumPartThisProc = Particles.size();
  MPI_Gather(&NumPartThisProc, 1, MPI_HBT_INT, NumberParticlesPerRank.data(), 1, MPI_HBT_INT, root_node, world.Communicator);
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

  /* Allocate a range of files for each task to read. */
  HBTInt FirstFileIndex, LastFileIndex;
  AssignTasks(world.rank(), world.size(), Header.NumberOfFiles, FirstFileIndex, LastFileIndex);

  /* Allocate sufficient space in each task to hold the incoming particle information. */
  {
    HBTInt NumberParticlesInRank = std::accumulate(NumberParticlesPerFile.begin() + FirstFileIndex, NumberParticlesPerFile.begin() + LastFileIndex, 0);
    Particles.resize(NumberParticlesInRank);
  }

  for (int i = 0, ireader = 0; i < world.size(); i++, ireader++)
  {
    if (ireader == HBTConfig.MaxConcurrentIO)
    {
      ireader = 0;                     // reset reader count
      MPI_Barrier(world.Communicator); // wait for every thread to arrive.
    }
    if (world.rank() == i) // read
    {
      for (int FileIndex = FirstFileIndex; FileIndex < LastFileIndex; FileIndex++)
      {
        ReadSnapshot(FileIndex, Particles.data() + offset_file[FileIndex] - offset_file[FirstFileIndex]);
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
  GetNumberParticlesPerRank(world, Particles);

  /* Identify how haloes are partitioned across MPI ranks. */
  HaloPartitioner_t HaloPartitioner;
  if (world.rank() == root_node)
  {
    HaloPartitioner.Fill(AllHaloSizes, NumberParticlesPerRank);
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

  /* The rest of the particles have no host, so initialise them to NullGroupId */
  for (auto particle_it = Particles.begin() + local_halo_offsets[local_halo_sizes.size()]; particle_it != Particles.end(); ++particle_it)
    particle_it->HostId = NullGroupId;

  /* Sanity check */
  HBTInt np_tot = 0;
  MPI_Reduce(&np, &np_tot, 1, MPI_HBT_INT, MPI_SUM, root_node, world.Communicator);
  if (world.rank() == root_node)
    assert(np_tot == TotalNumberGroupParticles);
}

/* Group up particles in the same rank that share the same HostId, and create
 * a new halo segment within for Halos. */
void Gadget4Reader_t::CreateHaloSegments(const std::vector<Particle_t> &SnapshotParticles, std::vector<Halo_t> &Halos)
{
  /* Create a copy of the Particle snapshot, from which we will create the halo
   * segments. */
  std::vector<Particle_t> GroupParticles(SnapshotParticles);

  /* Sort particles by their FoF membership. */
  std::sort(GroupParticles.begin(), GroupParticles.end(), CompParticleHost);
  if (!GroupParticles.empty())
  {
    assert(GroupParticles.back().HostId <= NullGroupId);
    assert(GroupParticles.front().HostId >= 0);
  }

  struct HaloLen_t
  {
    HBTInt haloid;
    HBTInt np;
    HaloLen_t(){};
    HaloLen_t(HBTInt haloid, HBTInt np) : haloid(haloid), np(np)
    {
    }
  };
  std::vector<HaloLen_t> HaloLen;

  HBTInt curr_host_id = -1;
  for (auto &&p : GroupParticles)
  {
    /* NullGroupId comes last, so we can break the loop because there are no
     * particles in FoF groups. */
    if (p.HostId == NullGroupId)
      break;

    /* We have reached a new FoF group. Create a new segment for it. */
    if (p.HostId != curr_host_id)
    {
      curr_host_id = p.HostId;
      HaloLen.emplace_back(curr_host_id, 1);
    }
    /* We are still in the same FoF group segment, grow its size. */
    else
      HaloLen.back().np++;
  }

  /* Allocate space for each halo segment we have, including to store its
   * associated particles. */
  Halos.resize(HaloLen.size());
  for (HBTInt i = 0; i < Halos.size(); i++)
  {
    Halos[i].HaloId = HaloLen[i].haloid;
    Halos[i].Particles.resize(HaloLen[i].np);
  }

  /* Copy the particles into the same host. */
  auto p_in = GroupParticles.begin();
  for (auto &&h : Halos)
  {
    for (auto &&p : h.Particles)
    {
      p = *p_in;
      ++p_in;
    }
  }
}

/* Creates and populates the Halos of the current snapshot, which for GADGET4 is
 * done using the snapshot particles. */
void Gadget4Reader_t::LoadGroups(MpiWorker_t &world, const ParticleSnapshot_t &partsnap, std::vector<Halo_t> &Halos)
{
  SetSnapshot(partsnap.GetSnapshotId());

  /* We group up particles in the same rank that share the same HostId into
   * individual halo segments. */
  CreateHaloSegments(partsnap.Particles, Halos);

  /* NOTE: by communicating GADGET4 particles before creating halo segments, I
   * have indirectly created many more halo segments than what would be required
   * if I had not communicated them. This is because particles in GADGET4
   * outputs are grouped by FoF, but communicating particles groups them by a
   * hash function intended to balance particle load. */

  /* We have halo segments of the same FoF group split across MPI ranks. We need
   * to gather each unique FoF group in the same rank and merge each segment. */
  ExchangeAndMerge(world, Halos);

  global_timer.Tick("halo_comms", world.Communicator);
  HBTConfig.GroupLoadedFullParticle = true;
}

bool IsGadget4Group(const string &GroupFileFormat)
{
  return GroupFileFormat.substr(0, 7) == "gadget4";
}
} // namespace Gadget4Reader
