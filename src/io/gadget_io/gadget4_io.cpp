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

/* IMPORTANT TODO: handle the fact that INTs and HBT_INTs are mixed up in some
 * MPI communications and array/vector declarations. */

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

void Gadget4Reader_t::CompileSnapshotParticleOffsets(int NumberFiles)
{
  /* First we load how many particles there are per file. */
  NumberParticlesPerFile.resize(NumberFiles);
  NumberParticlesPerTypePerFile.resize(NumberFiles);

  for (int ifile = 0; ifile < NumberFiles; ifile++)
  {
    std::array<int, TypeMax> NumberParticlesPerTypeThisFile{};

    string filename;
    GetFileName(ifile, filename);
    hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    GetParticleCountInFile(file, NumberParticlesPerTypeThisFile.data());
    H5Fclose(file);

    /* Store number of particles per type in each file. */
    NumberParticlesPerTypePerFile[ifile] = NumberParticlesPerTypeThisFile;

    /* Store total number of particles in each file. */
    NumberParticlesPerFile[ifile] = std::accumulate(NumberParticlesPerTypeThisFile.begin(), NumberParticlesPerTypeThisFile.end(), 0);
  }

  /* Now store the offsets. */
  OffsetParticlesPerFile.resize(NumberFiles);
  OffsetParticlesPerTypePerFile.resize(NumberFiles);

  HBTInt ParticleOffset = 0;
  std::array<int, TypeMax> ParticleOffsetPerType{};
  for (int ifile = 0; ifile < NumberFiles; ifile++)
  {
    OffsetParticlesPerFile[ifile] = ParticleOffset;
    OffsetParticlesPerTypePerFile[ifile] = ParticleOffsetPerType;

    ParticleOffset += NumberParticlesPerFile[ifile];
    for(int PartType = 0; PartType < TypeMax; PartType++)
      ParticleOffsetPerType[PartType] += NumberParticlesPerTypePerFile[ifile][PartType];
  }
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
void Gadget4Reader_t::ReadSnapshot(int FirstFileIndex, int CurrentFileIndex, Particle_t *ParticlesInFile)
{
  /* Open the simulation output to read */
  string filename;
  GetFileName(CurrentFileIndex, filename);
  hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

  /* Load each particle type. */
  for (int PartType = 0; PartType < TypeMax; PartType++)
  {
    int NumberParticlesThisTypeToRead = NumberParticlesPerTypePerFile[CurrentFileIndex][PartType];

    /* Nothing to read. */
    if (NumberParticlesThisTypeToRead == 0)
      continue;

    /* Calculate where the current particle type starts in the Particle vector. */
    HBTInt OffsetPreviousParticleTypes = 0;
    for(int DummyPartType = 0; DummyPartType < PartType;  DummyPartType++)
      OffsetPreviousParticleTypes += NumberParticlesPerTypeThisRank[DummyPartType];

    /* Locate the start of the section where the upcoming particles will be placed. */
    auto ParticlesThisTypeToRead = ParticlesInFile             \
                                 + OffsetPreviousParticleTypes \
                                 + OffsetParticlesPerTypePerFile[CurrentFileIndex][PartType] - OffsetParticlesPerTypePerFile[FirstFileIndex][PartType];

    stringstream grpname;
    grpname << "PartType" << PartType;

    if (!H5Lexists(file, grpname.str().c_str(), H5P_DEFAULT))
      continue;

    hid_t particle_data = H5Gopen2(file, grpname.str().c_str(), H5P_DEFAULT);
    check_id_size(particle_data);

    /* Particle cordinates */
    {
      std::vector<HBTxyz> x(NumberParticlesThisTypeToRead);
      ReadDataset(particle_data, "Coordinates", H5T_HBTReal, x.data());
      if (HBTConfig.PeriodicBoundaryOn)
      {
#pragma omp parallel for
        for (int i = 0; i < NumberParticlesThisTypeToRead; i++)
          for (int j = 0; j < 3; j++)
            x[i][j] = position_modulus(x[i][j], Header.BoxSize);
      }
#pragma omp parallel for
      for (int i = 0; i < NumberParticlesThisTypeToRead; i++)
        copyHBTxyz(ParticlesThisTypeToRead[i].ComovingPosition, x[i]);
    }

    /* Particle velocities */
    {
      /* Get the a-scale factor dependence on velocity. */
      HBTReal aexp;
      ReadAttribute(particle_data, "Velocities", "a_scaling", H5T_HBTReal, &aexp);

      vector<HBTxyz> v(NumberParticlesThisTypeToRead);
      ReadDataset(particle_data, "Velocities", H5T_HBTReal, v.data());

#pragma omp parallel for
      for (int i = 0; i < NumberParticlesThisTypeToRead; i++)
        for (int j = 0; j < 3; j++)
          ParticlesThisTypeToRead[i].PhysicalVelocity[j] = v[i][j] * pow(Header.ScaleFactor, aexp);
    }

    /* Particle IDs */
    {
      vector<HBTInt> id(NumberParticlesThisTypeToRead);
      ReadDataset(particle_data, "ParticleIDs", H5T_HBTInt, id.data());
#pragma omp parallel for
      for (int i = 0; i < NumberParticlesThisTypeToRead; i++)
        ParticlesThisTypeToRead[i].Id = id[i];
    }

    /* Particle masses */
    if (Header.mass[PartType] == 0)
    {
      vector<HBTReal> m(NumberParticlesThisTypeToRead);
      ReadDataset(particle_data, "Masses", H5T_HBTReal, m.data());
#pragma omp parallel for
      for (int i = 0; i < NumberParticlesThisTypeToRead; i++)
        ParticlesThisTypeToRead[i].Mass = m[i];
    }
    else
    {
#pragma omp parallel for
      for (int i = 0; i < NumberParticlesThisTypeToRead; i++)
        ParticlesThisTypeToRead[i].Mass = Header.mass[PartType];
    }

    /* Particle internal energies */
#ifndef DM_ONLY
#ifdef HAS_THERMAL_ENERGY
    if (PartType == 0)
    {
      vector<HBTReal> u(NumberParticlesThisTypeToRead);
      ReadDataset(particle_data, "InternalEnergy", H5T_HBTReal, u.data());
#pragma omp parallel for
      for (int i = 0; i < NumberParticlesThisTypeToRead; i++)
        ParticlesThisTypeToRead[i].InternalEnergy = u[i];
    }
    else // Zero out internal energy for non-gas particles
    {
#pragma omp parallel for
      for (int i = 0; i < NumberParticlesThisTypeToRead; i++)
        ParticlesThisTypeToRead[i].InternalEnergy = 0.0;
    }
#endif // HAS_THERMAL_ENERGY

    /* Particle types */
    {
      ParticleType_t t = static_cast<ParticleType_t>(PartType);
#pragma omp parallel for
      for (int i = 0; i < NumberParticlesThisTypeToRead; i++)
        ParticlesThisTypeToRead[i].Type = t;
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

  /* The root MPI rank reads snapshot header and how many particles of each type
   * we have in each file.*/
  if (world.rank() == root_node)
  {
    ReadHeader(0, Header);
    CompileSnapshotParticleOffsets(Header.NumberOfFiles);
  }

  /* Tell other MPI ranks about the information the root MPI rank just loaded. */
  MPI_Bcast(&Header, 1, MPI_Gadget4Header_t, root_node, world.Communicator);
  world.SyncContainer(NumberParticlesPerFile, MPI_HBT_INT, root_node);
  world.SyncContainer(OffsetParticlesPerFile, MPI_HBT_INT, root_node);

  /* Need a custom MPI type here. TODO: Make this MPI dtype creation into its
   * own function. */
  MPI_Datatype MPI_HBT_ARRAY;
  MPI_Type_contiguous(TypeMax, MPI_INT, &MPI_HBT_ARRAY);
  MPI_Type_commit(&MPI_HBT_ARRAY);
  world.SyncContainer(NumberParticlesPerTypePerFile, MPI_HBT_ARRAY, root_node);
  world.SyncContainer(OffsetParticlesPerTypePerFile, MPI_HBT_ARRAY, root_node);
  MPI_Type_free(&MPI_HBT_ARRAY);

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
     * element is an array of size TypeMax. */
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

/* Gathers in the root MPI rank how many particles each MPI rank has. Also, how
 * many of each type it has loaded. */
void Gadget4Reader_t::GetNumberParticlesPerRank(MpiWorker_t &world, const std::vector<Particle_t> &Particles)
{
  if (world.rank() == root_node)
  {
    NumberParticlesPerRank.resize(world.size());
    NumberParticlesPerTypePerRank.resize(world.size());
  }

  /* First the total number of particles per rank. */
  HBTInt NumberParticlesThisRank = Particles.size();
  MPI_Gather(&NumberParticlesThisRank, 1, MPI_HBT_INT, NumberParticlesPerRank.data(), 1, MPI_HBT_INT, root_node, world.Communicator);

  /* Now the number of particles per type.*/
  MPI_Datatype MPI_HBT_ARRAY;
  MPI_Type_contiguous(TypeMax, MPI_HBT_INT, &MPI_HBT_ARRAY);
  MPI_Type_commit(&MPI_HBT_ARRAY);
  MPI_Gather(NumberParticlesPerTypeThisRank.data(), 1, MPI_HBT_ARRAY, NumberParticlesPerTypePerRank.data(), 1, MPI_HBT_ARRAY, root_node, world.Communicator);
  MPI_Type_free(&MPI_HBT_ARRAY);
}

/* Identifies how haloes are partitioned into segments across MPI ranks for in
 * GADGET4 group outputs. */
struct HaloPartitioner_t
{

  typedef std::vector<HBTInt> CountBuffer_t;
  std::vector<CountBuffer_t> HaloSizesOnProc;
  std::vector<HBTInt> RankFirstHalo, RankLastHalo;

  void Fill(vector<HBTInt> HaloSizes, vector<HBTInt> ProcLen)
  {
    /* Number of ranks, what the first and last halo number is for each rank, and
     * corresponding halo particles. */
    int NumberRanks = ProcLen.size();
    RankFirstHalo.resize(NumberRanks);
    RankLastHalo.resize(NumberRanks);

    /* Offset of halo particle numbers */
    vector<HBTInt> HaloOffsets;
    HBTInt NumPartInHalos = CompileOffsets(HaloSizes, HaloOffsets);
    HaloOffsets.push_back(NumPartInHalos);

    /* Offset of MPI rank particle numbers */
    vector<HBTInt> ProcOffsets;
    HBTInt NumPartInProcs = CompileOffsets(ProcLen, ProcOffsets);
    ProcOffsets.push_back(NumPartInProcs);

    /* Sanity check: there cannot be more particles in haloes than total
     * particles. We can have equal numbers, e.g. if there are no particles of
     * a given type. */
    assert(NumPartInHalos <= NumPartInProcs);

    RankFirstHalo[0] = 0; // First segment of first MPI rank is always halo 0.
    int iproc = 1;

    /* We iterate over MPI ranks to find what the HostId of its first and last
     * halo segment is. We make use of the fact that GADGET4 groups are sorted
     * in ascending group number. */
    if (NumberRanks > 1)
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
            RankFirstHalo[iproc] = ihalo;
          }
          else // Continuation of the last halo segment of the previous rank.
          {
            RankFirstHalo[iproc] = ihalo - 1;
          }

          /* The value of iproc will keep on increasing until we find the last
           * segment of a halo that is split across MPI ranks. */
          RankLastHalo[iproc - 1] = ihalo - 1;
          iproc++;
        }
      }
    }

    /* Last segment of last MPI rank with at least one halo segment is always the
     * last halo. */
    RankLastHalo[iproc - 1] = HaloSizes.size() - 1;

    /* The remainder of the processors have no haloes (i.e. HostId = -1) */
    for (; iproc < NumberRanks; iproc++)
    {
      RankFirstHalo[iproc] = -1;
      RankLastHalo[iproc] = -1;
    }

    /* Get particle size of each local halo segment. */
    HaloSizesOnProc.resize(NumberRanks);
#pragma omp parallel for
    for (iproc = 0; iproc < NumberRanks; iproc++)
    {

      HBTInt first_halo = RankFirstHalo[iproc];
      HBTInt last_halo = RankLastHalo[iproc];
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

/* Loads particle properties in parallel and places particles of the same PartType
 * in a contiguous section of Particles */
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

  /* Fill particle vector with dummy Type values, which will be used for a sanity
   * check before this function completes. */
#ifndef NDEBUG
  for (HBTInt part_i = 0; part_i < Particles.size(); part_i++)
    Particles[part_i] = -1;
#endif

  /* Compute how many particles of each type this rank will read. Used to
   * position particles at the correct location in the particle vector.*/
  {
    NumberParticlesPerTypeThisRank.resize(TypeMax);
    for (int FileIndex = FirstFileIndex; FileIndex < LastFileIndex; FileIndex++)
      for (int PartType = 0; PartType < TypeMax; PartType++)
        NumberParticlesPerTypeThisRank[PartType] += NumberParticlesPerTypePerFile[FileIndex][PartType];
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
        /* We pass the whole vector particle vector. We handle offsetting within
         * this function call. */
        ReadSnapshot(FirstFileIndex, FileIndex, Particles.data());
      }
    }
  }

  /* Sanity check: Every particle entry should have been filled in, and particle
   * types should be together */
#ifndef NDEBUG
  for (int PartType = 0; PartType < TypeMax; PartType++)
  {
    /* Calculate where the current particle type starts and ends in the Particle vector. */
    HBTInt FirstIndex = 0;
    for(int particle = 0; particle < PartType;  particle++)
      FirstIndex += NumberParticlesPerTypeThisRank[particle];

    HBTInt LastIndex= std::accumulate(NumberParticlesPerTypeThisRank.begin(), NumberParticlesPerTypeThisRank.begin() + PartType + 1, 0);

    /* Check whether the type is as expected, which means that they have at
     * least been loaded in contiguous chunks. */
    for(HBTInt part_i = FirstIndex; part_i < LastIndex; part_i++)
      assert(Particles[part_i].Type == PartType);
  }
#endif
}

/* Load particle host halo. */
void Gadget4Reader_t::LoadParticleHosts(MpiWorker_t &world, vector<Particle_t> &Particles)
{
  /* For a sanity check at the end of this function. */
  HBTInt LocalNumberParticlesAssigned = 0;

  /* We default-initialise particles to be in no hosts. */
#pragma omp parallel for
  for (size_t part_i = 0; part_i < Particles.size(); part_i++)
    Particles[part_i].HostId = NullGroupId;

  /* We obtain the total length of haloes and how many particles each MPI rank
   * contains in the root MPI rank. */
  LoadHaloSizes(world);
  GetNumberParticlesPerRank(world, Particles);

  /* We iterate over particle types because we will use GroupLenType and the fact
   * that the same particle types are contiguous in Particles vector to assign
   * them the correct FoF group membership. */
  for(int PartType = 0; PartType < TypeMax; PartType++)
  {
    /* Identify how haloes are partitioned across MPI ranks. */
    HaloPartitioner_t HaloPartitioner;
    if (world.rank() == root_node)
    {
      /* We create two vectors for the current particle type, one holding how many
       * particles per FoF group and another per rank. It will be used to determine
       * which range of particles get assigned which FoF group. */
      std::vector<HBTInt> AllHaloSizesThisType(AllHaloSizesPerType.size());
      for(size_t halo_i = 0; halo_i < AllHaloSizesThisType.size(); halo_i++)
        AllHaloSizesThisType[halo_i] = AllHaloSizesPerType[halo_i][PartType];

      std::vector<HBTInt> NumberParticlesThisTypePerRank(NumberParticlesPerTypePerRank.size());
      for(int rank = 0; rank < world.size(); rank++)
        NumberParticlesThisTypePerRank[rank] = NumberParticlesPerTypePerRank[rank][PartType];

      HaloPartitioner.Fill(AllHaloSizesThisType, NumberParticlesThisTypePerRank);
    }

    /* Tell each MPI rank what the minimum and maximum halo number they contain */
    HBTInt FirstHaloThisRank, LastHaloThisRank;
    MPI_Scatter(HaloPartitioner.RankFirstHalo.data(), 1, MPI_HBT_INT, &FirstHaloThisRank, 1, MPI_HBT_INT, root_node,
                world.Communicator);
    MPI_Scatter(HaloPartitioner.RankLastHalo.data(), 1, MPI_HBT_INT, &LastHaloThisRank, 1, MPI_HBT_INT, root_node,
                world.Communicator);

    /* The root rank will send the size of each halo segment contained in non-root
     * MPI ranks. */
    std::vector<HBTInt> LocalHaloSizesThisType;
    if (world.rank() == root_node)
    {
      for (int rank = 0; rank < world.size(); rank++)
      {
        if (rank == root_node)
          LocalHaloSizesThisType = HaloPartitioner.HaloSizesOnProc[rank];
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
      LocalHaloSizesThisType.resize(nhalos);
      MPI_Recv(LocalHaloSizesThisType.data(), nhalos, MPI_HBT_INT, root_node, 0, world.Communicator, MPI_STATUS_IGNORE);
    }

    /* Retrieve the particle offsets of each halo segment, measured relative to
     * the first particle in each MPI rank. */
    std::vector<HBTInt> LocalHaloOffsetsThisType;
    HBTInt np = CompileOffsets(LocalHaloSizesThisType, LocalHaloOffsetsThisType);
    LocalHaloOffsetsThisType.push_back(np);

    /* We only want to iterate over particles of this type. */
    auto FirstParticleThisType = Particles.begin() + std::accumulate(NumberParticlesPerTypeThisRank.begin(), NumberParticlesPerTypeThisRank.begin() + PartType, 0);
    auto LastParticleThisType  = Particles.begin() + std::accumulate(NumberParticlesPerTypeThisRank.begin(), NumberParticlesPerTypeThisRank.begin() + PartType + 1, 0);

#pragma omp parallel for default(shared)
    for (size_t halo_i = 0; halo_i < LocalHaloSizesThisType.size(); halo_i++)
    {
      /* The range of particles that belong to the current halo segment. */
      auto start_particle = FirstParticleThisType + LocalHaloOffsetsThisType[halo_i];
      auto end_particle = FirstParticleThisType + LocalHaloOffsetsThisType[halo_i + 1];

      /* We assign a HostHalo to each particle based on its host. This needs to be done because
       * when we find host haloes we need up-to-date information for particles that are part of
       * subhaloes. */
      for (auto particle_it = start_particle; particle_it != end_particle; ++particle_it)
      {
        particle_it->HostId = halo_i + FirstHaloThisRank;
        LocalNumberParticlesAssigned++;
      }
    }
  }

  /* Sanity check */
  HBTInt TotalNumberParticlesAssigned = 0;
  MPI_Reduce(&LocalNumberParticlesAssigned, &TotalNumberParticlesAssigned, 1, MPI_HBT_INT, MPI_SUM, root_node, world.Communicator);

  /* TODO: Add loading of TotalNumberGroupParticles to enable assert, since we
   * are currently loading it at a later stage and this will trigger a memory
   * error. */
  if (world.rank() == root_node)
    assert(TotalNumberParticlesAssigned == TotalNumberGroupParticles);
}

/* Group up particles in the same rank that share the same HostId, and create
 * a new halo segment within Halos. */
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
