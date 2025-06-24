using namespace std;
#include <iostream>
#include <numeric>
#include <cassert>
#include <assert.h>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <list>
#include <sstream>
#include <string>
#include <typeinfo>
#include <stdexcept>
#include <glob.h>
#include "../../config_parser.h"
#include "../../hdf_wrapper.h"
#include "../../mymath.h"
#include "../../snapshot.h"
#include "bonsai_io.h"
#include "../exchange_and_merge.h"

void create_TipsyHeader_MPI_type(MPI_Datatype &dtype)
{
  /*to create the struct data type for communication*/
  TipsyHeader_t p;
#define NumAttr 18
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

  RegisterAttr(Time, MPI_DOUBLE, 1);
  RegisterAttr(NumPart, MPI_INT, 1);
  RegisterAttr(NumDimensions, MPI_INT, 1);
  RegisterAttr(NumSPHPart, MPI_INT, 1);
  RegisterAttr(NumDarkMatterParticles, MPI_INT, 1);
  RegisterAttr(NumStellarParticles, MPI_INT, 1);
  RegisterAttr(Version, MPI_INT, 1);

#undef RegisterAttr
  assert(i <= NumAttr);

  MPI_Type_create_struct(i, blockcounts, offsets, oldtypes, &dtype);
  MPI_Type_create_resized(dtype, (MPI_Aint)0, extent, &dtype);
  MPI_Type_commit(&dtype);
#undef NumAttr
}

/* Uses glob to return the file names that match the string pattern */
std::vector<std::string> GlobVector(const std::string &pattern)
{
  glob_t GlobResult;
  glob(pattern.c_str(), GLOB_TILDE, NULL, &GlobResult);

  std::vector<std::string> files;
  for(int i = 0; i<GlobResult.gl_pathc; i++)
  {
    files.push_back(string(GlobResult.gl_pathv[i]));
  }
  globfree(&GlobResult);

  return files;
}

void BonsaiSimReader_t::SetSnapshotFileName(int snapshotId)
{
  if(HBTConfig.SnapshotNameList.empty())
  {
    /* Since BONSAI names files according to output time, we will use glob to
     * get their path. */
    stringstream formatter;
    formatter << HBTConfig.SnapshotPath << "/" << HBTConfig.SnapshotFileBase << "_*";

    std::vector<std::string> SnapshotNameList = GlobVector(formatter.str());
    std::sort(SnapshotNameList.begin(), SnapshotNameList.end());
    SnapshotName = SnapshotNameList[snapshotId];
  }
  else
  {
    SnapshotName = HBTConfig.SnapshotNameList[snapshotId];
  }
}

void BonsaiSimReader_t::GetSnapshotFileName(std::string &filename)
{
  filename = SnapshotName;
}

void BonsaiSimReader_t::OpenFile(std::ifstream &file)
{
  /* Get the file name */
  string filename;
  GetSnapshotFileName(filename);

  /* Open in binary mode */
  file.open(filename, std::ios::in | std::ios::binary);

  if (file.fail() < 0)
  {
    cout << "Failed to open file: " << filename << "\n";
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
}

void BonsaiSimReader_t::ReadHeader(TipsyHeader_t &Header)
{
  /* Need to pass by reference to handle binary inputs */
  std::ifstream file;
  OpenFile(file);

  /* Read tipsy header */
  file.read((char*)&Header, sizeof(Header));
}

void BonsaiSimReader_t::GetParticleCountInFile(hid_t file, int np[])
{
  int NumPartTypes;
  ReadAttribute(file, "Header", "NumPartTypes", H5T_NATIVE_INT, &NumPartTypes);
  int NumPart_ThisFile[NumPartTypes]; // same size as attributes in the file
  ReadAttribute(file, "Header", "NumPart_ThisFile", H5T_NATIVE_INT, NumPart_ThisFile);
  for (int i = 0; i < TypeMax; i++)
    if (i < NumPartTypes)
      np[i] = NumPart_ThisFile[i];
    else
      np[i] = 0;
#ifdef DM_ONLY
  for (int i = 0; i < TypeMax; i++)
    if (i != TypeDM)
      np[i] = 0;
#endif
}

HBTInt BonsaiSimReader_t::CompileFileOffsets(int nfiles)
{
  HBTInt offset = 0;
  np_file.reserve(nfiles);
  offset_file.reserve(nfiles);

  for (int ifile = 0; ifile < nfiles; ifile++)
  {
    /* We already have the particle count from the header (single file). */
    offset_file.push_back(offset);
    HBTInt np = Header.NumPart;

    np_file.push_back(np);
    offset += np;
  }
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

void BonsaiSimReader_t::ReadSnapshot(int ifile, Particle_t *ParticlesInFile, HBTInt file_start, HBTInt file_count)
{
  hid_t file = OpenFile(ifile);
  vector<int> np_this(TypeMax);
  vector<HBTInt> offset_this(TypeMax);
  GetParticleCountInFile(file, np_this.data());
  CompileOffsets(np_this, offset_this);

  HBTReal boxsize = Header.BoxSize;
  auto ParticlesToRead = ParticlesInFile;
  for (int itype = 0; itype < TypeMax; itype++)
  {

    // Find the range of offsets in the file for particles of this type
    HBTInt type_first_offset = offset_this[itype];
    HBTInt type_last_offset = type_first_offset + np_this[itype] - 1;

    // Find the range of offsets in the file we actually want to read
    HBTInt read_first_offset = file_start;
    HBTInt read_last_offset = file_start + file_count - 1;

    // The overlap of these two ranges contains the particles we will read now.
    HBTInt i1 = type_first_offset;
    if (read_first_offset > i1)
      i1 = read_first_offset;
    HBTInt i2 = type_last_offset;
    if (read_last_offset < i2)
      i2 = read_last_offset;

    // Compute range of particles of this type to read from this file
    HBTInt read_offset = i1 - offset_this[itype];
    HBTInt read_count = i2 - i1 + 1;
    if (read_count <= 0)
      continue;
    assert(read_offset >= 0);
    assert(read_offset + read_count <= np_this[itype]);

    // Open the HDF5 group for this type
    stringstream grpname;
    grpname << "PartType" << itype;
    hid_t particle_data = H5Gopen2(file, grpname.str().c_str(), H5P_DEFAULT);
    check_id_size(particle_data);

    const hsize_t chunksize = 10 * 1024 * 1024;

    // Positions
    {
      // Check that positions are comoving
      HBTReal aexp;
      ReadAttribute(particle_data, "Coordinates", "a-scale exponent", H5T_HBTReal, &aexp);
      if (aexp != 1.0)
      {
        cout << "Can't handle Coordinates with a-scale exponent != 1\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
      }

      // Read data in chunks to minimize memory overhead
      for (hsize_t offset = 0; offset < read_count; offset += chunksize)
      {
        // Read the next chunk
        hsize_t count = read_count - offset; // number left to read
        if (count > chunksize)
          count = chunksize;
        vector<HBTxyz> x(count);
        ReadPartialDataset(particle_data, "Coordinates", H5T_HBTReal, x.data(), offset + read_offset, count);
        // Box wrap if necessary
        if (HBTConfig.PeriodicBoundaryOn)
        {
          for (hsize_t i = 0; i < count; i++)
            for (int j = 0; j < 3; j++)
              x[i][j] = position_modulus(x[i][j], boxsize);
        }
        // Store the particle positions
        for (hsize_t i = 0; i < count; i += 1)
          for (int j = 0; j < 3; j += 1)
            ParticlesToRead[offset + i].ComovingPosition[j] = x[i][j];
      }
    }

    // Velocities
    {
      HBTReal aexp;
      ReadAttribute(particle_data, "Velocities", "a-scale exponent", H5T_HBTReal, &aexp);

      // Read data in chunks to minimize memory overhead
      for (hsize_t offset = 0; offset < read_count; offset += chunksize)
      {
        // Read the next chunk
        hsize_t count = read_count - offset;
        if (count > chunksize)
          count = chunksize;
        vector<HBTxyz> v(count);
        ReadPartialDataset(particle_data, "Velocities", H5T_HBTReal, v.data(), offset + read_offset, count);
        // Convert units and store the particle velocities
        for (hsize_t i = 0; i < count; i += 1)
          for (int j = 0; j < 3; j += 1)
            ParticlesToRead[offset + i].PhysicalVelocity[j] = v[i][j] * pow(Header.ScaleFactor, aexp);
      }
    }

    // Ids
    {
      for (hsize_t offset = 0; offset < read_count; offset += chunksize)
      {
        hsize_t count = read_count - offset;
        if (count > chunksize)
          count = chunksize;
        vector<HBTInt> id(count);
        ReadPartialDataset(particle_data, "ParticleIDs", H5T_HBTInt, id.data(), offset + read_offset, count);
        for (hsize_t i = 0; i < count; i += 1)
          ParticlesToRead[offset + i].Id = id[i];
      }
    }

    // Masses
    {
      HBTReal aexp;
      std::string name;
      if (itype == 5)
        name = "DynamicalMasses";
      else
        name = "Masses";
      ReadAttribute(particle_data, name.c_str(), "a-scale exponent", H5T_HBTReal, &aexp);
      for (hsize_t offset = 0; offset < read_count; offset += chunksize)
      {
        hsize_t count = read_count - offset;
        if (count > chunksize)
          count = chunksize;
        vector<HBTReal> m(count);
        ReadPartialDataset(particle_data, name.c_str(), H5T_HBTReal, m.data(), offset + read_offset, count);
        for (hsize_t i = 0; i < count; i += 1)
          ParticlesToRead[offset + i].Mass = m[i] * pow(Header.ScaleFactor, aexp);
      }
    }

#ifndef DM_ONLY
    // internal energy
#ifdef HAS_THERMAL_ENERGY
    if (itype == 0)
    {
      HBTReal aexp;
      ReadAttribute(particle_data, "InternalEnergies", "a-scale exponent", H5T_HBTReal, &aexp);
      for (hsize_t offset = 0; offset < read_count; offset += chunksize)
      {
        hsize_t count = read_count - offset;
        if (count > chunksize)
          count = chunksize;
        vector<HBTReal> u(count);
        ReadPartialDataset(particle_data, "InternalEnergies", H5T_HBTReal, u.data(), offset + read_offset, count);
        for (hsize_t i = 0; i < count; i += 1)
          ParticlesToRead[offset + i].InternalEnergy = u[i] * pow(Header.ScaleFactor, aexp);
      }
    }
    else
    {
      // Zero out internal energy for non-gas particles
      for (hsize_t offset = 0; offset < read_count; offset += 1)
        ParticlesToRead[offset].InternalEnergy = 0.0;
    }
#endif
    { // type
      ParticleType_t t = static_cast<ParticleType_t>(itype);
      for (hsize_t i = 0; i < read_count; i++)
        ParticlesToRead[i].Type = t;
    }
#endif

    // Hostid
    {
      for (hsize_t offset = 0; offset < read_count; offset += chunksize)
      {
        hsize_t count = read_count - offset;
        if (count > chunksize)
          count = chunksize;
        vector<HBTInt> id(count);
        ReadPartialDataset(particle_data, "FOFGroupIDs", H5T_HBTInt, id.data(), offset + read_offset, count);
        for (hsize_t i = 0; i < count; i += 1)
          ParticlesToRead[offset + i].HostId = id[i];
      }
    }

    // Advance to next particle type
    H5Gclose(particle_data);
    ParticlesToRead += read_count;
  }
  H5Fclose(file);
}

void BonsaiSimReader_t::ReadGroupParticles(int ifile, Particle_t *ParticlesInFile, HBTInt file_start, HBTInt file_count,
                                          bool FlagReadParticleId)
{
  hid_t file = OpenFile(ifile);
  vector<int> np_this(TypeMax);
  vector<HBTInt> offset_this(TypeMax);
  GetParticleCountInFile(file, np_this.data());
  CompileOffsets(np_this, offset_this);

  HBTReal boxsize = Header.BoxSize;
  auto ParticlesToRead = ParticlesInFile;
  for (int itype = 0; itype < TypeMax; itype++)
  {
    // Find the range of offsets in the file for particles of this type
    HBTInt type_first_offset = offset_this[itype];
    HBTInt type_last_offset = type_first_offset + np_this[itype] - 1;

    // Find the range of offsets in the file we actually want to read
    HBTInt read_first_offset = file_start;
    HBTInt read_last_offset = file_start + file_count - 1;

    // The overlap of these two ranges contains the particles we will read now.
    HBTInt i1 = type_first_offset;
    if (read_first_offset > i1)
      i1 = read_first_offset;
    HBTInt i2 = type_last_offset;
    if (read_last_offset < i2)
      i2 = read_last_offset;

    // Compute range of particles of this type to read from this file
    HBTInt read_offset = i1 - offset_this[itype];
    HBTInt read_count = i2 - i1 + 1;
    if (read_count <= 0)
      continue;
    assert(read_offset >= 0);
    assert(read_offset + read_count <= np_this[itype]);

    // Open the HDF5 group for this particle type
    stringstream grpname;
    grpname << "PartType" << itype;
    hid_t particle_data = H5Gopen2(file, grpname.str().c_str(), H5P_DEFAULT);

    const hsize_t chunksize = 10 * 1024 * 1024;

    if (FlagReadParticleId)
    {

      // Positions
      {
        // Check that positions are comoving
        HBTReal aexp;
        ReadAttribute(particle_data, "Coordinates", "a-scale exponent", H5T_HBTReal, &aexp);
        if (aexp != 1.0)
        {
          cout << "Can't handle Coordinates with a-scale exponent != 1\n";
          MPI_Abort(MPI_COMM_WORLD, 1);
        }

        // Read data in chunks to minimize memory overhead
        for (hsize_t offset = 0; offset < read_count; offset += chunksize)
        {
          // Read the next chunk
          hsize_t count = read_count - offset;
          if (count > chunksize)
            count = chunksize;
          vector<HBTxyz> x(count);
          ReadPartialDataset(particle_data, "Coordinates", H5T_HBTReal, x.data(), offset + read_offset, count);
          // Box wrap if necessary
          if (HBTConfig.PeriodicBoundaryOn)
          {
            for (hsize_t i = 0; i < count; i++)
              for (int j = 0; j < 3; j++)
                x[i][j] = position_modulus(x[i][j], boxsize);
          }
          // Store the particle positions
          for (hsize_t i = 0; i < count; i += 1)
            for (int j = 0; j < 3; j += 1)
              ParticlesToRead[offset + i].ComovingPosition[j] = x[i][j];
        }
      }

      // Velocities
      {
        HBTReal aexp;
        ReadAttribute(particle_data, "Velocities", "a-scale exponent", H5T_HBTReal, &aexp);

        // Read data in chunks to minimize memory overhead
        for (hsize_t offset = 0; offset < read_count; offset += chunksize)
        {
          // Read the next chunk
          hsize_t count = read_count - offset;
          if (count > chunksize)
            count = chunksize;
          vector<HBTxyz> v(count);
          ReadPartialDataset(particle_data, "Velocities", H5T_HBTReal, v.data(), offset + read_offset, count);
          // Convert units and store the particle velocities
          for (hsize_t i = 0; i < count; i += 1)
            for (int j = 0; j < 3; j += 1)
              ParticlesToRead[offset + i].PhysicalVelocity[j] = v[i][j] * pow(Header.ScaleFactor, aexp);
        }
      }

      // Ids
      {
        for (hsize_t offset = 0; offset < read_count; offset += chunksize)
        {
          hsize_t count = read_count - offset;
          if (count > chunksize)
            count = chunksize;
          vector<HBTInt> id(count);
          ReadPartialDataset(particle_data, "ParticleIDs", H5T_HBTInt, id.data(), offset + read_offset, count);
          for (hsize_t i = 0; i < count; i += 1)
            ParticlesToRead[offset + i].Id = id[i];
        }
      }

      // Masses
      {
        HBTReal aexp;
        std::string name;
        if (itype == 5)
          name = "DynamicalMasses";
        else
          name = "Masses";
        ReadAttribute(particle_data, name.c_str(), "a-scale exponent", H5T_HBTReal, &aexp);
        for (hsize_t offset = 0; offset < read_count; offset += chunksize)
        {
          hsize_t count = read_count - offset;
          if (count > chunksize)
            count = chunksize;
          vector<HBTReal> m(count);
          ReadPartialDataset(particle_data, name.c_str(), H5T_HBTReal, m.data(), offset + read_offset, count);
          for (hsize_t i = 0; i < count; i += 1)
            ParticlesToRead[offset + i].Mass = m[i] * pow(Header.ScaleFactor, aexp);
        }
      }

#ifndef DM_ONLY
      // internal energy
#ifdef HAS_THERMAL_ENERGY
      if (itype == 0)
      {
        HBTReal aexp;
        ReadAttribute(particle_data, "InternalEnergies", "a-scale exponent", H5T_HBTReal, &aexp);
        for (hsize_t offset = 0; offset < read_count; offset += chunksize)
        {
          hsize_t count = read_count - offset;
          if (count > chunksize)
            count = chunksize;
          vector<HBTReal> u(count);
          ReadPartialDataset(particle_data, "InternalEnergies", H5T_HBTReal, u.data(), offset + read_offset, count);
          for (hsize_t i = 0; i < count; i += 1)
            ParticlesToRead[offset + i].InternalEnergy = u[i] * pow(Header.ScaleFactor, aexp);
        }
      }
#endif
      { // type
        ParticleType_t t = static_cast<ParticleType_t>(itype);
        for (int i = 0; i < read_count; i++)
          ParticlesToRead[i].Type = t;
      }
#endif
    }

    // Hostid
    {
      for (hsize_t offset = 0; offset < read_count; offset += chunksize)
      {
        hsize_t count = read_count - offset;
        if (count > chunksize)
          count = chunksize;
        vector<HBTInt> id(count);
        ReadPartialDataset(particle_data, "FOFGroupIDs", H5T_HBTInt, id.data(), offset + read_offset, count);
        for (hsize_t i = 0; i < count; i += 1)
          ParticlesToRead[offset + i].HostId = id[i];
      }
    }
    // Advance to next particle type
    ParticlesToRead += read_count;
    H5Gclose(particle_data);
  }
  H5Fclose(file);
}

void BonsaiSimReader_t::LoadSnapshot(MpiWorker_t &world, int snapshotId, vector<Particle_t> &Particles,
                                    Cosmology_t &Cosmology)
{

  MPI_Barrier(world.Communicator);

  // Decide how many ranks per node read simultaneously
  int nr_nodes = (world.size() / world.MaxNodeSize);
  int nr_reading = HBTConfig.MaxConcurrentIO / nr_nodes;
  if (nr_reading < 1)
    nr_reading = 1; // Always at least one per node

  SetSnapshotFileName(snapshotId);

  const int root = 0;
  if (world.rank() == root)
  {
    ReadHeader(Header);
    CompileFileOffsets(1);
  }

  MPI_Bcast(&Header, 1, MPI_TipsyHeader_t, root, world.Communicator);
  world.SyncContainer(np_file, MPI_HBT_INT, root);
  world.SyncContainer(offset_file, MPI_HBT_INT, root);

  /* Update the tree parameters based on the softenings */
  HBTConfig.TreeNodeResolution = HBTConfig.SofteningHalo * 0.1;
  HBTConfig.TreeNodeResolutionHalf = HBTConfig.TreeNodeResolution / 2.;

  /* Only dealing with virial units and non-cosmological simulations  */
  PhysicalConst::G = 1;
  PhysicalConst::H0 = 0;

  /* Since H0 = 0, there is no evolution in the Hubble Parameter. We set these
   * values to -1 to signal they are dummy variables. */
  Cosmology.Set_Hz(-1);
  Cosmology.Set(-1, -1, -1);

  // Decide how many particles this MPI rank will read
  HBTInt np_total = accumulate(np_file.begin(), np_file.end(), (HBTInt)0);
  HBTInt np_local = np_total / world.size();
  if (world.rank() < (np_total % world.size()))
    np_local += 1;
#ifndef NDEBUG
  HBTInt np_check;
  MPI_Allreduce(&np_local, &np_check, 1, MPI_HBT_INT, MPI_SUM, world.Communicator);
  assert(np_check == np_total);
#endif

  // Determine offset to the first and last particle this rank will read
  HBTInt local_first_offset;
  MPI_Scan(&np_local, &local_first_offset, 1, MPI_HBT_INT, MPI_SUM, world.Communicator);
  local_first_offset -= np_local;
  HBTInt local_last_offset = local_first_offset + np_local - 1;
  assert(local_first_offset >= 0);
  assert(local_last_offset < np_total);

  // Allocate storage for the particles
  Particles.resize(np_local);

  // Allow a limited number of ranks per node to read simultaneously
  int reads_done = 0;
  for (int rank_within_node = 0; rank_within_node < world.MaxNodeSize; rank_within_node += 1)
  {
    if (rank_within_node == world.NodeRank)
    {

      // Loop over all files
      HBTInt particle_offset = 0;
      for (int file_nr = 0; file_nr < Header.NumberOfFiles; file_nr += 1)
      {

        // Determine global offset of first particle to read from this file:
        // This is the larger of the offset of the first particle in the file
        // and the offset of the first particle this rank is to read.
        HBTInt i1 = offset_file[file_nr];
        if (local_first_offset > i1)
          i1 = local_first_offset;

        // Determine global offset of last particle to read from this file:
        // This is the smaller of the offset to the last particle in this file
        // and the offset of the last particle this rank is to read.
        HBTInt i2 = offset_file[file_nr] + np_file[file_nr] - 1;
        if (local_last_offset < i2)
          i2 = local_last_offset;

        if (i2 >= i1)
        {
          // We have particles to read from this file.
          HBTInt file_start = i1 - offset_file[file_nr]; // Offset to first particle to read
          HBTInt file_count = i2 - i1 + 1;               // Number of particles to read
          assert(file_count > 0);
          assert(file_start >= 0);
          assert(file_start + file_count <= np_file[file_nr]);
          ReadSnapshot(file_nr, Particles.data() + particle_offset, file_start, file_count);
          particle_offset += file_count;
        }
      }                                    // Next file
      assert(particle_offset == np_local); // Check we read the expected number of particles
      reads_done += 1;
    }
    if (rank_within_node % nr_reading == nr_reading - 1)
      MPI_Barrier(world.Communicator);
  } // Next MPI rank within the node

  // Every rank should have executed the reading code exactly once
  assert(reads_done == 1);

  global_timer.Tick("snap_io", world.Communicator);

#ifdef SNAPSHOT_IO_TEST
  // For testing: dump the snapshot to a new set of files
  // Generate test file name for this MPI  rank
  stringstream formatter1;
  formatter1 << HBTConfig.SubhaloPath << "/" << setw(3) << setfill('0') << snapshotId << "/"
             << "test_" << setw(3) << setfill('0') << snapshotId << "." << world.rank() << ".hdf5";
  string tfilename = formatter1.str();
  // Create array of coordinates
  double *pos = (double *)malloc(3 * sizeof(double) * np_local);
  for (size_t i = 0; i < np_local; i += 1)
  {
    pos[3 * i + 0] = Particles[i].ComovingPosition[0];
    pos[3 * i + 1] = Particles[i].ComovingPosition[1];
    pos[3 * i + 2] = Particles[i].ComovingPosition[2];
  }
  // Create array of IDs
  long long *ids = (long long *)malloc(sizeof(long long) * np_local);
  for (size_t i = 0; i < np_local; i += 1)
    ids[i] = Particles[i].Id;
  // Create array of types
  int *type = (int *)malloc(sizeof(int) * np_local);
  for (size_t i = 0; i < np_local; i += 1)
    type[i] = Particles[i].Type;

  // Create the file
  hid_t tfile = H5Fcreate(tfilename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  // Write out the data
  hsize_t ndims;
  hsize_t dims[2];
  ndims = 2;
  dims[0] = np_local;
  dims[1] = 3;
  writeHDFmatrix(tfile, pos, "Coordinates", ndims, dims, H5T_NATIVE_DOUBLE);
  ndims = 1;
  writeHDFmatrix(tfile, ids, "ParticleIDs", ndims, dims, H5T_NATIVE_LLONG);
  writeHDFmatrix(tfile, type, "Types", ndims, dims, H5T_NATIVE_INT);

  // Tidy up
  H5Fclose(tfile);
  free(pos);
  free(ids);
  free(type);

#endif
}

inline bool CompParticleHost(const Particle_t &a, const Particle_t &b)
{
  return a.HostId < b.HostId;
}

void BonsaiSimReader_t::LoadGroups(MpiWorker_t &world, int snapshotId, vector<Halo_t> &Halos)
{ // read in particle properties at the same time, to avoid particle look-up at later stage.
  SetSnapshotFileName(snapshotId);

  // Decide how many ranks per node read simultaneously
  int nr_nodes = (world.size() / world.MaxNodeSize);
  int nr_reading = HBTConfig.MaxConcurrentIO / nr_nodes;
  if (nr_reading < 1)
    nr_reading = 1; // Always at least one per node

  const int root = 0;
  if (world.rank() == root)
  {
    ReadHeader(Header);
    CompileFileOffsets(1);
  }
  MPI_Bcast(&Header, 1, MPI_TipsyHeader_t, root, world.Communicator);
  world.SyncContainer(np_file, MPI_HBT_INT, root);
  world.SyncContainer(offset_file, MPI_HBT_INT, root);

  // Decide how many particles this MPI rank will read
  HBTInt np_total = accumulate(np_file.begin(), np_file.end(), (HBTInt)0);
  HBTInt np_local = np_total / world.size();
  if (world.rank() < (np_total % world.size()))
    np_local += 1;
#ifndef NDEBUG
  HBTInt np_check;
  MPI_Allreduce(&np_local, &np_check, 1, MPI_HBT_INT, MPI_SUM, world.Communicator);
  assert(np_check == np_total);
#endif

  // Determine offset to the first and last particle this rank will read
  HBTInt local_first_offset;
  MPI_Scan(&np_local, &local_first_offset, 1, MPI_HBT_INT, MPI_SUM, world.Communicator);
  local_first_offset -= np_local;
  HBTInt local_last_offset = local_first_offset + np_local - 1;
  assert(local_first_offset >= 0);
  assert(local_last_offset < np_total);

  // Allocate storage for the particles
  vector<Particle_t> ParticleHosts;
  ParticleHosts.resize(np_local);

  bool FlagReadId = true;

  // Allow a limited number of ranks per node to read simultaneously
  int reads_done = 0;
  for (int rank_within_node = 0; rank_within_node < world.MaxNodeSize; rank_within_node += 1)
  {
    if (rank_within_node == world.NodeRank)
    {

      // Loop over all files
      HBTInt particle_offset = 0;
      for (int file_nr = 0; file_nr < Header.NumberOfFiles; file_nr += 1)
      {

        // Determine global offset of first particle to read from this file:
        // This is the larger of the offset of the first particle in the file
        // and the offset of the first particle this rank is to read.
        HBTInt i1 = offset_file[file_nr];
        if (local_first_offset > i1)
          i1 = local_first_offset;

        // Determine global offset of last particle to read from this file:
        // This is the smaller of the offset to the last particle in this file
        // and the offset of the last particle this rank is to read.
        HBTInt i2 = offset_file[file_nr] + np_file[file_nr] - 1;
        if (local_last_offset < i2)
          i2 = local_last_offset;

        if (i2 >= i1)
        {
          // We have particles to read from this file.
          HBTInt file_start = i1 - offset_file[file_nr]; // Offset to first particle to read
          HBTInt file_count = i2 - i1 + 1;               // Number of particles to read
          assert(file_count > 0);
          assert(file_start >= 0);
          assert(file_start + file_count <= np_file[file_nr]);
          ReadGroupParticles(file_nr, ParticleHosts.data() + particle_offset, file_start, file_count, FlagReadId);
          particle_offset += file_count;
        }
      }                                    // Next file
      assert(particle_offset == np_local); // Check we read the expected number of particles
      reads_done += 1;
    }
    if (rank_within_node % nr_reading == nr_reading - 1)
      MPI_Barrier(world.Communicator);
  } // Next MPI rank within the node

  // Every rank should have executed the reading code exactly once
  assert(reads_done == 1);

  global_timer.Tick("halo_io", world.Communicator);

  // #define HALO_IO_TEST
#ifdef HALO_IO_TEST
  //
  // For testing: dump the snapshot to a new set of files
  //
  // Generate test file name for this MPI  rank
  stringstream formatter1;
  formatter1 << HBTConfig.SubhaloPath << "/" << setw(3) << setfill('0') << snapshotId << "/"
             << "test_halo_" << setw(3) << setfill('0') << snapshotId << "." << world.rank() << ".hdf5";
  string tfilename = formatter1.str();
  // Create array of coordinates
  double *pos = (double *)malloc(3 * sizeof(double) * np_local);
  for (size_t i = 0; i < np_local; i += 1)
  {
    pos[3 * i + 0] = ParticleHosts[i].ComovingPosition[0];
    pos[3 * i + 1] = ParticleHosts[i].ComovingPosition[1];
    pos[3 * i + 2] = ParticleHosts[i].ComovingPosition[2];
  }
  // Create array of IDs
  long long *ids = (long long *)malloc(sizeof(long long) * np_local);
  for (size_t i = 0; i < np_local; i += 1)
    ids[i] = ParticleHosts[i].Id;
  // Create array of types
  int *type = (int *)malloc(sizeof(int) * np_local);
  for (size_t i = 0; i < np_local; i += 1)
    type[i] = ParticleHosts[i].Type;
  // Create array of group indexes
  int *fofnr = (int *)malloc(sizeof(int) * np_local);
  for (size_t i = 0; i < np_local; i += 1)
    fofnr[i] = ParticleHosts[i].HostId;

  // Create the file
  hid_t tfile = H5Fcreate(tfilename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  // Write out the data
  hsize_t ndims;
  hsize_t dims[2];
  ndims = 2;
  dims[0] = np_local;
  dims[1] = 3;
  writeHDFmatrix(tfile, pos, "Coordinates", ndims, dims, H5T_NATIVE_DOUBLE);
  ndims = 1;
  writeHDFmatrix(tfile, ids, "ParticleIDs", ndims, dims, H5T_NATIVE_LLONG);
  writeHDFmatrix(tfile, type, "Types", ndims, dims, H5T_NATIVE_INT);
  writeHDFmatrix(tfile, fofnr, "FoFNr", ndims, dims, H5T_NATIVE_INT);

  // Tidy up
  H5Fclose(tfile);
  free(pos);
  free(ids);
  free(type);
  free(fofnr);
  //
  // END OF TEST CODE
  //
#endif

  // Sort particles by host
  sort(ParticleHosts.begin(), ParticleHosts.end(), CompParticleHost);
  if (!ParticleHosts.empty())
  {
    assert(ParticleHosts.front().HostId >= 0);                 // min haloid>=0
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
  vector<HaloLen_t> HaloLen;

  HBTInt curr_host_id = -1;
  for (auto &&p : ParticleHosts)
  {
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
  auto p_in = ParticleHosts.begin();
  for (auto &&h : Halos)
  {
    for (auto &&p : h.Particles)
    {
      p = *p_in;
      ++p_in;
    }
  }

  VectorFree(ParticleHosts);

  ExchangeAndMerge(world, Halos);

  global_timer.Tick("halo_comms", world.Communicator);

  HBTConfig.GroupLoadedFullParticle = true;
}

bool IsBonsaiSimGroup(const string &GroupFileFormat)
{
  return GroupFileFormat.substr(0, 6) == "bonsai";
}