/* IO for Arepo hdf5 data.
 *
 * Currently only designed and tested for DM-only simulations; support for
 * multiple particle types to be added.
 *
 * To use this IO, in the config file, set SnapshotFormat and GroupFileFormat to
 * arepo.
 *
 * The groups loaded are already filled with particle properties, and the halos
 * are distributed to processors according to the CoM of each halo. */

#ifndef AREPO_IO_INCLUDED
#define AREPO_IO_INCLUDED
#include "../../hdf_wrapper.h"
#include "../../halo.h"
#include "../../mpi_wrapper.h"

namespace ArepoReader
{

struct ArepoHeader_t
{
  int NumberOfFiles;
  double BoxSize;
  double ScaleFactor;
  double OmegaM0;
  double OmegaLambda0;
  double mass[TypeMax];
  int npart[TypeMax];
  HBTInt npartTotal[TypeMax];

  double MassInMsunh;
  double LengthInMpch;
  double VelInKmS;

  double DM_comoving_softening;
  double DM_maximum_physical_softening;
};

void create_ArepoHeader_MPI_type(MPI_Datatype &dtype);

class ArepoReader_t
{
  int SnapshotId;
  std::string SnapshotName;
  void SetSnapshot(int snapshotId);

  /* Constants. TODO: extend these constants to the rest of the I/O readers. */
  const HBTInt NullGroupId = std::numeric_limits<HBTInt>::max();
  const int root_node = 0;

  /* File paths and number of subfiles. */
  int ReadGroupFileCounts(int ifile);
  std::string GetGroupFileName(int ifile);
  std::string GetSnapshotFileName(int ifile);

  /* Header and miscellaneous options*/
  ArepoHeader_t Header;
  ArepoHeader_t ReadHeader(int ifile);

  /* Particle count information. */
  std::vector<HBTInt> NumberParticlesPerRank;
  std::vector<HBTInt> NumberParticlesPerFile;
  std::vector<HBTInt> OffsetParticlesPerFile;
  std::vector<HBTInt> NumberParticlesPerTypeThisRank;
  std::vector<std::array<HBTInt, TypeMax>> NumberParticlesPerTypePerRank;
  std::vector<std::array<HBTInt, TypeMax>> NumberParticlesPerTypePerFile;
  std::vector<std::array<HBTInt, TypeMax>> OffsetParticlesPerTypePerFile;

  void CompileSnapshotParticleOffsets(int NumberFiles);
  void GetParticleCountInFile(hid_t file, HBTInt np[]);
  void GetNumberParticlesPerRank(MpiWorker_t &world, const std::vector<Particle_t> &Particles);

  void LoadParticleProperties(MpiWorker_t &world, std::vector<Particle_t> &Particles);
  void ReadSnapshot(int FirstFileIndex, int CurrentFileIndex, Particle_t *ParticlesOfThisType);

  /* FoF group information. */
  HBTInt TotalNumberGroupParticles;

  /* These two vectors are only significant in root MPI rank. */
  std::vector<HBTInt> AllHaloSizes;
  std::vector<std::array<HBTInt, TypeMax>> AllHaloSizesPerType;


  /* Loading FOF halo information. */
  void LoadHaloSizes(MpiWorker_t &world);
  HBTInt CompileGroupFileOffsets(std::vector<HBTInt> &nhalo_per_groupfile, std::vector<HBTInt> &offsethalo_per_groupfile);
  void ReadGroupLen(int ifile, HBTInt *buf);
  void ReadGroupLenPerType(int ifile, std::array<HBTInt, TypeMax> *buf);
  HBTInt ReadGroupFileTotParticles(int ifile);
  void LoadParticleHosts(MpiWorker_t &world, std::vector<Particle_t> &Particles);
  void CreateHaloSegments(const std::vector<Particle_t> &SnapshotParticles, std::vector<Halo_t> &Halos);

  /* Custom MPI datatypes for easier communication. */
  MPI_Datatype MPI_HBT_INT_ARRAY, MPI_AREPO_HEADER;
  MPI_Datatype INT_ARRAY_MPI_datatype();
  MPI_Datatype AREPO_HEADER_MPI_datatype();

public:

  ArepoReader_t()
  {
    MPI_HBT_INT_ARRAY = INT_ARRAY_MPI_datatype();
    MPI_AREPO_HEADER = AREPO_HEADER_MPI_datatype();
  }

  ~ArepoReader_t()
  {
    My_Type_free(&MPI_HBT_INT_ARRAY);
    My_Type_free(&MPI_AREPO_HEADER);
  }

  void LoadGroups(MpiWorker_t &world, const ParticleSnapshot_t &partsnap, std::vector<Halo_t> &Halos);
  void LoadSnapshot(MpiWorker_t &world, int snapshotId, std::vector<Particle_t> &Particles, Cosmology_t &Cosmology);

};

extern bool IsArepoGroup(const std::string &GroupFileFormat);
} // namespace ArepoReader

#endif // AREPO_IO_INCLUDED