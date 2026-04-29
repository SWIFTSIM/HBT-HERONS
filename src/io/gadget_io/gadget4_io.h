/* IO for GADGET4 hdf data.
 *
 * Currently only designed and tested for DM-only simulations; support for
 * multiple particle types to be added.
 *
 * To use this IO, in the config file, set SnapshotFormat and GroupFileFormat to
 * gadget4.
 *
 * The groups loaded are already filled with particle properties, and the halos
 * are distributed to processors according to the CoM of each halo. */

#ifndef GADGET4_IO_INCLUDED
#define GADGET4_IO_INCLUDED
#include "../../hdf_wrapper.h"
#include "../../halo.h"
#include "../../mpi_wrapper.h"

namespace Gadget4Reader
{
struct Gadget4Header_t
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

void create_Gadget4Header_MPI_type(MPI_Datatype &dtype);

class Gadget4Reader_t
{
  int SnapshotId;
  string SnapshotName;
  const HBTInt NullGroupId = std::numeric_limits<HBTInt>::max();
  const int root_node = 0;

  /* Header and misc */
  Gadget4Header_t Header;
  void ReadHeader(int ifile, Gadget4Header_t &header);

  /* Particle information. */
  std::vector<HBTInt> NumberParticlesPerRank;
  std::vector<HBTInt> NumberParticlesPerFile;
  std::vector<HBTInt> NumberParticlesPerTypeThisRank;

  vector<HBTInt> offset_file;
  HBTInt CompileFileOffsets(int nfiles);
  void GetFileName(int ifile, string &filename);
  void SetSnapshot(int snapshotId);
  void GetParticleCountInFile(hid_t file, int np[]);
  void ReadSnapshot(int ifile, Particle_t *ParticlesInFile);
  void LoadParticleProperties(MpiWorker_t &world, vector<Particle_t> &Particles);

  /* FoF group information. */
  HBTInt TotalNumberGroupParticles;

  /* These two vectors are only significant in root MPI rank. */
  std::vector<HBTInt> AllHaloSizes;
  std::vector<std::array<HBTInt, TypeMax>> AllHaloSizesPerType;

  void GetNumberParticlesPerRank(MpiWorker_t &world, const std::vector<Particle_t> &Particles);
  void LoadHaloSizes(MpiWorker_t &world);
  int ReadGroupFileCounts(int ifile);
  void GetGroupFileName(int ifile, string &filename);

  void ReadGroupLen(int ifile, HBTInt *buf);
  void ReadGroupLenPerType(int ifile, std::array<HBTInt, TypeMax> *buf);

  HBTInt ReadGroupFileTotParticles(int ifile);
  HBTInt CompileGroupFileOffsets(vector<HBTInt> &nhalo_per_groupfile, vector<HBTInt> &offsethalo_per_groupfile);
  void LoadParticleHosts(MpiWorker_t &world, vector<Particle_t> &Particles);
  void ReadGroupLengthPerParticleType(int ifile, HBTInt *buf);
  void CreateHaloSegments(const std::vector<Particle_t> &SnapshotParticles, std::vector<Halo_t> &Halos);

  MPI_Datatype MPI_Gadget4Header_t;

public:

  Gadget4Reader_t()
  {
    create_Gadget4Header_MPI_type(MPI_Gadget4Header_t);
  }
  ~Gadget4Reader_t()
  {
    My_Type_free(&MPI_Gadget4Header_t);
  }

  void LoadSnapshot(MpiWorker_t &world, int snapshotId, vector<Particle_t> &Particles, Cosmology_t &Cosmology);
  void LoadGroups(MpiWorker_t &world, const ParticleSnapshot_t &partsnap, vector<Halo_t> &Halos);
};

extern bool IsGadget4Group(const string &GroupFileFormat);
} // namespace Gadget4Reader

#endif // GADGET4_IO_INCLUDED