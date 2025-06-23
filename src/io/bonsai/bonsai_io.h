/* IO for BonsaiSim data.
 *
 * To specify a list of snapshot, list the snapshot directories (one per line) in snapshotlist.txt and place it under
 * your subhalo output directory.
 *
 * To use this IO, in the config file, set SnapshotFormat to Bonsaisim, and set GroupFileFormat to Bonsaisim or
 * Bonsaisim_particle_index.
 *
 * The groups loaded are already filled with particle properties, and the halos are distributed to processors according
 * to the CoM of each halo.
 */

#ifndef BONSAI_IO_INCLUDED
#define BONSAI_IO_INCLUDED
#include "../../halo.h"
#include "../../hdf_wrapper.h"
#include "../../mpi_wrapper.h"

struct BonsaiSimHeader_t
{
  int NumberOfFiles;
  double BoxSize;
  double ScaleFactor;
  double OmegaM0;
  double OmegaLambda0;
  double Hz;
  double h;
  double mass[TypeMax];
  int npart[TypeMax];
  HBTInt npartTotal[TypeMax];
  double MassInMsunh;
  double LengthInMpch;
  double VelInKmS;
  HBTInt NullGroupId;
  double DM_comoving_softening;
  double DM_maximum_physical_softening;
};

void create_BonsaiSimHeader_MPI_type(MPI_Datatype &dtype);

class BonsaiSimReader_t
{
  string SnapshotName;
  vector<HBTInt> np_file;
  vector<HBTInt> offset_file;
  BonsaiSimHeader_t Header;
  void OpenFile(std::ifstream &file);
  void ReadHeader(int ifile, BonsaiSimHeader_t &header);
  void ReadUnits(HBTReal &MassInMsunh, HBTReal &LengthInMpch, HBTReal &VelInKmS);
  HBTInt CompileFileOffsets(int nfiles);
  void ReadSnapshot(int ifile, Particle_t *ParticlesInFile, HBTInt file_start, HBTInt file_count);
  void ReadGroupParticles(int ifile, Particle_t *ParticlesInFile, HBTInt file_start, HBTInt file_count,
                          bool FlagReadParticleId);
  void GetSnapshotFileName(std::string &filename);
  void SetSnapshotFileName(int snapshotId);
  void GetParticleCountInFile(hid_t file, int np[]);

  MPI_Datatype MPI_BonsaiSimHeader_t;

public:
  BonsaiSimReader_t()
  {
    create_BonsaiSimHeader_MPI_type(MPI_BonsaiSimHeader_t);
  }
  ~BonsaiSimReader_t()
  {
    MPI_Type_free(&MPI_BonsaiSimHeader_t);
  }
  void LoadSnapshot(MpiWorker_t &world, int snapshotId, vector<Particle_t> &Particles, Cosmology_t &Cosmology);
  void LoadGroups(MpiWorker_t &world, int snapshotId, vector<Halo_t> &Halos);
};

extern bool IsBonsaiSimGroup(const string &GroupFileFormat);
#endif
