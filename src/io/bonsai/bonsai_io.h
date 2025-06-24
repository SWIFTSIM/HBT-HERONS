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

struct TipsyHeader_t {
  double Time;
  int NumPart;
  int NumDimensions;
  int NumSPHPart;
  int NumDarkMatterParticles;
  int NumStellarParticles;
  int Version;
};

struct TipsyDarkMatterParticle_t {
  float Mass;
  float Position[3];
  float Velocity[3];

  /* If I do long long _ID instead of int _ID[2], reading the binary file
   * results in bogus data. */
  int _ID[2];
  uint64_t GetID() const {
      return *reinterpret_cast<const uint64_t*>(_ID);
  }
};

struct BonsaiGroupParticle_t {
  long long HostHaloID;
};

void create_TipsyHeader_MPI_type(MPI_Datatype &dtype);

class BonsaiSimReader_t
{
  vector<HBTInt> np_file;
  vector<HBTInt> offset_file;
  TipsyHeader_t Header;

  void ReadHeader(TipsyHeader_t &header);
  HBTInt CompileFileOffsets(int nfiles);

  /* All things snapshot particle information . */
  string SnapshotName;
  void OpenFile(std::ifstream &file);
  void GetSnapshotFileName(std::string &filename);
  void SetSnapshotFileName(int snapshotId);
  void ReadSnapshot(int ifile, Particle_t *ParticlesInFile, HBTInt file_start, HBTInt file_count);

  /* All things FoF groups. */
  string GroupFileName;
  void ReadGroupHeader(TipsyHeader_t &header);
  void OpenGroupFile(std::ifstream &file);
  void GetGroupFileName(std::string &filename);
  void SetGroupFileName(int snapshotId);
  void ReadGroupParticles(int ifile, Particle_t *ParticlesInFile, HBTInt file_start, HBTInt file_count);

  MPI_Datatype MPI_TipsyHeader_t;

public:
  BonsaiSimReader_t()
  {
    create_TipsyHeader_MPI_type(MPI_TipsyHeader_t);
  }
  ~BonsaiSimReader_t()
  {
    MPI_Type_free(&MPI_TipsyHeader_t);
  }
  void LoadSnapshot(MpiWorker_t &world, int snapshotId, vector<Particle_t> &Particles, Cosmology_t &Cosmology);
  void LoadGroups(MpiWorker_t &world, int snapshotId, vector<Halo_t> &Halos);

  /* To copy over the information of snapshot particle information */
  vector<Particle_t> ParticleHosts;
};

extern bool IsBonsaiSimGroup(const string &GroupFileFormat);
#endif
