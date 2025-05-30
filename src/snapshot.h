#ifndef SNAPSHOT_H_INCLUDED
#define SNAPSHOT_H_INCLUDED

#include <assert.h>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <unordered_map>

#include "config_parser.h"
#include "datatypes.h"
#include "hash.h"
#include "mpi_wrapper.h"
#include "mymath.h"
#include "snapshot_number.h"

struct Cosmology_t
{
  HBTReal OmegaM0;
  HBTReal OmegaLambda0;
  HBTReal ScaleFactor;

  // current Hubble param in internal units
  // we try to read this in if possible, otherwise derive it
  HBTReal Hz = -1.0;

  void Set_Hz(double hz)
  {
    Hz = hz;
  }

  void Set(double scalefactor, double omega0, double omegaLambda0)
  {
    OmegaM0 = omega0;
    OmegaLambda0 = omegaLambda0;
    ScaleFactor = scalefactor;

    if (Hz == -1.0)
    {
      double Hratio;
      Hratio = sqrt(omega0 / (scalefactor * scalefactor * scalefactor) +
                         (1 - omega0 - omegaLambda0) / (scalefactor * scalefactor) +
                         omegaLambda0);
      Hz = Hratio * PhysicalConst::H0;
    }
  }
};

struct RadMassVel_t
{
  HBTReal r, m, v;
  RadMassVel_t(){};
  RadMassVel_t(HBTReal r, HBTReal m) : r(r), m(m)
  {
  }
  RadMassVel_t(HBTReal r, HBTReal m, HBTReal v) : r(r), m(m), v(v)
  {
  }
};

struct Particle_t
{
  HBTInt Id;
  HBTxyz ComovingPosition;
  HBTvel PhysicalVelocity;
  HBTMassType Mass;
#ifndef DM_ONLY
#ifdef HAS_THERMAL_ENERGY
  HBTReal InternalEnergy;
#endif
  ParticleType_t Type;
#endif
  HBTInt HostId;
  void create_MPI_type(MPI_Datatype &dtype);
  Particle_t(){};
  Particle_t(HBTInt id) : Id(id)
  {
  }
  bool operator==(const Particle_t &other) const
  {
    return Id == other.Id;
  }
  HBTxyz GetPhysicalVelocity() const
  {
    HBTxyz vel;
    copyXYZ(vel, PhysicalVelocity);
    return vel;
  }
  int IsTracer() const
  {
#ifdef DM_ONLY
    return 1;
#else
    return (1 << Type) & HBTConfig.TracerParticleBitMask;
#endif
  }
  int DoNotSubsample() const
  {
#ifdef DM_ONLY
    return 0;
#else
    return (1 << Type) & HBTConfig.DoNotSubsampleParticleBitMask;
#endif
  }
};
extern ostream &operator<<(ostream &o, Particle_t &p);

class Snapshot_t : public SnapshotNumber_t
{
public:
  Cosmology_t Cosmology;
  //   Snapshot_t()=default;
  virtual HBTInt size() const = 0;
  virtual HBTInt GetId(const HBTInt index) const
  {
    return index;
  }
  virtual const HBTxyz &GetComovingPosition(const HBTInt index) const = 0;
  virtual const HBTxyz GetPhysicalVelocity(const HBTInt index) const = 0;
  virtual HBTReal GetMass(const HBTInt index) const = 0;
  virtual HBTReal GetInternalEnergy(HBTInt index) const
  {
    return 0.;
  }
  void SphericalOverdensitySize(float &Mvir, float &Rvir, HBTReal VirialFactor, const vector<HBTReal> &RSorted,
                                HBTReal ParticleMass) const;
  void SphericalOverdensitySize(float &Mvir, float &Rvir, HBTReal VirialFactor, const vector<RadMassVel_t> &prof) const;
  void SphericalOverdensitySize2(float &Mvir, float &Rvir, HBTReal VirialFactor, const vector<HBTReal> &RSorted,
                                 HBTReal ParticleMass) const;
  void RelativeVelocity(const HBTxyz &targetPos, const HBTxyz &targetVel, const HBTxyz &refPos, const HBTxyz &refVel,
                        HBTxyz &relativeVel) const;
};

inline void Snapshot_t::RelativeVelocity(const HBTxyz &targetPos, const HBTxyz &targetVel, const HBTxyz &refPos,
                                         const HBTxyz &refVel, HBTxyz &relativeVel) const
{
  HBTxyz dx;
  HBTxyz &dv = relativeVel;
  for (int j = 0; j < 3; j++)
  {
    dx[j] = targetPos[j] - refPos[j];
    if (HBTConfig.PeriodicBoundaryOn)
      dx[j] = NEAREST(dx[j]);
    dv[j] = targetVel[j] - refVel[j];
    dv[j] += Cosmology.Hz * Cosmology.ScaleFactor * dx[j];
  }
}

class SnapshotView_t : public Snapshot_t
{
public:
  HBTInt *Ids;
  HBTInt N;
  Snapshot_t &Snapshot;
  SnapshotView_t(vector<HBTInt> &ids, Snapshot_t &fullsnapshot)
    : Ids(ids.data()), N(ids.size()), Snapshot(fullsnapshot), Snapshot_t(fullsnapshot){};
  SnapshotView_t(VectorView_t<HBTInt> &ids, Snapshot_t &fullsnapshot)
    : Ids(ids.data()), N(ids.size()), Snapshot(fullsnapshot), Snapshot_t(fullsnapshot){};
  SnapshotView_t(HBTInt *ids, HBTInt n, Snapshot_t &fullsnapshot)
    : Ids(ids), N(n), Snapshot(fullsnapshot), Snapshot_t(fullsnapshot){};
  void ReSize(HBTInt n)
  {
    N = n;
  }
  HBTInt size() const
  {
    return N;
  }
  HBTInt GetId(HBTInt i) const
  {
    return Snapshot.GetId(Ids[i]);
  }
  HBTReal GetMass(HBTInt i) const
  {
    return Snapshot.GetMass(Ids[i]);
  }
  const HBTxyz GetPhysicalVelocity(HBTInt i) const
  {
    return Snapshot.GetPhysicalVelocity(Ids[i]);
  }
  const HBTxyz &GetComovingPosition(HBTInt i) const
  {
    return Snapshot.GetComovingPosition(Ids[i]);
  }
};

class ParticleSnapshot_t : public Snapshot_t
{
  typedef vector<HBTInt> IndexList_t;

  FlatIndexTable_t<HBTInt, HBTInt> FlatHash;
  MappedIndexTable_t<HBTInt, HBTInt> MappedHash;
  IndexTable_t<HBTInt, HBTInt> *ParticleHash;

  void ExchangeParticles(MpiWorker_t &world);
  vector<HBTInt> PartitionParticles(MpiWorker_t &world);

public:
  vector<Particle_t> Particles;
  HBTInt NumberOfParticlesOnAllNodes;

  ParticleSnapshot_t()
    : Snapshot_t(), Particles(), ParticleHash(), MappedHash(), FlatHash(), NumberOfParticlesOnAllNodes(0)
  {
    if (HBTConfig.ParticleIdNeedHash)
      ParticleHash = &MappedHash;
    else
      ParticleHash = &FlatHash;
  }
  ParticleSnapshot_t(MpiWorker_t &world, int snapshot_index, bool fill_particle_hash = true) : ParticleSnapshot_t()
  {
    Load(world, snapshot_index, fill_particle_hash);
  }
  ~ParticleSnapshot_t()
  {
    Clear(); // not necessary
  }
  void FillParticleHash();
  void ClearParticleHash();
  void ClearParticles();

  HBTInt size() const;
  HBTInt GetId(HBTInt index) const;
  HBTInt GetIndex(HBTInt particle_id) const;
  HBTInt GetIndex(Particle_t &particle) const;
  template <class ParticleIdList_t>
  void GetIndices(ParticleIdList_t &particles) const;
  const HBTxyz &GetComovingPosition(HBTInt index) const;
  const HBTxyz GetPhysicalVelocity(HBTInt index) const;
  HBTReal GetMass(HBTInt index) const;
  HBTReal GetInternalEnergy(HBTInt index) const;
  ParticleType_t GetParticleType(HBTInt index) const;

  void Load(MpiWorker_t &world, int snapshot_index, bool fill_particle_hash = true);
  void Clear();

  /* For use when including particles that split. */
#ifndef DM_ONLY
  std::unordered_map<HBTInt, HBTInt> ParticleSplitMap;
#endif

  void AveragePosition(HBTxyz &CoM, const HBTInt Particles[], HBTInt NumPart) const;
  void AverageVelocity(HBTxyz &CoV, const HBTInt Particles[], HBTInt NumPart) const;

  template <class Halo_T>
  void ExchangeHalos(MpiWorker_t &world, vector<Halo_T> &InHalos, vector<Halo_T> &OutHalos,
                     MPI_Datatype MPI_Halo_Shell_Type) const;
};
inline HBTInt ParticleSnapshot_t::size() const
{
  return Particles.size();
}
inline HBTInt ParticleSnapshot_t::GetId(HBTInt index) const
{
  return Particles[index].Id;
}
inline HBTInt ParticleSnapshot_t::GetIndex(HBTInt particle_id) const
{
  return ParticleHash->GetIndex(particle_id);
}
inline HBTInt ParticleSnapshot_t::GetIndex(Particle_t &particle) const
{
  return ParticleHash->GetIndex(particle.Id);
}
inline const HBTxyz &ParticleSnapshot_t::GetComovingPosition(HBTInt index) const
{
  return Particles[index].ComovingPosition;
}
inline const HBTxyz ParticleSnapshot_t::GetPhysicalVelocity(HBTInt index) const
{
  return Particles[index].GetPhysicalVelocity();
}
inline HBTReal ParticleSnapshot_t::GetMass(HBTInt index) const
{
  return Particles[index].Mass;
}
inline HBTReal ParticleSnapshot_t::GetInternalEnergy(HBTInt index) const
{
#if !defined(DM_ONLY) && defined(HAS_THERMAL_ENERGY)
  return Particles[index].InternalEnergy;
#else
  return 0.;
#endif
}
inline ParticleType_t ParticleSnapshot_t::GetParticleType(HBTInt index) const
{
#ifdef DM_ONLY
  return TypeDM;
#else
  return Particles[index].Type;
#endif
}

extern double AveragePosition(HBTxyz &CoM, const Particle_t Particles[], HBTInt NumPart);
extern double AverageVelocity(HBTxyz &CoV, const Particle_t Particles[], HBTInt NumPart);
#endif
