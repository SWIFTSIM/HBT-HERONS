// TODO: unify the reference frame for specificProperties...
#include <algorithm>
#include <iostream>
#include <new>
#include <omp.h>

#include "datatypes.h"
#include "gravity_tree.h"
#include "snapshot_number.h"
#include "subhalo.h"

struct ParticleEnergy_t
{
  float Energy;
  HBTInt ParticleIndex;
};

inline bool CompareEnergy(const ParticleEnergy_t &a, const ParticleEnergy_t &b)
{
  return (a.Energy < b.Energy);
}

inline bool IsBound(const ParticleEnergy_t &a)
{
  return (a.Energy < 0);
}

inline bool IsNotSubsampleParticleType(const Particle_t &a)
{
  return a.DoNotSubsample();
}

/* Separates unbound particles from bound particles by placing them at the end
 * of Elist */
static HBTInt RemoveUnboundParticles(vector<ParticleEnergy_t> &Elist, const size_t NumPart)
{
  /* Separate bound from unbound particles. Stable partition since we want to 
   * keep original relative ordering (bound particles that cannot be subsampled
   * will stay at the beginning). */
  auto iter = std::stable_partition(Elist.begin(), Elist.begin() + NumPart, IsBound);

  /* Equals the number of bound particles*/
  return iter - Elist.begin();
}

class EnergySnapshot_t : public Snapshot_t
{
  HBTInt GetParticle(HBTInt i) const
  {
    return Elist[i].ParticleIndex;
  }

  HBTReal MassFactor;

public:
  ParticleEnergy_t *Elist;
  typedef vector<Particle_t> ParticleList_t;
  HBTInt N;
  const ParticleList_t &Particles;

  EnergySnapshot_t(ParticleEnergy_t *e, HBTInt n, const ParticleList_t &particles, const Snapshot_t &epoch)
    : Elist(e), N(n), Particles(particles), MassFactor(1.)
  {
    Cosmology = epoch.Cosmology;
  };

  HBTInt size() const
  {
    return N;
  }

  HBTInt GetId(HBTInt i) const
  {
    return Particles[GetParticle(i)].Id;
  }

  /* Sets how much more massive subsampled particles are. */
  void SetMassUpscaleFactor(HBTReal factor)
  {
    MassFactor = factor;
  }

  /* Returns the (scaled) mass of a particle. It will be larger than the true
   * mass if the particles are being subsampled. */
  HBTReal GetMass(HBTInt i) const
  {
    if(IsNotSubsampleParticleType(Particles[GetParticle(i)]))
      return Particles[GetParticle(i)].Mass;
    else
      return Particles[GetParticle(i)].Mass * MassFactor;
  }

  HBTReal GetInternalEnergy(HBTInt i) const
  {
#ifdef HAS_THERMAL_ENERGY
    return Particles[GetParticle(i)].InternalEnergy;
#else
    return 0.;
#endif
  }

  HBTReal GetPotentialEnergy(HBTInt i, const HBTxyz &refPos, const HBTxyz &refVel) const
  {
    // Load the total binding energy of the particle, then remove the thermal
    // and kinetic terms so we can return the potential energy
    HBTReal E = Elist[i].Energy;
#ifdef UNBIND_WITH_THERMAL_ENERGY
    E -= GetInternalEnergy(i);
#endif
    const HBTxyz &x = GetComovingPosition(i);
    const HBTxyz v = GetPhysicalVelocity(i);
    double dx[3], dv[3];
    for (int j = 0; j < 3; j++)
    {
      dx[j] = x[j] - refPos[j];
      if (HBTConfig.PeriodicBoundaryOn)
        dx[j] = NEAREST(dx[j]);
      dx[j] *= Cosmology.ScaleFactor; // physical
      dv[j] = v[j] - refVel[j] + Cosmology.Hz * dx[j];
      E -= 0.5 * dv[j] * dv[j];
    }
    return E;
  }

  const HBTxyz GetPhysicalVelocity(HBTInt i) const
  {
    return Particles[GetParticle(i)].GetPhysicalVelocity();
  }

  const HBTxyz &GetComovingPosition(HBTInt i) const
  {
    return Particles[GetParticle(i)].ComovingPosition;
  }

  /* Computes the total mass of the first NumPart most bound particles. */
  double SumUpMass(HBTInt NumPart)
  {
    double msum = 0;
#pragma omp parallel for reduction(+ : msum) if (NumPart > 100)
    for (HBTInt i = 0; i < NumPart; i++)
    {
      msum += GetMass(i);
    }
    return msum;
  }

  /* Mass-weighted average velocity */
  double AverageVelocity(HBTxyz &CoV, HBTInt NumPart)
  {
    HBTInt i;
    double svx = 0, svy = 0, svz = 0, msum = 0.;

#pragma omp parallel for reduction(+ : msum, svx, svy, svz) if (NumPart > 100)
    for (i = 0; i < NumPart; i++)
    {
      HBTReal m = GetMass(i);
      const HBTxyz v = GetPhysicalVelocity(i);
      msum += m;
      svx += v[0] * m;
      svy += v[1] * m;
      svz += v[2] * m;
    }

    CoV[0] = svx / msum;
    CoV[1] = svy / msum;
    CoV[2] = svz / msum;
    return msum;
  }

  /* Mass-weighted average position*/
  double AveragePosition(HBTxyz &CoM, HBTInt NumPart)
  {
    HBTInt i;
    double sx = 0, sy = 0, sz = 0, msum = 0;

    double origin[3];
    if (HBTConfig.PeriodicBoundaryOn)
      for (int j = 0; j < 3; j++)
        origin[j] = GetComovingPosition(0)[j];

#pragma omp parallel for reduction(+ : msum, sx, sy, sz) if (NumPart > 100)
    for (i = 0; i < NumPart; i++)
    {
      HBTReal m = GetMass(i);
      const HBTxyz &x = GetComovingPosition(i);
      msum += m;
      if (HBTConfig.PeriodicBoundaryOn)
      {
        sx += NEAREST(x[0] - origin[0]) * m;
        sy += NEAREST(x[1] - origin[1]) * m;
        sz += NEAREST(x[2] - origin[2]) * m;
      }
      else
      {
        sx += x[0] * m;
        sy += x[1] * m;
        sz += x[2] * m;
      }
    }
    sx /= msum;
    sy /= msum;
    sz /= msum;
    if (HBTConfig.PeriodicBoundaryOn)
    {
      sx += origin[0];
      sy += origin[1];
      sz += origin[2];
    }
    CoM[0] = sx;
    CoM[1] = sy;
    CoM[2] = sz;
    return msum;
  }

  void AverageKinematics(float &SpecificPotentialEnergy, float &SpecificKineticEnergy, float SpecificAngularMomentum[3],
                         HBTInt NumPart, const HBTxyz &refPos, const HBTxyz &refVel)
  /*obtain specific potential, kinetic energy, and angular momentum for the first NumPart particles
   * all quantities are physical

   * Note there is a slight inconsistency in the energy since they were calculated from the previous unbinding loop, but
   the refVel has been updated.
   */
  {
    if (NumPart <= 1)
    {
      SpecificPotentialEnergy = 0.;
      SpecificKineticEnergy = 0.;
      SpecificAngularMomentum[0] = SpecificAngularMomentum[1] = SpecificAngularMomentum[2] = 0.;
      return;
    }
    double E = 0., K = 0., AMx = 0., AMy = 0., AMz = 0., M = 0.;
#pragma omp parallel for reduction(+ : E, K, AMx, AMy, AMz, M) if (NumPart > 100)
    for (HBTInt i = 0; i < NumPart; i++)
    {
      HBTReal m = GetMass(i);
      E += Elist[i].Energy * m;
#ifdef UNBIND_WITH_THERMAL_ENERGY
      E -= GetInternalEnergy(i) * m;
#endif
      const HBTxyz &x = GetComovingPosition(i);
      const HBTxyz v = GetPhysicalVelocity(i);
      double dx[3], dv[3];
      for (int j = 0; j < 3; j++)
      {
        dx[j] = x[j] - refPos[j];
        if (HBTConfig.PeriodicBoundaryOn)
          dx[j] = NEAREST(dx[j]);
        dx[j] *= Cosmology.ScaleFactor; // physical
        dv[j] = v[j] - refVel[j] + Cosmology.Hz * dx[j];
        K += dv[j] * dv[j] * m;
      }
      AMx += (dx[1] * dv[2] - dx[2] * dv[1]) * m;
      AMy += (dx[2] * dv[0] - dx[0] * dv[2]) * m;
      AMz += (dx[0] * dv[1] - dx[1] * dv[0]) * m;
      M += m;
    }
    E /= M;
    K *= 0.5 / M;
    SpecificPotentialEnergy = E - K;
    SpecificKineticEnergy = K;
    SpecificAngularMomentum[0] = AMx / M;
    SpecificAngularMomentum[1] = AMy / M;
    SpecificAngularMomentum[2] = AMz / M;
  }
};

inline void RefineBindingEnergyOrder(EnergySnapshot_t &ESnap, HBTInt Size, GravityTree_t &tree, HBTxyz &RefPos,
                                     HBTxyz &RefVel)
{ // reorder the first Size particles according to their self-binding energy
  auto &Elist = ESnap.Elist;
  auto &Particles = ESnap.Particles;
  tree.Build(ESnap, Size);
  vector<ParticleEnergy_t> Einner(Size);
#pragma omp parallel if (Size > 100)
  {
#pragma omp for
    for (HBTInt i = 0; i < Size; i++)
    {
      HBTInt index = Elist[i].ParticleIndex;
      Einner[i].ParticleIndex = i;
      Einner[i].Energy = tree.BindingEnergy(Particles[index].ComovingPosition, Particles[index].GetPhysicalVelocity(), RefPos,
                                       RefVel, Particles[index].Mass);
    }
#pragma omp single
    sort(Einner.begin(), Einner.end(), CompareEnergy);
#pragma omp for
    for (HBTInt i = 0; i < Size; i++)
    {
      Einner[i] = Elist[Einner[i].ParticleIndex];
    }
#pragma omp for
    for (HBTInt i = 0; i < Size; i++)
    {
      Elist[i] = Einner[i];
    }
  }
}

/* Finds the factor by which the mass of particles need to be multiplied after
 * subsampling, to ensure mass conservation. */
HBTReal GetMassUpscaleFactor(const EnergySnapshot_t &ESnap, const HBTInt &Nlast, const HBTReal &Mlast, const HBTInt &MaxSampleSize, const HBTInt &Nunsample)
{
  /* Compute the total mass of particles that are not subsampled, to subtract
   * contribution from Mlast and hence get total true mass of subsampled 
   * particles. */
  HBTReal Munsampled = 0;
  if(Nunsample > 0)
  {
#pragma omp parallel for if (Nunsample > 100) reduction(+:Munsampled)
    for (HBTInt i = 0; i < Nunsample; i++)
    {
      Munsampled += ESnap.GetMass(i);
    }
  }

  /* This is the mass value that we need to convert when doing subsampling. */
  HBTReal MsubsampleTrue = Mlast - Munsampled;

  /* Total mass of the subsampled particle set. */
  HBTReal Msubsample = 0;
#pragma omp parallel for if (MaxSampleSize > 100) reduction(+:Msubsample)
  for (HBTInt i = Nunsample; i < (Nunsample + MaxSampleSize); i++)
  {
    Msubsample += ESnap.GetMass(i);
  }

  return MsubsampleTrue / Msubsample;
}

/* Randomly shuffles the particles whose type are eligible to be subsampled 
 * during unbinding. Particles types that are not eligible will be placed at the
 * start of the vector. */
HBTInt PrepareParticlesForSubsampling(vector<Particle_t> &Particles)
{
  /* We first partition the particle vector into those which we will 
   * subsample and those which we will not. */
  auto iter = std::partition(Particles.begin(), Particles.end(), IsNotSubsampleParticleType);

  /* Amount of particles that will not be subsampled. */
  HBTInt Nunsample = iter - Particles.begin();

  /* Shuffle the particles that can be subsampled. */
  std::random_shuffle(Particles.begin() + Nunsample, Particles.end());

  return Nunsample;
}

/* Counts how many bound particles are not eligible to be subsampled. */
HBTInt CountUnsampledParticles(const vector<ParticleEnergy_t> &Elist, const vector<Particle_t> &Particles, const HBTInt &OldNsubsample)
{
  HBTInt NewNunsample = 0;
  for (HBTInt i = 0; i < OldNsubsample; i++)
  {
    const auto &p = Particles[Elist[i].ParticleIndex];

    /* Particles that cannot be sampled are always at the beginning. Thus we can
     * exit the loop if we find a particle that we can subsample. */
    if(!p.DoNotSubsample())
      break;

    NewNunsample++;
  }
  return NewNunsample;
}

void Subhalo_t::Unbind(const Snapshot_t &epoch)
{ // the reference frame (pos and vel) should already be initialized before unbinding.

  /* We skip already existing orphans */
  if (!Particles.size())
  {
    Nbound = Particles.size();
    CountParticles();
    GetCorePhaseSpaceProperties();

    /* No bound particles, hence zero binding/potential energies will be saved */
    if (HBTConfig.SaveBoundParticleBindingEnergies)
      ParticleBindingEnergies.clear();
    if (HBTConfig.SaveBoundParticlePotentialEnergies)
      ParticlePotentialEnergies.clear();

    return;
  }

  /* We only expect (potentially) resolved subhaloes to make it here, or masked
   * out subhaloes (which should have at least one tracer particle if everything
   * is working correctly) */
  assert(Particles.size() >= 1);

  HBTInt MaxSampleSize = HBTConfig.MaxSampleSizeOfPotentialEstimate;
  bool RefineMostBoundParticle = (MaxSampleSize > 0 && HBTConfig.RefineMostBoundParticle);
  HBTReal BoundMassPrecision = HBTConfig.BoundMassPrecision;

  /* Need to initialise here, since orphans/disrupted objects do not call the
   * function used to set the value of TracerIndex (CountParticleTypes). This
   * prevents accessing entries beyond the corresponding particle array. */
  SetTracerIndex(0);

  /* Variables that store the centre of mass reference frame, which is used during
   * unbinding. In the first iteration, the centre of mass reference frame is that
   * of all particles associated to the subhalo in the previous output. Subsequent
   * iterations use the centre of mass frame of particles identified as bound. */
  HBTxyz OldRefPos, OldRefVel;
  auto &RefPos = ComovingAveragePosition;
  auto &RefVel = PhysicalAverageVelocity;

  /* Starting number of bound particles we include all particles because we do
   * not know which ones are bound or unbound at this point. */
  Nbound = Particles.size();

  GravityTree_t tree;
  tree.Reserve(Particles.size());

  /* Shuffle the particle vector, which we use as our basis of random 
   * subsamples. */
  HBTInt Nunsample = 0;
  if (MaxSampleSize > 0 && Nbound > MaxSampleSize)
    Nunsample = PrepareParticlesForSubsampling(Particles);

  /* This vector stores the original ordering of particles, and will later store
   * their binding energies. */
  vector<ParticleEnergy_t> Elist(Particles.size());
  for (HBTInt i = 0; i < Particles.size(); i++)
    Elist[i].ParticleIndex = i;

  EnergySnapshot_t ESnap(Elist.data(), Elist.size(), Particles, epoch);

  /* Associated mass of the starting number of particles. */
  Mbound = ESnap.SumUpMass(Nbound);

  /* Used to determine when iterative unbinding has converged to within the 
   * specified accuracy. */
  HBTInt Nlast;
  HBTReal Mlast;

  /* Iteratively unbind until we find that the subhalo bound particle number (or 
   * mass) has either converged or the subhalo has disrupted. */
  bool CorrectionLoop = false;
  while (true)
  {
    /* Correct the binding energies of bound particles due to removing particles
     * from the bound set in the previous iteration. */
    if (CorrectionLoop)
    {
      /* Compute the kinetic energy of the new centre of mass relative to the
       * old centre of mass frame. */
      HBTxyz RefVelDiff;
      epoch.RelativeVelocity(OldRefPos, OldRefVel, RefPos, RefVel, RefVelDiff);
      HBTReal dK = 0.5 * VecNorm(RefVelDiff);

      /* The tree will only contain particles that were identified as unbound in 
       * the last unbinding iteration. */
      EnergySnapshot_t ESnapCorrection(&Elist[Nbound], Nlast - Nbound, Particles,
                                       epoch);
      tree.Build(ESnapCorrection);

#pragma omp parallel for if (Nlast > 100)
      for (HBTInt i = 0; i < Nbound; i++)
      {
        HBTInt index = Elist[i].ParticleIndex;
        auto &x = Particles[index].ComovingPosition;
        auto v = Particles[index].GetPhysicalVelocity();
        HBTxyz OldVel;
        epoch.RelativeVelocity(x, v, OldRefPos, OldRefVel, OldVel);

        /* Remove the potential contribution of unbound particles and update the
         * kinetic energy based on new reference frame. */
        Elist[i].Energy += VecDot(OldVel, RefVelDiff) + dK - tree.EvaluatePotential(x, 0);
      }
      Nlast = Nbound;
    }
    else
    {
      Nlast = Nbound;
      Mlast = Mbound;
      HBTInt np_tree = Nlast;

      /* If we subsample, then we need to upscale the masses of particles that
       * will be subsampled when doing potential calculations. */
      if ((MaxSampleSize > 0) && (Nlast > (MaxSampleSize + Nunsample)))
      {
        np_tree = MaxSampleSize + Nunsample;

        ESnap.SetMassUpscaleFactor(1.); /* To get true particle mass */
        HBTReal MassUpscaleFactor = GetMassUpscaleFactor(ESnap, Nlast, Mlast, MaxSampleSize, Nunsample);
        ESnap.SetMassUpscaleFactor(MassUpscaleFactor);
      }

      tree.Build(ESnap, np_tree);
#pragma omp parallel for if (Nlast > 100)
      for (HBTInt i = 0; i < Nlast; i++)
      {
        /* Non-zero masses for particles in the tree because we need to remove
         * their self-gravity. */
        HBTReal particle_mass = (i < np_tree) ? ESnap.GetMass(i) : 0.;

        HBTInt index = Elist[i].ParticleIndex;
        Elist[i].Energy = tree.BindingEnergy(Particles[index].ComovingPosition, 
                                             Particles[index].GetPhysicalVelocity(),
                                             RefPos, RefVel, particle_mass);
#ifdef UNBIND_WITH_THERMAL_ENERGY
        Elist[i].Energy += Particles[index].InternalEnergy;
#endif
      }

      /* This ensures we use the true particle mass if there is no subsampling 
       * in the next unbinding iteration. */
      ESnap.SetMassUpscaleFactor(1.);
    }

    /* If we disable stripping, the we do no update Nbound nor Nunsample, and 
     * they retain the values assigned before the first unbinding iteration. */
#ifndef NO_STRIPPING
    Nbound = RemoveUnboundParticles(Elist, Nlast); // TODO: parallelize this.
    Nunsample = CountUnsampledParticles(Elist, Particles, Nunsample);
#endif

    // Count the number of bound tracer particles
#ifdef DM_ONLY
    // All particles are tracers in DMO runs
    HBTInt Nbound_tracers = Nbound;
#else
    HBTInt Nbound_tracers = 0;
    for (HBTInt i = 0; i < Nbound; i += 1)
    {
      const auto &p = Particles[Elist[i].ParticleIndex];
      if (p.IsTracer())
        Nbound_tracers += 1;
      if (Nbound_tracers >= HBTConfig.MinNumTracerPartOfSub)
        break; // We found enough, so no need to continue
    }
#endif

    /* Subhalo has disrupted */
    if ((Nbound < HBTConfig.MinNumPartOfSub) || (Nbound_tracers < HBTConfig.MinNumTracerPartOfSub))
    {
      /* Store when it disrupted. */
      if (IsAlive())
        SnapshotIndexOfDeath = epoch.GetSnapshotIndex();

      /* The most bound positions of the new orphan were found when updating
       * every subhalo particles. Copy over to the comoving ones. For future
       * outputs, we will rely on UpdateMostBoundPosition instead */
      copyHBTxyz(ComovingAveragePosition, ComovingMostBoundPosition);
      copyHBTxyz(PhysicalAverageVelocity, PhysicalMostBoundVelocity);

      /* Do not allow the orphan to have any particles, so they can be subject to
       * unbinding in their parent. The particle array will be updated after this
       * subhalo has been done, when truncating the source. */
      Nbound = 0;
      Mbound = 0;

      break;
    }

    /* Sort the particles that were found to be newly unbound. */
    sort(Elist.begin() + Nbound, Elist.begin() + Nlast, CompareEnergy);
    HBTInt Ndiff = Nlast - Nbound;
    if (Ndiff < Nbound)
    {
      if (MaxSampleSize <= 0 || Ndiff < MaxSampleSize)
      {
        CorrectionLoop = true;
        copyHBTxyz(OldRefPos, RefPos);
        copyHBTxyz(OldRefVel, RefVel);
      }
    }

    /* The centre of mass frame is updated here */
    Mbound = ESnap.AverageVelocity(PhysicalAverageVelocity, Nbound);
    ESnap.AveragePosition(ComovingAveragePosition, Nbound);

    /* We have converged according to the user specified threshold. */
    if (Nbound >= Nlast * BoundMassPrecision)
    {
      /* Since we have a resolved subhalo, we reset the death and sink 
       * information. */
      if (!IsAlive())
        SnapshotIndexOfDeath = SpecialConst::NullSnapshotId;
      if (IsTrapped())
      {
        SnapshotIndexOfSink = SpecialConst::NullSnapshotId;
        SinkTrackId = SpecialConst::NullTrackId;
      }

      /* Sort the bound particles in binding energy */
      sort(Elist.begin(), Elist.begin() + Nbound, CompareEnergy);

      /* We need to refine the most bound particle, as subsampling large subhaloes will lead to
       * incorrect ordering of binding energies. Hence, the most bound particle before this step
       * may not be the true most bound particle. */
      if (RefineMostBoundParticle && Nbound > MaxSampleSize)
      {
        /* If the number of bound particles is large, the number of particles used in this step scales with Nbound.
         * Using too few particles without this scaling would not result in a better centering. This is because it
         * would be limited to the (MaxSampleSize / Nbound) fraction of most bound particles, whose ranking can be
         * extremely sensitive to the randomness used during unbinding. */
        HBTInt SampleSizeCenterRefinement =
          max(MaxSampleSize, static_cast<HBTInt>(HBTConfig.BoundFractionCenterRefinement * Nbound));

        /* Computes self-binding energy of the SampleSizeCenterRefinement most bound
         * particles, and sorts them according to their binding energy. */
        // NOTE: Elist.Energies is not updated, so if we save binding energies
        // in the output, the values will be "out of order" since they reflect 
        // the subsampled energy estimate.
        RefineBindingEnergyOrder(ESnap, SampleSizeCenterRefinement, tree, RefPos, RefVel);
      }

      /* Replaces the original particle array with a copy where bound and unbound
       * particles are partitioned, and the bound particles are sorted in binding
       * energy. */
      // todo: optimize this with in-place permutation, to avoid mem alloc and copying.
      ParticleList_t p(Particles.size());
      for (HBTInt i = 0; i < Particles.size(); i++)
      {
        p[i] = Particles[Elist[i].ParticleIndex];
        Elist[i].ParticleIndex = i; // update particle index in Elist as well.
      }
      Particles.swap(p);

      /* Update the most bound coordinate. Note that for resolved subhaloes,
       * this is not necessarily a tracer particle. */
      copyHBTxyz(ComovingMostBoundPosition, Particles[0].ComovingPosition);
      copyHBTxyz(PhysicalMostBoundVelocity, Particles[0].GetPhysicalVelocity());
      break;
    }
  }

  /* Computes the specific potential and kinetic energy of the bound subhalo, as
   * well as its specific angular momentum. */
  ESnap.AverageKinematics(SpecificSelfPotentialEnergy, SpecificSelfKineticEnergy, SpecificAngularMomentum, Nbound,
                          RefPos, RefVel);

  /* For orphans, this function call only sets it MboundType and NboundType equal to 0. For resolved objects, it
   * updates those fields, as well as the index of the most bound tracer particle.*/
  CountParticleTypes();

  /* At this stage we know the updated TracerIndex, so if we are bound we should
   * update the most bound ID. */
  if (IsAlive())
    MostBoundParticleId = Particles[GetTracerIndex()].Id;

  /* Get centre of mass position and velocity of the most bound particles, which
   * is later used to determine if this subhalo overlaps in phase-space with another. */
  GetCorePhaseSpaceProperties();

  /* Store the binding energy information to save later */
  if (HBTConfig.SaveBoundParticleBindingEnergies)
  {
    ParticleBindingEnergies.resize(Nbound);
#pragma omp parallel for if (Nbound > 100)
    for (HBTInt i = 0; i < Nbound; i++)
      ParticleBindingEnergies[i] = Elist[i].Energy;
  }

  /* Store the potential energy information to save later */
  if (HBTConfig.SaveBoundParticlePotentialEnergies)
  {
    ParticlePotentialEnergies.resize(Nbound);
#pragma omp parallel for if (Nbound > 100)
    for (HBTInt i = 0; i < Nbound; i++)
      ParticlePotentialEnergies[i] = ESnap.GetPotentialEnergy(i, RefPos, RefVel);
  }
}

void Subhalo_t::RecursiveUnbind(SubhaloList_t &Subhalos, const Snapshot_t &snap)
{
  /* Unbind all subhaloes that are nested deeper in the hierarchy of the current
   * one. */
  for (HBTInt i = 0; i < NestedSubhalos.size(); i++)
  {
    /* One of the children of the current subhalo */
    auto subid = NestedSubhalos[i];
    auto &subhalo = Subhalos[subid];
    subhalo.RecursiveUnbind(Subhalos, snap);

    /* The unbound particles of the child we just subjected to unbinding are
     * accreted to the source of the current subhalo. */
    Particles.insert(Particles.end(), subhalo.Particles.begin() + subhalo.Nbound, subhalo.Particles.end());
  }

  /* Unbind the current subhalo */
  Unbind(snap);

  /* Check if any of the subgroups deeper in this tree's hierarchy merge with it.
   * We update the particle list and the entries of the merged subhaloes if the
   * option is enabled. */
  bool HasExperiencedMerger = MergeRecursively(Subhalos, snap, *this);

  /* We need to subject the subhalo to unbinding once more, as it has accreted
   * new particles as a result of mergers. */
  if (HBTConfig.MergeTrappedSubhalos && HasExperiencedMerger)
    Unbind(snap);

  /* We are now sure about which particles are bound to this subhalo, so we can
   * safely pass the unbound ones to its parent and truncate the source.*/
}

void Subhalo_t::TruncateSource()
{
  HBTInt Nsource;
  if (Nbound <= 1)
    Nsource = Nbound;
  else
    Nsource = Nbound * HBTConfig.SourceSubRelaxFactor;
  if (Nsource > Particles.size())
    Nsource = Particles.size();
  Particles.resize(Nsource);
}

void SubhaloSnapshot_t::RefineParticles()
{ // it's more expensive to build an exclusive list. so do inclusive here.
  // TODO: ensure the inclusive unbinding is stable (contaminating particles from big subhaloes may hurdle the unbinding

#ifdef INCLUSIVE_MASS
#pragma omp parallel for schedule(dynamic, 1) if (ParallelizeHaloes)
  for (HBTInt subid = 0; subid < Subhalos.size(); subid++)
  {
    Subhalos[subid].Unbind(*this);
    Subhalos[subid].TruncateSource();
  }
#else
  HBTInt NumHalos = MemberTable.SubGroups.size();
#pragma omp parallel for schedule(dynamic, 1) if (ParallelizeHaloes)
  for (HBTInt haloid = 0; haloid < NumHalos; haloid++)
  {
    auto &subgroup = MemberTable.SubGroups[haloid];
    if (subgroup.size() == 0)
      continue;
    // add new satellites to central's NestedSubhalos
    auto &central = Subhalos[subgroup[0]];
    auto &nests = central.NestedSubhalos;
    auto old_membercount = nests.size();
    auto &heads = MemberTable.SubGroupsOfHeads[haloid];
    // update central member list (append other heads except itself)
    nests.insert(nests.end(), heads.begin() + 1, heads.end());
    central.RecursiveUnbind(Subhalos, *this);
    nests.resize(old_membercount); // restore old satellite list
  }
// unbind field subs
#pragma omp parallel
  {
    HBTInt NumField = MemberTable.SubGroups[-1].size();
#pragma omp for schedule(dynamic, 1) nowait
    for (HBTInt i = 0; i < NumField; i++)
    {
      HBTInt subid = MemberTable.SubGroups[-1][i];
      Subhalos[subid].Unbind(*this);
    }
    // unbind new-born subs
    HBTInt NumSubOld = MemberTable.AllMembers.size(), NumSub = Subhalos.size();
#pragma omp for schedule(dynamic, 1)
    for (HBTInt i = NumSubOld; i < NumSub; i++)
    {
      Subhalos[i].Unbind(*this);
    }
#pragma omp for schedule(dynamic, 1)
    for (HBTInt i = 0; i < NumSub; i++)
      Subhalos[i].TruncateSource();
  }
#endif
}