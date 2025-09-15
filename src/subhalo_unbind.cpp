#include <new>
#include <omp.h>
#include <chrono>
#include <iostream>
#include <algorithm>

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

  /* Index in the Particles vector of the ith most bound particle. */
  HBTInt GetParticle(HBTInt i) const
  {
    return Elist[i].ParticleIndex;
  }

  vector<HBTReal> MassUpscaleFactor;

public:
  ParticleEnergy_t *Elist;
  typedef vector<Particle_t> ParticleList_t;
  HBTInt N;
  const ParticleList_t &Particles;

  EnergySnapshot_t(ParticleEnergy_t *e, HBTInt n, const ParticleList_t &particles, const Snapshot_t &epoch)
    : Elist(e), N(n), Particles(particles), MassUpscaleFactor(TypeMax, 1.)
  {
    Cosmology = epoch.Cosmology;
  };

  /* Used within the gravitational tree construction. */
  HBTInt size() const
  {
    return N;
  }

  /* Particle ID of the ith most bound particle. */
  HBTInt GetId(HBTInt i) const
  {
    return Particles[GetParticle(i)].Id;
  }

#ifndef DM_ONLY
  /* Particle type of the ith most bound particle. */
  HBTInt GetType(HBTInt i) const
  {
    return Particles[GetParticle(i)].Type;
  }
#endif

  /* Scaled mass of the ith most bound particle. It will be larger than the true
   * mass if the particles are being subsampled. */
  HBTReal GetMass(HBTInt i) const
  {
#ifdef DM_ONLY
    return Particles[GetParticle(i)].Mass * MassUpscaleFactor[1];
#else
    if(IsNotSubsampleParticleType(Particles[GetParticle(i)]))
      return Particles[GetParticle(i)].Mass;
    else
      return Particles[GetParticle(i)].Mass * MassUpscaleFactor[GetType(i)];
#endif
  }

  /* Internal energy of the ith most bound particle. */
  HBTReal GetInternalEnergy(HBTInt i) const
  {
#ifdef HAS_THERMAL_ENERGY
    return Particles[GetParticle(i)].InternalEnergy;
#else
    return 0.;
#endif
  }

  /* Potential energy of the ith most bound particle. */
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

  /* Physical velocity of the ith most bound particle. */
  const HBTxyz GetPhysicalVelocity(HBTInt i) const
  {
    return Particles[GetParticle(i)].GetPhysicalVelocity();
  }

  /* Comoving position of the ith most bound particle. */
  const HBTxyz &GetComovingPosition(HBTInt i) const
  {
    return Particles[GetParticle(i)].ComovingPosition;
  }

  /* Updates the location and velocity of the centre of mass of the NumPart most
   * bound particles. */
  double CentreOfMassReferenceFrame(HBTxyz &CentreOfMassPosition, HBTxyz &CentreOfMassVelocity, HBTInt NumPart)
  {
    double TotalMass = 0, WeightedPosition[3] = {0,0,0}, WeightedVelocity[3] = {0,0,0};

    /* Need a reference point to centre if we have periodic boundary conditions */
    double ReferencePosition[3] = {0,0,0};
    if (HBTConfig.PeriodicBoundaryOn)
      for(int dimension = 0; dimension < 3; dimension++)
        ReferencePosition[dimension] = GetComovingPosition(0)[dimension];

#pragma omp parallel for reduction(+ : TotalMass, WeightedPosition[:3], WeightedVelocity[:3]) if (NumPart > 100)
    for (HBTInt i = 0; i < NumPart; i++)
    {
      const HBTReal mass = GetMass(i);
      const HBTxyz position = GetComovingPosition(i);
      const HBTxyz velocity = GetPhysicalVelocity(i);

      TotalMass += mass;

      for(int dimension = 0; dimension < 3; dimension++)
      {
        WeightedPosition[dimension] += HBTConfig.PeriodicBoundaryOn ? \
                                       NEAREST(position[dimension] - ReferencePosition[dimension]) * mass \
                                     : position[dimension] * mass;
        WeightedVelocity[dimension] += velocity[dimension] * mass;
      }
    }

    /* Remove mass factor */
    for(int dimension = 0;  dimension < 3; dimension++)
    {
      WeightedPosition[dimension] /= TotalMass;
      WeightedVelocity[dimension] /= TotalMass;
    }

    /* Express centre of mass in box coordinates. */
    if (HBTConfig.PeriodicBoundaryOn)
      for(int dimension = 0; dimension < 3; dimension++)
        WeightedPosition[dimension] += ReferencePosition[dimension];

    for(int dimension = 0;  dimension < 3; dimension++)
    {
      CentreOfMassPosition[dimension] = WeightedPosition[dimension];
      CentreOfMassVelocity[dimension] = WeightedVelocity[dimension];
    }

    return TotalMass;
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

  /* Accumulates the total mass of particles between the Nstart and Nfinish most
   * bound particle. */
   void AccumulateMass(const HBTInt &Nstart, const HBTInt &Nfinish, float &Mass)
   {
      Mass = 0;
#pragma omp parallel for reduction(+:Mass) if ((Nfinish - Nstart) > 1000)
      for (HBTInt i = Nstart; i < Nfinish; i++)
        Mass += GetMass(i);
   }

#ifndef DM_ONLY
  /* Accumulates the total mass per particle type between the Nstart and Nfinish
   * most bound particles. */
   void AccumulateMass(const HBTInt &Nstart, const HBTInt &Nfinish, float MassPerType[])
   {
      for(int i = 0; i < TypeMax; i++)
        MassPerType[i] = 0;

#pragma omp parallel for reduction(+:MassPerType[:TypeMax]) if ((Nfinish - Nstart) > 1000)
      for (HBTInt i = Nstart; i < Nfinish; i++)
        MassPerType[GetType(i)] += GetMass(i);
   }

  /* Accumulates the total mass of particles and per particle type between the
   * Nstart and Nfinish most bound particle. */
   void AccumulateMass(const HBTInt &Nstart, const HBTInt &Nfinish, float &Mass, float MassPerType[])
   {
      /* Compute mass per type */
      AccumulateMass(Nstart, Nfinish, MassPerType);

      /* Sum over all types to get total */
      Mass = 0;
      for(int i = 0; i < TypeMax; i++)
        Mass += MassPerType[i];
   }
#endif

  /* Finds the factor by which the mass of particles need to be multiplied after
   * subsampling, to ensure mass conservation. */
  std::vector<HBTReal> GetMassUpscaleFactor(const HBTInt &Nlast, const HBTReal &Mlast, const HBTInt &MaxSampleSize, const HBTInt &Nunsample)
  {
    ResetMassUpscaleFactor(); /* To get true particle mass */

    /* Compute the total mass of particles that are not subsampled, to subtract
     * contribution from Mlast and hence get total true mass of subsampled
     * particles. */
    HBTReal Munsampled = 0;
    if(Nunsample > 0)
    {
#pragma omp parallel for if (Nunsample > 100) reduction(+:Munsampled)
      for (HBTInt i = 0; i < Nunsample; i++)
        Munsampled += GetMass(i);
    }

    /* This is the mass value that we need to convert when doing subsampling. */
    HBTReal MsubsampleTrue = Mlast - Munsampled;

    /* Total mass of the subsampled particle set. */
    HBTReal Msubsample = 0;
#pragma omp parallel for if (MaxSampleSize > 100) reduction(+:Msubsample)
    for (HBTInt i = Nunsample; i < (Nunsample + MaxSampleSize); i++)
    {
      Msubsample += GetMass(i);
    }

    return std::vector<HBTReal> (TypeMax, MsubsampleTrue / Msubsample);
  }

#ifndef DM_ONLY
  /* Finds the factor by which the mass of particles need to be multiplied after
   * subsampling, to ensure mass conservation. Contrary to GetMassUpscaleFactor,
   * we define the upscale factor based on a particle type level. */
  std::vector<HBTReal> GetMassUpscaleFactorPerParticleType(const HBTInt &Nlast, const HBTReal &Mlast, const HBTInt &MaxSampleSize, const HBTInt &Nunsample)
  {
    ResetMassUpscaleFactor(); /* To get true particle mass */

    /* Mass of particles in the subsampled set */
    float MassPerTypeSubsample[TypeMax];
    AccumulateMass(Nunsample, Nunsample + MaxSampleSize, MassPerTypeSubsample);

    /* Compute mass of particles that are not in the subsampled set. */
    float MassPerTypeTotal[TypeMax];
    AccumulateMass(Nunsample + MaxSampleSize, Nlast, MassPerTypeTotal);

    /* Add both together to get **actual** total mass */
    for(int type = 0; type < TypeMax; type++)
      MassPerTypeTotal[type] += MassPerTypeSubsample[type];

    /* Ratio of both gives the mass upscale factor for each particle type.
     * We do not set values with zero particle types, but those particle types
     * would not contribute to potential since by definition there are none in
     * the subsampled set. */
    std::vector<HBTReal> MassUpscaleFactor = std::vector<HBTReal>(TypeMax, 0);
    for(int type = 0; type < TypeMax; type++)
      if(MassPerTypeSubsample[type])
        MassUpscaleFactor[type] = MassPerTypeTotal[type] / MassPerTypeSubsample[type];

    /* Ensure total mass conservation if we are missing a particle type from
     * our subsampled set. */
    float TotalScaledMass = 0, TotalTrueMass = 0;
    for(int type = 0; type < TypeMax; type++)
    {
      TotalTrueMass += MassPerTypeTotal[type];
      TotalScaledMass += MassPerTypeSubsample[type] * MassUpscaleFactor[type];
    }

    for(int type = 0; type < TypeMax; type++)
      MassUpscaleFactor[type] *= TotalTrueMass / TotalScaledMass;

    return MassUpscaleFactor;
  }
#endif

  /* Makes it so GetMass(i) returns the true mass of the ith most bound
   * particle. */
  void ResetMassUpscaleFactor()
  {
    std::fill(MassUpscaleFactor.begin(), MassUpscaleFactor.end(), 1.);
  }

  /* Upscales the masses of particles if they are subsampled for potential
   * calculation purposes, or leaves them as is if the subhalo is small.
   * Returns the number of particles that will be gravity sources. */
  HBTInt SetMassUpscaleFactor(const HBTInt &Nlast, const HBTReal &Mlast, const HBTInt &MaxSampleSize, const HBTInt &Nunsample)
  {
    /* Subsampling disabled by user or the subhalo does not cross the threshold
     * to subsample. Use all currently bound particles in tree. */
    if((MaxSampleSize == 0) || (Nlast < (MaxSampleSize + Nunsample)))
      return Nlast;

#ifdef DM_ONLY
    MassUpscaleFactor = GetMassUpscaleFactor(Nlast, Mlast, MaxSampleSize, Nunsample);
#else
    if(HBTConfig.PotentialEstimateUpscaleMassesPerType)
      MassUpscaleFactor = GetMassUpscaleFactorPerParticleType(Nlast, Mlast, MaxSampleSize, Nunsample);
    else
      MassUpscaleFactor = GetMassUpscaleFactor(Nlast, Mlast, MaxSampleSize, Nunsample);
#endif

    return MaxSampleSize + Nunsample;
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

#ifdef MEASURE_UNBINDING_TIME
  /* We use default -1 for times because not all times will be valid, e.g.
   * for orphans there is no CentreRefinement. */
  NumberUnbindingIterations = 0;
  StartSubhalo = -1;
  StartUnbinding = -1;
  StartCentreRefinement = -1;
  StartPhaseSpace = -1;
  EndSubhalo = -1;
  Timer_t SubhaloTimer;
  SubhaloTimer.Tick("Start");
#endif

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

#ifdef MEASURE_UNBINDING_TIME
    SubhaloTimer.Tick("PhaseSpace");
    for (int i = 1; i < SubhaloTimer.Size(); i++)
    {
      double time = std::chrono::duration<double>(SubhaloTimer.tickers[i].time_since_epoch()).count();
      // SubhaloTimer.tickers[i].time_since_epoch()).count();
      // double time = 0;
      // chrono::duration_cast<chrono::duration<double>>(SubhaloTimer.tickers[i].time_since_epoch().count());

      if (SubhaloTimer.names[i] == "Start")
        StartSubhalo = time;
      if (SubhaloTimer.names[i] == "Unbinding")
        StartUnbinding = time;
      if (SubhaloTimer.names[i] == "CentreRefinement")
        StartCentreRefinement = time;
      if (SubhaloTimer.names[i] == "PhaseSpace")
        StartPhaseSpace = time;
      if (SubhaloTimer.names[i] == "End")
        EndSubhalo = time;
    }
#endif

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

  /* To access particles in binding energy order and most methods related to
   * unbinding */
  EnergySnapshot_t ESnap(Elist.data(), Elist.size(), Particles, epoch);

  /* Associated mass of the starting number of particles. */
#ifdef DM_ONLY
  ESnap.AccumulateMass(0, Nbound, Mbound);
#else
  ESnap.AccumulateMass(0, Nbound, Mbound, MboundType);
#endif

  /* Used to determine when iterative unbinding has converged to within the
   * specified accuracy. */
  HBTInt Nlast;
  HBTReal Mlast;

  /* Iteratively unbind until we find that the subhalo bound particle number (or
   * mass) has either converged or the subhalo has disrupted. */
  bool CorrectionLoop = false;
  while (true)
  {
#ifdef MEASURE_UNBINDING_TIME
    NumberUnbindingIterations++;
#endif

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

      /* Checks whether subhalo will be subsampled, and if so, scales the masses
       * of subsampled particles to conserve mass. */
      HBTInt NGravitySources = ESnap.SetMassUpscaleFactor(Nlast, Mlast, MaxSampleSize, Nunsample);
      tree.Build(ESnap, NGravitySources);

#pragma omp parallel for if (Nlast > 100)
      for (HBTInt i = 0; i < Nlast; i++)
      {
        /* Non-zero masses for particles in the tree because we need to remove
         * their self-gravity. */
        HBTReal particle_mass = (i < NGravitySources) ? ESnap.GetMass(i) : 0.;

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
      ESnap.ResetMassUpscaleFactor();
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
#ifdef MEASURE_UNBINDING_TIME
      SubhaloTimer.Tick("Unbinding");
#endif

      /* Store when it disrupted. */
      if (IsAlive())
        SnapshotOfDeath = epoch.GetSnapshotId();

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

    /* The centre of mass frame is updated here. We return bound mass because we
     * already need to calculate it, so we save another loop. */
    Mbound = ESnap.CentreOfMassReferenceFrame(ComovingAveragePosition, PhysicalAverageVelocity, Nbound);

    /* We have converged according to the user specified threshold. */
    if (Nbound >= Nlast * BoundMassPrecision)
    {
#ifdef MEASURE_UNBINDING_TIME
      SubhaloTimer.Tick("Unbinding");
#endif

      /* Since we have a resolved subhalo, we reset the death and sink
       * information. */
      if (!IsAlive())
      {
        SnapshotOfDeath = SpecialConst::NullSnapshotId;
        DescendantTrackId = SpecialConst::NullTrackId;
      }
      if (IsTrapped())
      {
        SnapshotOfSink = SpecialConst::NullSnapshotId;
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

#ifdef MEASURE_UNBINDING_TIME
      SubhaloTimer.Tick("CentreRefinement");
#endif

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

#ifdef MEASURE_UNBINDING_TIME
  SubhaloTimer.Tick("PhaseSpace");
  for (int i = 1; i < SubhaloTimer.Size(); i++)
  {
    // double time = chrono::duration_cast<chrono::duration<double>>(SubhaloTimer.tickers[i]);
    double time = std::chrono::duration<double>(SubhaloTimer.tickers[i].time_since_epoch()).count();

    if (SubhaloTimer.names[i] == "Start")
      StartSubhalo = time;
    if (SubhaloTimer.names[i] == "Unbinding")
      StartUnbinding = time;
    if (SubhaloTimer.names[i] == "CentreRefinement")
      StartCentreRefinement = time;
    if (SubhaloTimer.names[i] == "PhaseSpace")
      StartPhaseSpace = time;
  }
#endif // MEASURE_UNBINDING_TIME

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

void SubhaloSnapshot_t::RefineParticles(MpiWorker_t &world)
{ // it's more expensive to build an exclusive list. so do inclusive here.
  // TODO: ensure the inclusive unbinding is stable (contaminating particles from big subhaloes may hurdle the unbinding

  Timer_t ImbalanceTimer;
  ImbalanceTimer.Tick("start");

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
#endif // INCLUSIVE_MASS

  ImbalanceTimer.Tick("end");

  PrintTimeImbalanceStatistics(world, ImbalanceTimer);
  PrintSubhaloStatistics(world);
}

/* Prints the maximum time taken by a rank, and its associated imbalance. */
void SubhaloSnapshot_t::PrintTimeImbalanceStatistics(MpiWorker_t &world, Timer_t Timer)
{
  double LocalElapsedTime = Timer.GetSeconds();
  std::vector<double> AllElapsedTime(world.size(), 0);
  MPI_Allgather(&LocalElapsedTime, 1, MPI_DOUBLE, AllElapsedTime.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD);

  /* Compute average run time and maximum time across all ranks */
  double AverageTime = 0, MaxTime=0;
  for(auto &TimePerRank:AllElapsedTime)
  {
    AverageTime += TimePerRank;
    MaxTime = std::max(TimePerRank, MaxTime);
  }
  AverageTime /= world.size();

  if (world.rank() == 0)
    std::cout << "Took " << MaxTime << " seconds. Maximum imbalance across ranks was " << MaxTime / AverageTime << "." << std::endl;
}

/* Print information about how many subhaloes have been sunk, disrupted, newly
 * idenfied and failed subhaloes. */
void SubhaloSnapshot_t::PrintSubhaloStatistics(MpiWorker_t &world)
{
   HBTInt LocalSunkSubhaloes = 0, LocalDisruptedSubhaloes = 0, LocalNewSubhaloes = 0, LocalFakeSubhaloes = 0;
#pragma omp parallel for reduction(+ : LocalSunkSubhaloes, LocalDisruptedSubhaloes, LocalNewSubhaloes, LocalFakeSubhaloes) if (Subhalos.size() > 1000)
   for(size_t subhalo_index = 0;  subhalo_index < Subhalos.size(); subhalo_index++)
   {
      Subhalo_t &sub = Subhalos[subhalo_index];

      /* Sunk subhaloes */
      LocalSunkSubhaloes += (sub.SnapshotOfDeath == GetSnapshotId()) \
                          & (sub.SnapshotOfSink  == GetSnapshotId());

      /* Disrupted subhaloes */
      LocalDisruptedSubhaloes += (sub.SnapshotOfBirth < GetSnapshotId())  \
                               & (sub.SnapshotOfDeath == GetSnapshotId()) \
                               & (sub.SnapshotOfSink  == SpecialConst::NullSnapshotId);

      /* A subhalo that was never self-bound. */
      LocalFakeSubhaloes += (sub.SnapshotOfBirth == GetSnapshotId()) \
                          & (sub.SnapshotOfDeath == GetSnapshotId());

      /* A subhalo that was just identified as self-bound for the first time. */
      LocalNewSubhaloes += (sub.SnapshotOfBirth == GetSnapshotId()) \
                         & (sub.SnapshotOfDeath == SpecialConst::NullSnapshotId);
   }

  /* Gather across ranks */
  HBTInt TotalSunkSubhaloes = 0, TotalDisruptedSubhaloes = 0, TotalNewSubhaloes = 0, TotalFakeSubhaloes = 0;
  MPI_Allreduce(&LocalSunkSubhaloes     , &TotalSunkSubhaloes     , 1, MPI_HBT_INT, MPI_SUM, world.Communicator);
  MPI_Allreduce(&LocalDisruptedSubhaloes, &TotalDisruptedSubhaloes, 1, MPI_HBT_INT, MPI_SUM, world.Communicator);
  MPI_Allreduce(&LocalNewSubhaloes      , &TotalNewSubhaloes      , 1, MPI_HBT_INT, MPI_SUM, world.Communicator);
  MPI_Allreduce(&LocalFakeSubhaloes     , &TotalFakeSubhaloes     , 1, MPI_HBT_INT, MPI_SUM, world.Communicator);

  if(world.rank() == 0)
  {
    std::cout << "    Number of merged subhaloes = " << TotalSunkSubhaloes << std::endl;
    std::cout << "    Number of disrupted subhaloes = " << TotalDisruptedSubhaloes << std::endl;
    std::cout << "    Number of newly identified subhaloes = " << TotalNewSubhaloes << std::endl;
    std::cout << "    Number of FOF groups without any self-bound subhaloes = " << TotalFakeSubhaloes << std::endl;
  }
}