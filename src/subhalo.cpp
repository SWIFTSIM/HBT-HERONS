#include <algorithm>
#include <iostream>
#include <new>
#include <numeric>
#include <omp.h>
#include <limits>

#include "datatypes.h"
#include "snapshot_number.h"
#include "subhalo.h"
#include "particle_exchanger.h"
/*
void MemberShipTable_t::SubIdToTrackId(const SubhaloList_t& Subhalos)
{
   for(HBTInt i=0;i<AllMembers.size();i++)
    AllMembers[i]=Subhalos[AllMembers[i]].TrackId;
}
void MemberShipTable_t::TrackIdToSubId(SubhaloList_t& Subhalos)
{
cout<<"Warning: TrackIdToSubId ToBe fully Implemented!\n";
// exit(1);
}
*/

void ExchangeSubHalos(MpiWorker_t &world, vector<Subhalo_t> &InHalos, vector<Subhalo_t> &OutHalos,
                      MPI_Datatype MPI_Halo_Shell_Type, const ParticleSnapshot_t &snap)
{
  typedef typename vector<Subhalo_t>::iterator HaloIterator_t;
  typedef HaloParticleIterator_t<HaloIterator_t> ParticleIterator_t;
  typedef HaloNestIterator_t<HaloIterator_t> NestIterator_t;

  //   cout<<"Query particle..."<<flush;
  { // query particles
    ParticleExchanger_t<Subhalo_t> Exchanger(world, snap, InHalos);
    Exchanger.Exchange();
  }

  { // do not need to distribute subhalos now. will do it during host-finding anyway.
#pragma omp parallel for
    for (HBTInt i = 0; i < InHalos.size(); i++)
      InHalos[i].AverageCoordinates();
    InHalos.swap(OutHalos);
  }
  /*
  vector <IdRank_t>TargetRank(InHalos.size());
  DecideTargetProcessor(world.size(), InHalos, TargetRank);

  //distribute halo shells
    vector <int> SendHaloCounts(world.size(),0), RecvHaloCounts(world.size()), SendHaloDisps(world.size()),
  RecvHaloDisps(world.size()); sort(TargetRank.begin(), TargetRank.end(), CompareRank); vector <Subhalo_t>
  InHalosSorted(InHalos.size()); vector <HBTInt> InHaloSizes(InHalos.size()), InHaloNestSizes(InHalos.size());
    for(HBTInt haloid=0;haloid<InHalos.size();haloid++)
    {
      InHalosSorted[haloid]=move(InHalos[TargetRank[haloid].Id]);
      SendHaloCounts[TargetRank[haloid].Rank]++;
      InHaloSizes[haloid]=InHalosSorted[haloid].Particles.size();
      InHaloNestSizes[haloid]=InHalosSorted[haloid].NestedSubhalos.size();
    }
    MPI_Alltoall(SendHaloCounts.data(), 1, MPI_INT, RecvHaloCounts.data(), 1, MPI_INT, world.Communicator);
    CompileOffsets(SendHaloCounts, SendHaloDisps);
    HBTInt NumNewHalos=CompileOffsets(RecvHaloCounts, RecvHaloDisps);
    OutHalos.resize(OutHalos.size()+NumNewHalos);
    auto NewHalos=OutHalos.end()-NumNewHalos;
    MPI_Alltoallv(InHalosSorted.data(), SendHaloCounts.data(), SendHaloDisps.data(), MPI_Halo_Shell_Type, &NewHalos[0],
  RecvHaloCounts.data(), RecvHaloDisps.data(), MPI_Halo_Shell_Type, world.Communicator);
  //resize receivehalos
    vector <HBTInt> OutHaloSizes(NumNewHalos), OutHaloNestSizes(NumNewHalos);
    MPI_Alltoallv(InHaloSizes.data(), SendHaloCounts.data(), SendHaloDisps.data(), MPI_HBT_INT, OutHaloSizes.data(),
  RecvHaloCounts.data(), RecvHaloDisps.data(), MPI_HBT_INT, world.Communicator); MPI_Alltoallv(InHaloNestSizes.data(),
  SendHaloCounts.data(), SendHaloDisps.data(), MPI_HBT_INT, OutHaloNestSizes.data(), RecvHaloCounts.data(),
  RecvHaloDisps.data(), MPI_HBT_INT, world.Communicator); for(HBTInt i=0;i<NumNewHalos;i++)
    {
      NewHalos[i].Particles.resize(OutHaloSizes[i]);
      NewHalos[i].NestedSubhalos.resize(OutHaloNestSizes[i]);
    }

    {
    //distribute halo particles
    MPI_Datatype MPI_HBT_Particle;
    Particle_t().create_MPI_type(MPI_HBT_Particle);
    //create combined iterator for each bunch of haloes
    vector <ParticleIterator_t> InParticleIterator(world.size());
    vector <ParticleIterator_t> OutParticleIterator(world.size());
    for(int rank=0;rank<world.size();rank++)
    {
      InParticleIterator[rank].init(InHalosSorted.begin()+SendHaloDisps[rank],
  InHalosSorted.begin()+SendHaloDisps[rank]+SendHaloCounts[rank]);
      OutParticleIterator[rank].init(NewHalos+RecvHaloDisps[rank], NewHalos+RecvHaloDisps[rank]+RecvHaloCounts[rank]);
    }
    vector <HBTInt> InParticleCount(world.size(),0);
    for(HBTInt i=0;i<InHalosSorted.size();i++)
      InParticleCount[TargetRank[i].Rank]+=InHalosSorted[i].Particles.size();

    MyAllToAll<Particle_t, ParticleIterator_t, ParticleIterator_t>(world, InParticleIterator, InParticleCount,
  OutParticleIterator, MPI_HBT_Particle);

    MPI_Type_free(&MPI_HBT_Particle);
    }

    {
    //distribute nests
    //create combined iterator for each bunch of haloes
    vector <NestIterator_t> InNestIterator(world.size());
    vector <NestIterator_t> OutNestIterator(world.size());
    for(int rank=0;rank<world.size();rank++)
    {
      InNestIterator[rank].init(InHalosSorted.begin()+SendHaloDisps[rank],
  InHalosSorted.begin()+SendHaloDisps[rank]+SendHaloCounts[rank]);
      OutNestIterator[rank].init(NewHalos+RecvHaloDisps[rank], NewHalos+RecvHaloDisps[rank]+RecvHaloCounts[rank]);
    }
    vector <HBTInt> InNestCount(world.size(),0);
    for(HBTInt i=0;i<InHalosSorted.size();i++)
      InNestCount[TargetRank[i].Rank]+=InHalosSorted[i].NestedSubhalos.size();

    MyAllToAll<HBTInt, NestIterator_t, NestIterator_t>(world, InNestIterator, InNestCount, OutNestIterator,
  MPI_HBT_INT);
    }
    */
}

void SubhaloSnapshot_t::BuildMPIDataType()
{
  /*to create the struct data type for communication*/
  Subhalo_t p;
  const int MaxNumAttr = 50;
  MPI_Datatype oldtypes[MaxNumAttr];
  int blockcounts[MaxNumAttr];
  MPI_Aint offsets[MaxNumAttr], origin, extent;

  MPI_Get_address(&p, &origin);
  MPI_Get_address((&p) + 1, &extent); // to get the extent of s
  extent -= origin;

  int NumAttr = 0;
#define RegisterAttr(x, type, count)                                                                                   \
  {                                                                                                                    \
    MPI_Get_address(&(p.x), offsets + NumAttr);                                                                        \
    offsets[NumAttr] -= origin;                                                                                        \
    oldtypes[NumAttr] = type;                                                                                          \
    blockcounts[NumAttr] = count;                                                                                      \
    NumAttr++;                                                                                                         \
  }
  RegisterAttr(TrackId, MPI_HBT_INT, 1);
  RegisterAttr(Nbound, MPI_HBT_INT, 1);
  RegisterAttr(Mbound, MPI_FLOAT, 1);
#ifndef DM_ONLY
  RegisterAttr(NboundType, MPI_HBT_INT, TypeMax);
  RegisterAttr(MboundType, MPI_FLOAT, TypeMax);
#endif
  RegisterAttr(TracerIndex, MPI_HBT_INT, 1);
#ifdef CHECK_TRACER_INDEX
  RegisterAttr(TracerId, MPI_HBT_INT, 1);
#endif
  RegisterAttr(HostHaloId, MPI_HBT_INT, 1);
  RegisterAttr(Rank, MPI_HBT_INT, 1);
  RegisterAttr(Depth, MPI_INT, 1);
  RegisterAttr(LastMaxMass, MPI_FLOAT, 1);
  RegisterAttr(SnapshotIndexOfLastMaxMass, MPI_INT, 1);
  RegisterAttr(SnapshotIndexOfLastIsolation, MPI_INT, 1);
  RegisterAttr(SnapshotIndexOfBirth, MPI_INT, 1);
  RegisterAttr(SnapshotIndexOfDeath, MPI_INT, 1);
  RegisterAttr(SnapshotIndexOfSink, MPI_INT, 1);
  RegisterAttr(RmaxComoving, MPI_FLOAT, 1);
  RegisterAttr(RmaxComovingOfLastMaxVmax, MPI_FLOAT, 1);
  RegisterAttr(VmaxPhysical, MPI_FLOAT, 1);
  RegisterAttr(LastMaxVmaxPhysical, MPI_FLOAT, 1);
  RegisterAttr(SnapshotIndexOfLastMaxVmax, MPI_INT, 1);
  RegisterAttr(REncloseComoving, MPI_FLOAT, 1);
  RegisterAttr(RHalfComoving, MPI_FLOAT, 1);
  RegisterAttr(BoundR200CritComoving, MPI_FLOAT, 1);
  // RegisterAttr(R200MeanComoving, MPI_FLOAT, 1);
  // RegisterAttr(RVirComoving, MPI_FLOAT, 1);
  RegisterAttr(BoundM200Crit, MPI_FLOAT, 1);
  // RegisterAttr(M200Mean, MPI_FLOAT, 1);
  // RegisterAttr(MVir, MPI_FLOAT, 1);
  RegisterAttr(SpecificSelfPotentialEnergy, MPI_FLOAT, 1);
  RegisterAttr(SpecificSelfKineticEnergy, MPI_FLOAT, 1);
  RegisterAttr(SpecificAngularMomentum[0], MPI_FLOAT, 3);
#ifdef HAS_GSL
  RegisterAttr(InertialEigenVector[0], MPI_FLOAT, 9);
  RegisterAttr(InertialEigenVectorWeighted[0], MPI_FLOAT, 9);
#endif
  RegisterAttr(InertialTensor[0], MPI_FLOAT, 6);
  RegisterAttr(InertialTensorWeighted[0], MPI_FLOAT, 6);
  RegisterAttr(ComovingAveragePosition[0], MPI_HBT_REAL, 3);
  RegisterAttr(PhysicalAverageVelocity[0], MPI_HBT_REAL, 3);
  RegisterAttr(ComovingMostBoundPosition[0], MPI_HBT_REAL, 3);
  RegisterAttr(PhysicalMostBoundVelocity[0], MPI_HBT_REAL, 3);

  assert(offsets[NumAttr - 1] - offsets[NumAttr - 2] == sizeof(HBTReal) * 3); // to make sure HBTxyz is stored locally.

  RegisterAttr(MostBoundParticleId, MPI_HBT_INT, 1);
  RegisterAttr(SinkTrackId, MPI_HBT_INT, 1);
  RegisterAttr(DescendantTrackId, MPI_HBT_INT, 1);
  RegisterAttr(NestedParentTrackId, MPI_HBT_INT, 1);
#undef RegisterAttr
  assert(NumAttr <= MaxNumAttr);

  MPI_Type_create_struct(NumAttr, blockcounts, offsets, oldtypes, &MPI_HBT_SubhaloShell_t);
  MPI_Type_create_resized(MPI_HBT_SubhaloShell_t, (MPI_Aint)0, extent, &MPI_HBT_SubhaloShell_t);
  MPI_Type_commit(&MPI_HBT_SubhaloShell_t);
}

void SubhaloSnapshot_t::UpdateParticles(MpiWorker_t &world, const ParticleSnapshot_t &snapshot)
{
  /* We need to split particles here, since ExchangeSubHalos updates the particle information of
   * subhaloes based on the Particle IDs present in the Particle vector when it is called.*/
#ifndef DM_ONLY
  if (HBTConfig.ParticlesSplit)
  {
    UpdateSplitParticles(snapshot);
  }
#endif

  Cosmology = snapshot.Cosmology;
  SubhaloList_t LocalSubhalos;
  ExchangeSubHalos(world, Subhalos, LocalSubhalos, MPI_HBT_SubhaloShell_t, snapshot);
  Subhalos.swap(LocalSubhalos);
#pragma omp parallel for
  for (HBTInt i = 0; i < Subhalos.size(); i++)
    Subhalos[i].CountParticles();
}

#ifndef DM_ONLY
/* This function iterates over all local subhaloes, and adds newly spawned particles
 * to the respective Particle vector */
void SubhaloSnapshot_t::UpdateSplitParticles(const ParticleSnapshot_t &snapshot)
{
  for (std::size_t sub_number = 0; sub_number < Subhalos.size(); ++sub_number)
  {
    for (std::size_t particle_number = 0; particle_number < Subhalos[sub_number].Particles.size(); ++particle_number)
    {
      if (snapshot.ParticleSplitMap.count(Subhalos[sub_number].Particles[particle_number].Id) == 1)
      {
        const HBTInt new_id = snapshot.ParticleSplitMap.at(Subhalos[sub_number].Particles[particle_number].Id);

        Particle_t new_particle;
        new_particle.Id = new_id;

        Subhalos[sub_number].Particles.push_back(new_particle);
      }
    }
  }
}
#endif

void Subhalo_t::AverageCoordinates()
{
  /* We do not need to handle orphans here, since their most bound particle
   * position and velocity is set in UpdateMostBoundPosition*/
  if (Particles.size())
  {
    copyHBTxyz(ComovingMostBoundPosition, Particles[GetTracerIndex()].ComovingPosition);
    copyHBTxyz(PhysicalMostBoundVelocity, Particles[GetTracerIndex()].GetPhysicalVelocity());

    AveragePosition(ComovingAveragePosition, Particles.data(), Nbound);
    AverageVelocity(PhysicalAverageVelocity, Particles.data(), Nbound);
  }
}

inline bool CompProfRadius(const RadMassVel_t &a, const RadMassVel_t &b)
{
  return a.r < b.r;
}
inline bool CompProfVel(const RadMassVel_t &a, const RadMassVel_t &b)
{
  return a.v < b.v;
}

void Subhalo_t::CalculateProfileProperties(const Snapshot_t &epoch)
{
  /* to calculate the following density-profile related properties
   *
  HBTReal RmaxComoving;
  HBTReal RmaxComovingOfLastMaxVmax;
  HBTReal VmaxPhysical;
  HBTReal LastMaxVmaxPhysical;
  HBTInt SnapshotIndexOfLastMaxVmax; //the snapshot when it has the maximum Vmax, only considering past snapshots.

  HBTReal REncloseComoving;
  HBTReal RHalfComoving;

  HBTReal R200CritComoving;
  HBTReal R200MeanComoving;
  HBTReal RVirComoving;
  HBTReal M200Crit;
  HBTReal M200Mean;
  HBTReal MVir;
   */
  if (Nbound <= 1)
  {
    RmaxComoving = 0.;
    VmaxPhysical = 0.;
    REncloseComoving = 0.;
    RHalfComoving = 0.;
    BoundR200CritComoving = 0.;
    // 	R200MeanComoving=0.;
    // 	RVirComoving=0.;
    BoundM200Crit = 0.;
    // 	M200Mean=0.;
    // 	MVir=0.;
    return;
  }
  HBTReal VelocityUnit = PhysicalConst::G / epoch.Cosmology.ScaleFactor;

  const HBTxyz &cen = ComovingMostBoundPosition; // most-bound particle as center.

  vector<RadMassVel_t> prof(Nbound);
#pragma omp parallel if (Nbound > 100)
  {
#pragma omp for
    for (HBTInt i = 0; i < Nbound; i++)
    {
      prof[i].r = PeriodicDistance(cen, Particles[i].ComovingPosition);
      prof[i].m = Particles[i].Mass;
    }
#pragma omp single
    {
      sort(prof.begin(), prof.end(), CompProfRadius);
      double m_cum = 0.;
      for (auto &&p : prof)
        p.m = (m_cum += p.m);
    }
#pragma omp for
    for (HBTInt i = 0; i < Nbound; i++)
    {
      if (prof[i].r < HBTConfig.GetCurrentSoftening(epoch.Cosmology.ScaleFactor))
        prof[i].r = HBTConfig.GetCurrentSoftening(epoch.Cosmology.ScaleFactor); // resolution
      prof[i].v = prof[i].m / prof[i].r;                                        // v**2
    }
  }
  auto maxprof = max_element(prof.begin(), prof.end(), CompProfVel);
  RmaxComoving = maxprof->r;
  VmaxPhysical = sqrt(maxprof->v * VelocityUnit);
  RHalfComoving = prof[Nbound / 2].r;
  REncloseComoving = prof[Nbound - 1].r;

  epoch.SphericalOverdensitySize(BoundM200Crit, BoundR200CritComoving, 200., prof);

  if (VmaxPhysical >= LastMaxVmaxPhysical)
  {
    SnapshotIndexOfLastMaxVmax = epoch.GetSnapshotIndex();
    LastMaxVmaxPhysical = VmaxPhysical;
    RmaxComovingOfLastMaxVmax = RmaxComoving;
  }
}

void Subhalo_t::CalculateShape()
{
  if (Nbound <= 1)
  {
#ifdef HAS_GSL
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
      {
        InertialEigenVector[i][j] = 0.;
        InertialEigenVectorWeighted[i][j] = 0.;
      }
#endif
    for (auto &&I : InertialTensor)
      I = 0.;
    for (auto &&I : InertialTensorWeighted)
      I = 0.;
    return;
  }
  const HBTxyz &cen = ComovingMostBoundPosition; // most-bound particle as center.

  double Ixx = 0, Iyy = 0, Izz = 0, Ixy = 0, Ixz = 0, Iyz = 0;
  double Ixxw = 0, Iyyw = 0, Izzw = 0, Ixyw = 0, Ixzw = 0, Iyzw = 0;
#pragma omp parallel for reduction(+ : Ixx, Iyy, Izz, Ixy, Ixz, Iyz, Ixxw, Iyyw, Izzw, Ixyw, Ixzw,                     \
                                     Iyzw) if (Nbound > 100)
  for (HBTInt i = 1; i < Nbound; i++)
  {
    HBTReal m = Particles[i].Mass;
    const HBTxyz &pos = Particles[i].ComovingPosition;
    HBTReal dx = pos[0] - cen[0];
    HBTReal dy = pos[1] - cen[1];
    HBTReal dz = pos[2] - cen[2];
    if (HBTConfig.PeriodicBoundaryOn)
    {
      dx = NEAREST(dx);
      dy = NEAREST(dy);
      dz = NEAREST(dz);
    }
    HBTReal dx2 = dx * dx;
    HBTReal dy2 = dy * dy;
    HBTReal dz2 = dz * dz;
    Ixx += dx2 * m;
    Iyy += dy2 * m;
    Izz += dz2 * m;
    Ixy += dx * dy * m;
    Ixz += dx * dz * m;
    Iyz += dy * dz * m;

    HBTReal dr2 = dx2 + dy2 + dz2;
    dr2 /= m; // for mass weighting
    Ixxw += dx2 / dr2;
    Iyyw += dy2 / dr2;
    Izzw += dz2 / dr2;
    Ixyw += dx * dy / dr2;
    Ixzw += dx * dz / dr2;
    Iyzw += dy * dz / dr2;
  }
  InertialTensor[0] = Ixx;
  InertialTensor[1] = Ixy;
  InertialTensor[2] = Ixz;
  InertialTensor[3] = Iyy;
  InertialTensor[4] = Iyz;
  InertialTensor[5] = Izz;
  InertialTensorWeighted[0] = Ixxw;
  InertialTensorWeighted[1] = Ixyw;
  InertialTensorWeighted[2] = Ixzw;
  InertialTensorWeighted[3] = Iyyw;
  InertialTensorWeighted[4] = Iyzw;
  InertialTensorWeighted[5] = Izzw;
  for (auto &&I : InertialTensor)
    I /= Mbound;
  for (auto &&I : InertialTensorWeighted)
    I /= Mbound;
#ifdef HAS_GSL
  EigenAxis(Ixx, Ixy, Ixz, Iyy, Iyz, Izz, InertialEigenVector);
  EigenAxis(Ixxw, Ixyw, Ixzw, Iyyw, Iyzw, Izzw, InertialEigenVectorWeighted);
#endif
}

void Subhalo_t::CountParticleTypes()
{
#ifndef DM_ONLY

  /* Used for finding TracerIndex for large (Nbound > 100) tracks */
#pragma omp declare reduction(minPair:IndexParticleType_t : omp_out = firstIndex(omp_out, omp_in))                     \
  initializer(omp_priv = IndexParticleType_t(numeric_limits<HBTInt>::max(), -1))

  /* Initialise the value of the pair to identify specified tracer */
  IndexParticleType_t Tracer_Index_ParticleType(numeric_limits<HBTInt>::max(), -1);

  for (int itype = 0; itype < TypeMax; itype++)
  {
    NboundType[itype] = 0;
    MboundType[itype] = 0.;
  }
  if (Nbound > 100) // parallelize
  {
#pragma omp parallel
    {
      vector<HBTInt> nboundtype(TypeMax, 0);
      vector<float> mboundtype(TypeMax, 0.);
#pragma omp for reduction(minPair : Tracer_Index_ParticleType)
      for (HBTInt i = 0; i < Nbound; i++)
      {
        auto &p = Particles[i];
        // 		if(p.Id==SpecialConst::NullParticleId) continue;
        int itype = p.Type;
        nboundtype[itype]++;
        mboundtype[itype] += p.Mass;

        /* Check whether we found a collisionless particle for the first time,
         * indicating we have reached the tentative most bound collisionless
         * tracer */
        if (Tracer_Index_ParticleType.second == -1) // No tracer yet
          if (p.IsTracer())                         // Found tracer
          {
            Tracer_Index_ParticleType.first = i;
            Tracer_Index_ParticleType.second = itype;
          }
      }
#pragma omp critical
      for (int i = 0; i < TypeMax; i++)
      {
        NboundType[i] += nboundtype[i];
        MboundType[i] += mboundtype[i];
      }
    }
  }
  else
  {
    auto end = Particles.begin() + Nbound;
    for (auto it = Particles.begin(); it != end; ++it)
    {
      auto &p = *it;
      if (p.Id == SpecialConst::NullParticleId)
        continue;
      int itype = p.Type;
      NboundType[itype]++;
      MboundType[itype] += p.Mass;

      /* Check whether we found a collisionless particle for the first time,
       * indicating we have reached the tentative most bound collisionless
       * tracer */
      if (Tracer_Index_ParticleType.second == -1) // No tracer yet
        if (p.IsTracer())                         // Found tracer
        {
          Tracer_Index_ParticleType.first = it - Particles.begin();
          Tracer_Index_ParticleType.second = itype;
        }
    }
  }

  // If we found no tracer in this subgroup, default to the most bound particle.
  HBTInt index = (Tracer_Index_ParticleType.second == -1) ? 0 : Tracer_Index_ParticleType.first;
  SetTracerIndex(index);

#else
  // Always use the first particle in DMO runs
  SetTracerIndex(0);
#endif
}

HBTInt Subhalo_t::KickNullParticles()
{
#ifdef DM_ONLY
  return 0;
#else
  HBTInt np_old = Particles.size();
  auto it_begin = Particles.begin(), it_save = it_begin, it = it_begin, it_end = it_begin + Nbound;
  for (; it != it_end; ++it)
  {
    if (it->Id != SpecialConst::NullParticleId) // there will be consumed particles
    {
      if (it != it_save)
        *it_save = move(*it);
      ++it_save;
    }
  }
  Nbound = it_save - it_begin;

  it_end = Particles.end();
  for (; it != it_end; ++it) // unbound particles
  {
    if (it->Id != SpecialConst::NullParticleId) // there will be consumed particles
    {
      if (it != it_save)
        *it_save = move(*it);
      ++it_save;
    }
  }
  Particles.resize(it_save - it_begin);

  // Temporary fix to update TracerIndex
  CountParticleTypes();

  // if(it!=it_save) cout<<it-it_save<<" outof "<<np_old<<" particles consumed for track "<<TrackId<<"\n";
  return it - it_save;
#endif
}

/*update Mbound, NboundType, MboundType *
 * this function is called during unbinding, merger and BH consumption(UpdateParticles)*
 */
void Subhalo_t::CountParticles()
{
#ifdef DM_ONLY
  Mbound = 0.;
  for (HBTInt i = 0; i < Nbound; i++)
    Mbound += Particles[i].Mass;
#else
  for (auto &n : NboundType)
    n = 0;
  for (auto &m : MboundType)
    m = 0.;
  auto it = Particles.begin(), it_end = Particles.begin() + Nbound;
  for (; it != it_end; ++it)
  {
    int itype = it->Type;
    NboundType[itype]++;
    MboundType[itype] += it->Mass;
  }
  Mbound = accumulate(begin(MboundType), end(MboundType), (HBTReal)0.);
#endif
}

/*
  Return the IDs of the N most bound tracer particles.
  Will return non-tracers if there are not enough tracers.
*/
vector<HBTInt> Subhalo_t::GetMostBoundTracerIds(HBTInt n)
{

  // Allocate the output vector
  std::vector<HBTInt> Ids(0);
  Ids.reserve(n);

  // Store IDs of tracers
  for (int pass_nr = 0; pass_nr < 2; pass_nr += 1)
  {
    for (HBTInt i = 0; i < Nbound; i++)
    {
      const int is_tracer = Particles[i].IsTracer();
      if ((is_tracer && (pass_nr == 0)) || ((!is_tracer) && (pass_nr == 1)))
        Ids.push_back(Particles[i].Id);
      if (Ids.size() == n)
        break;
    }
    if (Ids.size() == n)
      break;
  }

  /* We should have found sufficient tracers*/
  assert(Ids.size() == n);

  return Ids;
}
