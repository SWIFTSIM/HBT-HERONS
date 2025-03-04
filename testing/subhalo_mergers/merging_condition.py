#!/bin/env python
#
# Compute quantities used to determine subhalo mergers in HBTplus
#

import os
import os.path

import numpy as np
import h5py
import unyt

import virgo.mpi.parallel_hdf5 as phdf5
import virgo.mpi.parallel_sort as psort
import virgo.util.partial_formatter as pf
import virgo.formats.swift

from mpi4py import MPI
comm = MPI.COMM_WORLD
comm_rank = comm.Get_rank()
comm_size = comm.Get_size()

nr_tracers = 20
tracer_types = (1, 4)
#tracer_types = (1,)


def destination_rank(arr, comm_size):
    """
    Assign elements of integer array arr to MPI ranks.
    """
    if arr.dtype.itemsize == 4:
        arr_view = arr.view(dtype=np.uint32)
        arr_hash = np.bitwise_xor(np.right_shift(arr_view, 16), arr_view)
        arr_hash *= 0x45d9f3b
        arr_hash = np.bitwise_xor(np.right_shift(arr_hash, 16), arr_hash)
        arr_hash *= 0x45d9f3b
        arr_hash = np.bitwise_xor(np.right_shift(arr_hash, 16), arr_hash)
    elif arr.dtype.itemsize == 8:
        arr_view = arr.view(dtype=np.uint64)
        arr_hash = np.bitwise_xor(arr_view, np.right_shift(arr_view, 30))
        arr_hash *= 0xbf58476d1ce4e5b9
        arr_hash = np.bitwise_xor(arr_hash, np.right_shift(arr_hash, 27))
        arr_hash *= 0x94d049bb133111eb
        arr_hash = np.bitwise_xor(arr_hash, np.right_shift(arr_hash, 31))
    else:
        raise RuntimeError("Unsupported data type: must be 4 or 8 bytes per element")
    return np.mod(arr_hash, comm_size).astype(int)


def message(s):
    if comm_rank == 0:
        print(s)


def read_snapshot(snapshot_format, membership_format, snap_nr):
    """
    Read particle data from the snapshot and membership files.
    
    Here we need position, velocity, bound rank and group number.
    Need to read and concatenate data from all particle types.
    """

    # First need to determine number of snapshot files
    if comm_rank == 0:
        filename = snapshot_format.format(file_nr=0, snap_nr=snap_nr)
        with h5py.File(filename, "r") as infile:
            nr_files = infile["Header"].attrs["NumFilesPerSnapshot"][0]
    else:
        nr_files = None
    nr_files = comm.bcast(nr_files)
    
    # Read the snapshot files
    formatter = pf.PartialFormatter()
    snap_filename = formatter.format(snapshot_format, snap_nr=snap_nr, file_nr=None)
    snap = phdf5.MultiFile(snap_filename, file_idx=np.arange(nr_files, dtype=int), comm=comm)
    pos = []
    vel = []
    mass = []
    ptype = []
    for type_nr in tracer_types:
        data = snap.read(("Coordinates", "Velocities", "Masses"), group=f"PartType{type_nr}", unpack=True)
        pos.append(data[0])
        vel.append(data[1])
        mass.append(data[2])
        ptype.append(np.ones(pos[-1].shape[0], dtype=np.int32) * type_nr)

    # Read the membership files
    memb_filename = formatter.format(membership_format, snap_nr=snap_nr, file_nr=None)
    memb = phdf5.MultiFile(memb_filename, file_idx=np.arange(nr_files, dtype=int), comm=comm)
    bound_rank = []
    groupnr_bound = []
    for type_nr in tracer_types:
        data = memb.read(("GroupNr_bound", "Rank_bound"), group=f"PartType{type_nr}", unpack=True)
        groupnr_bound.append(data[0])
        bound_rank.append(data[1])
    
    # Combine arrays
    pos = np.concatenate(pos)
    vel = np.concatenate(vel)
    mass = np.concatenate(mass)
    ptype = np.concatenate(ptype)
    groupnr_bound = np.concatenate(groupnr_bound)
    bound_rank = np.concatenate(bound_rank)

    # Discard particles not bound to any halo
    keep = (bound_rank >= 0)
    pos = pos[keep,:]
    vel = vel[keep,:]
    mass = mass[keep]
    ptype = ptype[keep]
    groupnr_bound = groupnr_bound[keep]
    bound_rank = bound_rank[keep]

    return pos, vel, mass, groupnr_bound, bound_rank, ptype


def read_hbt(hbt_basedir, snap_nr):
    """
    Read the HBT halo catalogue
    """

    hbt_filename = hbt_basedir + f"/{snap_nr:03d}/SubSnap_{snap_nr:03d}" + ".{file_nr}.hdf5"
    mf = phdf5.MultiFile(hbt_filename, file_nr_dataset="NumberOfFiles", comm=comm)
    subhalos = mf.read("Subhalos")
    return subhalos
    

def compute_merging_condition(snapshot_format, membership_format, hbt_basedir, output_dir, snap_nr):
    """
    Compute position and radius for each subhalo in phase space
    """

    # Read box size and units
    length_in_mpch = None
    vel_in_kms = None
    mass_in_msunh = None
    if comm_rank == 0:
        filename = f"{hbt_basedir}/Parameters.log"
        with open(filename) as f:
            for line in f:
                fields = line.rstrip().split()
                if len(fields) < 2:
                    continue
                if fields[0] == "LengthInMpch":
                    length_in_mpch = float(fields[1])
                if fields[0] == "VelInKmS":
                    vel_in_kms = float(fields[1])
                if fields[0] == "MassInMsunh":
                    mass_in_msunh = float(fields[1])
    length_in_mpch, vel_in_kms, mass_in_msunh = comm.bcast((length_in_mpch, vel_in_kms, mass_in_msunh))

    # Set up swift unit registry by reading a snapshot file
    reg = None
    boxsize = None
    if comm_rank == 0:
        filename = snapshot_format.format(file_nr=0, snap_nr=snap_nr)
        with h5py.File(filename, "r") as snap:
            reg = virgo.formats.swift.soap_unit_registry_from_snapshot(snap)
            boxsize = snap["Header"].attrs["BoxSize"][0]
    reg, boxsize = comm.bcast((reg, boxsize))
    #boxsize = boxsize * unyt.Unit("snap_length", registry=reg) # Get box size with units
    message(f"Box size = {boxsize}")

    # Determine HBT units using Swift's Mpc and Msun definitions
    hbt_length   = length_in_mpch*unyt.Unit("a*swift_mpc/h", registry=reg)
    hbt_velocity = vel_in_kms*unyt.Unit("km/s", registry=reg)
    hbt_mass     = mass_in_msunh*unyt.Unit("swift_msun/h", registry=reg)
    message(f"HBT length unit = {hbt_length}")
    message(f"HBT velocity unit = {hbt_velocity}")
    message(f"HBT mass unit = {hbt_mass}")
    
    # Read in the HBT catalogue
    subhalos = read_hbt(hbt_basedir, snap_nr)
    nr_subs = len(subhalos)
    nr_subs_tot = comm.allreduce(nr_subs)
    message(f"Total number of subhalos read = {nr_subs_tot}")
    
    # Assign HBT halos to MPI ranks according to hash of their index
    nr_local_subhalos = len(subhalos)
    local_subhalo_offset = comm.scan(nr_local_subhalos) - nr_local_subhalos
    subhalo_groupnr = np.arange(nr_local_subhalos, dtype=int) + local_subhalo_offset
    subhalo_rank = destination_rank(subhalo_groupnr, comm_size)
    message("Assigned destination ranks to subhalos")
    
    # Sort HBT halos by destination rank
    sort_key = subhalo_rank.copy()
    order = psort.parallel_sort(sort_key, return_index=True, comm=comm)
    subhalos = psort.fetch_elements(subhalos, order, comm=comm)
    subhalo_groupnr = psort.fetch_elements(subhalo_groupnr, order, comm=comm)
    subhalo_rank = psort.fetch_elements(subhalo_rank, order, comm=comm)
    message("Sorted subhalos")
    
    # Count subhalos to go to each rank
    nr_per_rank = np.bincount(subhalo_rank, minlength=comm_size)
    comm.Allreduce(MPI.IN_PLACE, nr_per_rank)
    nr_total_subhalos = comm.allreduce(nr_local_subhalos)
    assert np.sum(nr_per_rank) == nr_total_subhalos
    message("Counted subhalos per rank")

    # Repartition subhalos. After this they should all be on the assigned MPI rank.
    subhalos = psort.repartition(subhalos, ndesired=nr_per_rank, comm=comm)
    subhalo_groupnr = psort.repartition(subhalo_groupnr, ndesired=nr_per_rank, comm=comm)
    subhalo_rank = psort.repartition(subhalo_rank, ndesired=nr_per_rank, comm=comm)
    assert np.all(subhalo_rank==comm_rank)
    del subhalo_rank
    message("Repartitioned subhalos")
    
    # Read in the bound particles
    pos, vel, mass, groupnr_bound, bound_rank, ptype = read_snapshot(snapshot_format, membership_format, snap_nr)
    nr_parts = pos.shape[0]
    nr_parts_tot = comm.allreduce(nr_parts)
    message(f"Total number of particles read = {nr_parts_tot}")

    # Assign particles to MPI rank according to hash of their group index
    particle_rank = destination_rank(groupnr_bound, comm_size)
    
    # Sort particles by destination rank
    order = psort.parallel_sort(particle_rank, return_index=True, comm=comm)
    pos = psort.fetch_elements(pos, order, comm=comm)
    vel = psort.fetch_elements(vel, order, comm=comm)
    mass = psort.fetch_elements(mass, order, comm=comm)
    groupnr_bound = psort.fetch_elements(groupnr_bound, order, comm=comm)
    bound_rank = psort.fetch_elements(bound_rank, order, comm=comm)
    ptype = psort.fetch_elements(ptype, order, comm=comm)
    message("Sorted particles by destination")
    
    # Count particles to go to each rank
    nr_per_rank = np.bincount(particle_rank, minlength=comm_size)
    comm.Allreduce(MPI.IN_PLACE, nr_per_rank)
    nr_total_particles = comm.allreduce(nr_parts)
    assert np.sum(nr_per_rank) == nr_total_particles
    message("Counted particles per rank")

    # Repartition particles: should leave all on the same rank as their subhalo
    particle_rank = psort.repartition(particle_rank, ndesired=nr_per_rank, comm=comm)
    pos = psort.repartition(pos, ndesired=nr_per_rank, comm=comm)
    vel = psort.repartition(vel, ndesired=nr_per_rank, comm=comm)
    mass = psort.repartition(mass, ndesired=nr_per_rank, comm=comm)
    groupnr_bound = psort.repartition(groupnr_bound, ndesired=nr_per_rank, comm=comm)
    bound_rank = psort.repartition(bound_rank, ndesired=nr_per_rank, comm=comm)
    ptype = psort.repartition(ptype, ndesired=nr_per_rank, comm=comm)
    assert np.all(particle_rank==comm_rank)
    del particle_rank
    
    # Now need to locally sort particles by group membership and then bound rank
    sort_key_t = np.dtype([("groupnr_bound", groupnr_bound.dtype),
                           ("bound_rank", bound_rank.dtype)])
    sort_key = np.ndarray(len(ptype), dtype=sort_key_t)
    sort_key["groupnr_bound"] = groupnr_bound
    sort_key["bound_rank"] = bound_rank
    order = np.argsort(sort_key)
    pos = pos[order,:]
    vel = vel[order,:]
    mass = mass[order]
    groupnr_bound = groupnr_bound[order]
    bound_rank = bound_rank[order]
    ptype = ptype[order]
    message("Locally sorted particles")

    # Then, for each local subhalo find the corresponding range of particles
    first_particle_in_this_sub = np.searchsorted(groupnr_bound, subhalo_groupnr, side="left")
    first_particle_in_next_sub = np.searchsorted(groupnr_bound, subhalo_groupnr, side="right")
    nr_particles_in_sub = first_particle_in_next_sub - first_particle_in_this_sub
    message("Identified range of particles in each subhalo")

    # Sanity check: the number of particles found should equal the total
    # number of tracer type particles in the subhalo. Compute expected number
    # of tracers in each local subhalo.
    subhalo_nr_tracers = np.zeros_like(subhalos["Nbound"])
    if "NboundType" in subhalos.dtype.fields:
        # Hydro run, so sum tracer types only
        for tt in tracer_types:
            subhalo_nr_tracers += subhalos["NboundType"][:,tt]
    else:
        # DMO run, so assume all particles are tracers because we don't have NboundType
        subhalo_nr_tracers = subhalos["Nbound"]
    assert np.all(subhalo_nr_tracers == nr_particles_in_sub)
    message("Found expected number of tracer particles in each subhalo")

    # Check that we identified particles in halos correctly
    for i in range(len(subhalos)):
        i1 = first_particle_in_this_sub[i]
        i2 = first_particle_in_next_sub[i]
        assert np.all(groupnr_bound[i1:i2] == subhalo_groupnr[i])
    message("Subhalo bound particles have expected groupnr_bound")
    
    # Add units to particle arrays from the snapshot
    #pos = pos * unyt.Unit("snap_length", registry=reg)
    #vel = vel * unyt.Unit("snap_length/snap_time", registry=reg)
    #mass = mass * unyt.Unit("snap_mass", registry=reg)

    # Allocate output arrays
    core_mean_pos = np.zeros((len(subhalos),3), dtype=float)
    core_sigma_r  = np.zeros(len(subhalos), dtype=float)
    core_mean_vel = np.zeros((len(subhalos),3), dtype=float)
    core_sigma_v  = np.zeros(len(subhalos), dtype=float)
    
    # Loop over all local subhalos
    for i in range(len(subhalos)):
        i1 = first_particle_in_this_sub[i]
        i2 = min(first_particle_in_next_sub[i], i1+nr_tracers)
        if i2 > i1:
            # Find particles in this subhalo
            sub_part_pos  = pos[i1:i2,:].astype(float)
            sub_part_vel  = vel[i1:i2,:].astype(float)
            sub_part_mass = mass[i1:i2].astype(float)
            # Pick a reference point for this subhalo
            ref_pos = sub_part_pos[0,:]
            # Wrap particles to periodic copy closest to the reference point
            offset = (0.5*boxsize) - ref_pos # offset to place ref_pos at box centre
            sub_part_pos = ((sub_part_pos + offset) % boxsize) - offset
            # Get total mass
            mtot = np.sum(sub_part_mass, dtype=float)
            # Compute mass weighted mean position
            dx = sub_part_pos #- ref_pos
            sum_m_dx = np.sum(sub_part_mass[:,None]*dx, axis=0, dtype=float)
            pos_mean = sum_m_dx / mtot #+ ref_pos
            # Compute dispersion in position
            sum_m_dx2 = np.sum(sub_part_mass[:,None]*dx**2, axis=0, dtype=float)
            pos_disp = sum_m_dx2 / mtot
            pos_disp -= pos_mean**2
            # Compute mass weighted mean velocity
            sum_m_dv = np.sum(sub_part_mass[:,None]*sub_part_vel, axis=0, dtype=float)
            vel_mean = sum_m_dv / mtot
            # Compute dispersion in velocity
            sum_m_dv2 = np.sum(sub_part_mass[:,None]*sub_part_vel**2, axis=0, dtype=float)
            vel_disp = sum_m_dv2 / mtot
            vel_disp -= vel_mean**2
            # Wrap position back into the box
            pos_mean = pos_mean % boxsize
            # Store results
            core_mean_pos[i,:] = pos_mean
            core_sigma_r[i] = np.sqrt(np.sum(pos_disp))
            core_mean_vel[i,:] = vel_mean
            core_sigma_v[i] = np.sqrt(np.sum(vel_disp))

    message("Computed positions and dispersions")
            
    # Put output into input order
    order = psort.parallel_sort(subhalo_groupnr, return_index=True, comm=comm)
    core_mean_pos = psort.fetch_elements(core_mean_pos, order, comm=comm)
    core_sigma_r  = psort.fetch_elements(core_sigma_r, order, comm=comm)    
    core_mean_vel = psort.fetch_elements(core_mean_vel, order, comm=comm)
    core_sigma_v  = psort.fetch_elements(core_sigma_v, order, comm=comm)
    nbound        = psort.fetch_elements(subhalos["Nbound"], order, comm=comm)
    trackid       = psort.fetch_elements(subhalos["TrackId"], order, comm=comm)
    parent_trackid = psort.fetch_elements(subhalos["NestedParentTrackId"], order, comm=comm)
    
    message("Sorted results")
            
    # Write out the results
    filename = f"{output_dir}/merging_info_{snap_nr:03d}.hdf5"
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    with h5py.File(filename, "w", driver="mpio", comm=comm) as f:
        f["BoxSize"] = boxsize
        phdf5.collective_write(f, "CoreComovingPosition", core_mean_pos, comm)
        phdf5.collective_write(f, "CoreComovingSigmaR",   core_sigma_r, comm)
        phdf5.collective_write(f, "CorePhysicalVelocity", core_mean_vel, comm)
        phdf5.collective_write(f, "CorePhysicalSigmaV",   core_sigma_v, comm)
        phdf5.collective_write(f, "Nbound",               nbound, comm)
        phdf5.collective_write(f, "TrackId",              trackid, comm)
        phdf5.collective_write(f, "NestedParentTrackId",  parent_trackid, comm)

    comm.barrier()
    message("Done")

    
if __name__ == "__main__":

    from virgo.mpi.util import MPIArgumentParser

    parser = MPIArgumentParser(comm, description='Compute HBTplus merging criteria.')
    parser.add_argument('snapshot_format', help='Format string for snapshot filenames (using {snap_nr}, {file_nr})')
    parser.add_argument('membership_format', help='Format string for halo membership filenames (using {snap_nr}, {file_nr})')
    parser.add_argument('hbt_basedir', help='Location of the HBTplus output')
    parser.add_argument('output_dir', help='Where to write the output')
    parser.add_argument('snap_nr', type=int, help='Snapshot number to do')
    args = parser.parse_args()

    compute_merging_condition(**vars(args))
    
