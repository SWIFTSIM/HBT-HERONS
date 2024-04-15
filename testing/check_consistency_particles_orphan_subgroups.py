#!/bin/env python
from tqdm import tqdm
from mpi4py import MPI
comm = MPI.COMM_WORLD
comm_rank = comm.Get_rank()
comm_size = comm.Get_size()

import h5py
import numpy as np

import virgo.mpi.util
import virgo.mpi.parallel_hdf5 as phdf5
import virgo.mpi.parallel_sort as psort

def read_snapshot(snapshot_file, snap_nr, particle_ids):
    """
    Read particle properties for the specified particle IDs.
    Returns a dict of arrays.
    """

    # Datasets to pass through from the snapshot
    passthrough_datasets = ("FOFGroupIDs",)
    
    # Sub in the snapshot number
    from virgo.util.partial_formatter import PartialFormatter
    formatter = PartialFormatter()
    filenames = formatter.format(snapshot_file, snap_nr=snap_nr, file_nr=None)

    # Determine what particle types we have in the snapshot
    if comm_rank == 0:
        ptypes = []
        with h5py.File(filenames.format(file_nr=0), "r") as infile:
            nr_types = int(infile["Header"].attrs["NumPartTypes"])
            nr_parts = infile["Header"].attrs["NumPart_Total"]
            nr_parts_hw = infile["Header"].attrs["NumPart_Total_HighWord"]
            for i in range(nr_types):
                if nr_parts[i] > 0 or nr_parts_hw[i] > 0:
                    ptypes.append(i)
    else:
        ptypes = None
    ptypes = comm.bcast(ptypes)

    # Read the particle data from the snapshot
    particle_data = {"Type" : -np.ones(particle_ids.shape, dtype=np.int32)}
    mf = phdf5.MultiFile(filenames, file_nr_attr=("Header","NumFilesPerSnapshot"), comm=comm)
    for ptype in ptypes:

        if ptype == 6:
            continue # skip neutrinos
        
        # Read the IDs of this particle type
        if comm_rank == 0:
            print(f"Reading snapshot particle IDs for type {ptype}")

        snapshot_ids = mf.read(f"PartType{ptype}/ParticleIDs")

        # For each subhalo particle ID, find matching index in the snapshot (if any)
        ptr = psort.parallel_match(particle_ids, snapshot_ids, comm=comm)
        matched = (ptr>=0)

        # Loop over particle properties to pass through
        for name in passthrough_datasets:

            # Read this property from the snapshot
            snapshot_data = mf.read(f"PartType{ptype}/{name}")

            # Allocate output array, if we didn't already
            if name not in particle_data:
                shape = (len(particle_ids),)+snapshot_data.shape[1:]
                dtype = snapshot_data.dtype
                particle_data[name] = -np.ones(shape, dtype=dtype) # initialize to -1 = not found

            # Look up the value for each subhalo particle
            particle_data[name][matched,...] = psort.fetch_elements(snapshot_data, ptr[matched], comm=comm)
        
        # Also store the type of each matched particle
        particle_data["Type"][matched] = ptype

    # Should have matched all particles, except for merged black holes 
    # assert np.all(particle_data["Type"] >= 0)
        
    return particle_data
    

def read_particles(filenames, nr_local_subhalos):
    """
    Read in the particle IDs belonging to the subhalos on this MPI
    rank from the specified SubSnap files. Returns a single array with
    the concatenated IDs from all local subhalos in the order they
    appear in the SubSnap files.
    """

    # First determine how many subhalos are in each SubSnap file
    if comm_rank == 0:
        subhalos_per_file = []
        nr_files = 1
        file_nr = 0
        while file_nr < nr_files:
            with h5py.File(filenames.format(file_nr=file_nr), "r") as infile:
                subhalos_per_file.append(infile["Subhalos"].shape[0])
                nr_files = int(infile["NumberOfFiles"][...])
            file_nr += 1
    else:
        subhalos_per_file = None
        nr_files = None
    nr_files = comm.bcast(nr_files)
    subhalos_per_file = np.asarray(comm.bcast(subhalos_per_file), dtype=int)
    first_subhalo_in_file = np.cumsum(subhalos_per_file) - subhalos_per_file
    
    # Determine offset to first subhalo this rank reads
    first_local_subhalo = comm.scan(nr_local_subhalos) - nr_local_subhalos

    # Loop over all files
    particle_ids = []
    for file_nr in range(nr_files):

        # Find range of subhalos this rank read from this file
        i1 = first_local_subhalo - first_subhalo_in_file[file_nr]
        i2 = i1 + nr_local_subhalos
        i1 = max(0, i1)
        i2 = min(subhalos_per_file[file_nr], i2)

        # Read subhalo particle IDs, if there are any in this file for this rank
        if i2 > i1:
            with h5py.File(filenames.format(file_nr=file_nr), "r") as infile:
                particle_ids.append(infile["SubhaloParticles"][i1:i2])

    if len(particle_ids) > 0:
        # Combine arrays from different files
        particle_ids = np.concatenate(particle_ids)
        # Combine arrays from different subhalos
        particle_ids = np.concatenate(particle_ids)        
    else:
        # Some ranks may have read zero files
        particle_ids = None
    particle_ids = virgo.mpi.util.replace_none_with_zero_size(particle_ids, comm=comm)
    
    return particle_ids

def sort_hbt_output(basedir, hbt_nr, snap_nr, snapshot_file):

    """
    This reorganizes a set of HBT SubSnap files into a single file which
    contains one HDF5 dataset for each subhalo property. Subhalos are written
    in order of TrackId.

    Particle IDs in groups can be optionally copied to the output.
    """

    # Read in the input subhalos
    if comm_rank == 0:
        print(f"Testing HBTplus consistency of MostBoundParticleId between snapshot index {hbt_nr} and {hbt_nr + 1}.")

    #===========================================================================
    # Load catalogues for snapshot N 
    #===========================================================================
    
    # Make a format string for the filenames
    filenames = f"{basedir}/{hbt_nr:03d}/SubSnap_{hbt_nr:03d}" + ".{file_nr}.hdf5"
    if comm_rank ==0:
        print(f"Opening HBT catalogue: {filenames}", end=' --- ')

    # Open file and load
    mf = phdf5.MultiFile(filenames, file_nr_dataset="NumberOfFiles", comm=comm)
    subhalos_before = mf.read("Subhalos")

    field_names = list(subhalos_before.dtype.fields)

    # Convert array of structs to dict of arrays
    data_before = {}
    for name in field_names:
        data_before[name] = np.ascontiguousarray(subhalos_before[name])
    del subhalos_before

    if comm_rank == 0:
        print("DONE")

    #===========================================================================
    # Load catalogues for snapshot N + 1
    #===========================================================================
    
    # Make a format string for the filenames
    filenames = f"{basedir}/{hbt_nr + 1:03d}/SubSnap_{hbt_nr + 1:03d}" + ".{file_nr}.hdf5"
    if comm_rank == 0:
        print(f"Opening HBT catalogue: {filenames}", end=' --- ')

    mf = phdf5.MultiFile(filenames, file_nr_dataset="NumberOfFiles", comm=comm)
    subhalos_after = mf.read("Subhalos")

    # Keep the ones that are currently orphans, not formed through sinking. If
    # the last condition is not present, we'd get subgroups deemed to be bound in
    # snapshot N + 1, and hence their MostBoundParticleId would change to reflect 
    # this
    subhalos_after = np.array([sub for sub in subhalos_after \
                               if (sub['Nbound'] <= 1) & (sub['SnapshotIndexOfSink'] != (hbt_nr+1))]) 

    field_names = list(subhalos_after.dtype.fields)

    # Convert array of structs to dict of arrays
    data_after = {}
    for name in field_names:
        data_after[name] = np.ascontiguousarray(subhalos_after[name])
    del subhalos_after

    if comm_rank == 0:
        print("DONE")

    #===========================================================================
    # Sort catalogue N, and compare what their MostBoundParticleIds are compared
    # to catalogue N+1.
    #===========================================================================

    # Establish TrackId ordering for the subhalos
    order = psort.parallel_sort(data_before["TrackId"], return_index=True, comm=comm)

    # Sort the subhalo properties by TrackId
    for name in field_names:
        if name != "TrackId":
            data_before[name] = psort.fetch_elements(data_before[name], order, comm=comm)
    del order

    order = data_after['TrackId']
    mostboundids  = psort.fetch_elements(data_before['MostBoundParticleId'], order, comm=comm)

    # Check across all ranks across
    local_number_disagreements = np.sum((mostboundids != data_after['MostBoundParticleId']))
    total_number_disagreements = comm.allreduce(local_number_disagreements)
    total_number_checks = comm.allreduce(len(data_after['MostBoundParticleId']))
    
    if comm_rank == 0:
        print(f"{total_number_disagreements} out of {total_number_checks} orphans disagree.")                

if __name__ == "__main__":

    from virgo.mpi.util import MPIArgumentParser
    
    parser = MPIArgumentParser(comm, description="Reorganize HBTplus SubSnap outputs")
    parser.add_argument("basedir", type=str, help="Location of the HBTplus output")
    parser.add_argument("hbt_nr", type=int, help="Index of the HBT output to process")
    parser.add_argument("snap_nr", type=int, help="Index of the snapshot to process")
    parser.add_argument("--snapshot-file", type=str, help="Format string for snapshot files (f-string using {snap_nr}, {file_nr})")

    args = parser.parse_args()

    sort_hbt_output(**vars(args))