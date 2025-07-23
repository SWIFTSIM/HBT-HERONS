#!/bin/env python

from mpi4py import MPI
comm = MPI.COMM_WORLD
comm_rank = comm.Get_rank()
comm_size = comm.Get_size()

import os
os.environ["OPENBLAS_NUM_THREADS"] = "1"

import h5py
import numpy as np
import virgo.mpi.util
import virgo.mpi.gather_array as ga
import virgo.mpi.parallel_hdf5 as phdf5
import virgo.mpi.parallel_sort as psort

def load_group_catalogues(catalogue_path):
    """
    Reads the Swift FoF group catalogue and returns their IDs and particle
    sizes.

    Parameters
    ----------
    catalogue_path: str
        Path to the FoF group catalogues.

    Returns
    -------
    group_ids: np.ndarray
        Unique ID for each FoF group.
    group_sizes: np.ndarray
        Total particle size for each FoF group.
    """
    file =  phdf5.MultiFile(catalogue_path, file_nr_attr=("Header","NumFilesPerSnapshot"), comm=comm)
    group_sizes = file.read(f"Groups/Sizes")
    group_ids = file.read(f"Groups/GroupIDs")
    del file

    return group_ids, group_sizes

def read_particle_groups(snapshot_path, particle_type):
    """
    Reads the assigned FoF IDs ids of particles of a given particle type.

    Parameters
    ----------
    snapshot_path: str
        Path to the Swift particle snapshot.
    particle_type: int
        Particle type to read.

    Returns
    -------
    particle_group_ids: np.ndarray
        Host FoF group ID for each particle of the requested type.
    """
    file =  phdf5.MultiFile(snapshot_path, file_nr_attr=("Header","NumFilesPerSnapshot"), comm=comm)
    particle_group_ids = file.read(f"PartType{particle_type}/FOFGroupIDs_old")
    del file

    return particle_group_ids

def remove_fof_groups(snapshot_base_folder, snapshot_basename, snapshot_is_distributed, fof_base_folder, fof_is_distributed, snap_nr, size_threshold):
    """
    Makes all particles that belong to FoF groups whose size is
    lower than the provided threshold hostless.

    Parameters
    ----------
    snapshot_base_folder: str
        Base path where Swift particle snapshots are saved.
    snapshot_base_name: str
        Base name of Swift particle snapshots.
    snapshot_is_distributed: bool
        Whether multiple subfiles exist per snapshot (True) or not (False).
    fof_base_folder: str
        Base path where Swift FoF group catalogues are saved.
    fof_is_distributed: bool
        Whether multiple subfiles exist per fof group catalogue (True) or not (False).
    snap_nr: int
        Snapshot number to process.
    size_threshold: int
        Minimum size of FoF groups to remain in the particle data. Any FoF group
        whose size is below this threshold will be removed.
    """

    # Make a format string for the filenames
    snapshot_filenames = f"{snapshot_base_folder}/{snapshot_basename}_{snap_nr:04d}/{snapshot_basename}_{snap_nr:04d}" + (".{file_nr}.hdf5" if snapshot_is_distributed else ".hdf5")
    catalogue_filenames = f"{fof_base_folder}/fof_output_{snap_nr:04d}/fof_output_{snap_nr:04d}" + (".{file_nr}.hdf5" if fof_is_distributed else ".hdf5")

    if comm_rank == 0:
        print(f"Reading FoF catalogue for snapshot number {snap_nr}")
        print(f"Opening file: {catalogue_filenames}")

    # NOTE: we can compute the FoF group sizes and ids directly from the particle datai instead of calling
    # load_group_catalogues
    # Read FOF group ids to remove from catalogue
    local_group_ids, local_group_sizes = load_group_catalogues(catalogue_filenames)

    # Find total number of haloes that require fixing and the total
    total_nr_groups_to_remove = comm.allreduce(local_group_ids[local_group_sizes < size_threshold].sum())
    total_nr_groups = comm.allreduce(len(local_group_ids))

    if(total_nr_groups_to_remove == 0):
        if(comm_rank == 0):
            print("No FOF groups need to be removed. Exiting now.")
        return

    if comm_rank == 0:
        print(f"We will remove {total_nr_groups_to_remove} FOF groups from the snapshots. This is {total_nr_groups_to_remove / total_nr_groups_present * 100:.3f}% of the total.")

    # Before changing anything, make sure that we have renamed the original FOF dataset
    if comm_rank == 0:
        print("Renaming original FOF dataset to FOFGroupIDs_old.")

        # Get number of subfiles
        with h5py.File(snapshot_filenames.format(file_nr = 0)) as file:
            nr_subfiles = file["Header"].attrs["NumFilesPerSnapshot"][0]

        for subfile in range(nr_subfiles):
            with h5py.File(snapshot_filenames.format(file_nr = subfile),'a') as file:
                for particle_type in [0,1,4,5]:
                    # Check we will not overwrite anything by accident
                    assert(f'PartType{particle_type}/FOFGroupIDs_old' not in file)

                    group = file[f'PartType{particle_type}']
                    group.move('FOFGroupIDs','FOFGroupIDs_old')
    comm.barrier()

    # Get the group null id used in this simulation
    if comm_rank == 0:
        with h5py.File(snapshot_filenames.format(snap_nr=snap_nr,file_nr=0)) as file:
            null_group_id = int(file['Parameters'].attrs['FOF:group_id_default'])
    else:
        null_group_id = None

    null_group_id = comm.bcast(null_group_id, root=0)
    assert(null_group_id is not None)
    comm.barrier()

    # Sort out the sizes and the groups
    order = psort.parallel_sort(local_group_ids, return_index=True, comm=comm)
    sorted_group_ids   = psort.fetch_elements(local_group_ids, order, comm=comm)
    sorted_group_sizes = psort.fetch_elements(local_group_sizes, order, comm=comm)

    # Load the fof group ids of each particle type.
    if comm_rank == 0:
        print(f"Reading snapshot catalogue for snapshot number {snap_nr}")

    for particle_type in [0,1,4,5]:
        if comm_rank == 0:
            print(f"Doing PartType{particle_type}")

        # Read memberships
        data = {}
        data [f"FOFGroupIDs"] = read_particle_groups(snapshot_filenames, particle_type)
        comm.barrier()

        # Find those we need to remove
        index_to_reset = np.isin(data[f"FOFGroupIDs"],global_group_ids_to_remove)
        data[f"FOFGroupIDs"][index_to_reset] = null_group_id

        # Check how many particles of this type we have per file
        file = phdf5.MultiFile(snapshot_filenames, file_nr_attr=("Header","NumFilesPerSnapshot"), comm=comm)
        number_particles_per_file = file.get_elements_per_file(f"PartType{particle_type}/FOFGroupIDs_old")
        comm.barrier()

        # Save
        file.write(data, number_particles_per_file,snapshot_filenames,'a', group=f"PartType{particle_type}")

        comm.barrier()
        # del file

    comm.barrier()

    if comm_rank == 0:
        print("Copying attributes from old to new FOF dataset.")

        # Get number of subfiles
        with h5py.File(snapshot_filenames.format(file_nr = 0)) as file:
            nr_subfiles = file["Header"].attrs["NumFilesPerSnapshot"][0]

        for subfile in range(nr_subfiles):
            with h5py.File(snapshot_filenames.format(file_nr = subfile),'a') as file:
                for particle_type in [0,1,4,5]:
                    attributes = dict(file[f"PartType{particle_type}/FOFGroupIDs_old"].attrs)
                    for attr, value in attributes.items():
                        file[f"PartType{particle_type}/FOFGroupIDs"].attrs[attr] =  value

                    # Update the description of the old FOF group ids, to prevent confusion
                    file[f"PartType{particle_type}/FOFGroupIDs_old"].attrs['Description'] = "Old values of the Friends-Of-Friends ID membership of particles. Use the new dataset instead."

    comm.barrier()
    if comm_rank == 0:
        print("Done.")

if __name__ == "__main__":

    from virgo.mpi.util import MPIArgumentParser

    parser = MPIArgumentParser(comm, description="Remove FoF groups whose size is below the chosen threshold")
    parser.add_argument("snapshot_base_folder", type=str, help="Base path where particle snapshots are saved.")
    parser.add_argument("snapshot_base_name", type=str, help="Base name used for snapshots.")
    parser.add_argument("snapshot_is_distributed", type=bool, help="Whether multiple subfiles exist per snapshot.")
    parser.add_argument("fof_base_folder", type=str, help="Base path where FoF catalogues are saved.")
    parser.add_argument("fof_is_distributed", type=bool, help="Whether multiple subfiles exist per FOF group catalogue.")
    parser.add_argument("snap_nr", type=int, help="Snapshot number to processs.")
    parser.add_argument("size_threshold",  type=int, help="Minimum size threshold for a FoF group to not be removed.")

    args = parser.parse_args()

    update_fof_group_memberships(**vars(args))