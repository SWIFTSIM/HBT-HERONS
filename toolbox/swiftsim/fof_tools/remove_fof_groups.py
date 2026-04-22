#!/bin/env python

from mpi4py import MPI

comm = MPI.COMM_WORLD
comm_rank = comm.Get_rank()
comm_size = comm.Get_size()

import os

os.environ["OPENBLAS_NUM_THREADS"] = "1"

import h5py
import numpy as np
from virgo.mpi.gather_array import gather_array
import virgo.mpi.parallel_hdf5 as phdf5


def load_group_catalogues(catalogue_filenames):
    """
    Reads the Swift FoF group catalogue and returns their IDs and particle
    sizes.

    The FOF virtual files sometimes have NumFilesPerSnapshot != 1, so we
    explicity pass a list of the filenames to read

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
    file = phdf5.MultiFile(catalogue_filenames, comm=comm)
    group_sizes = file.read(f"Groups/Sizes")
    group_ids = file.read(f"Groups/GroupIDs")
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
    file = phdf5.MultiFile(
        snapshot_path, file_nr_attr=("Header", "NumFilesPerSnapshot"), comm=comm
    )
    particle_group_ids = file.read(f"PartType{particle_type}/FOFGroupIDs_old")

    return particle_group_ids


def copy_attrs(src_obj, dst_obj):
    for key, val in src_obj.attrs.items():
        dst_obj.attrs[key] = val


def copy_object(src_obj, dst_obj, src_filename, prefix=""):
    copy_attrs(src_obj, dst_obj)
    for name, item in src_obj.items():
        if isinstance(item, h5py.Dataset):
            shape = item.shape
            dtype = item.dtype
            layout = h5py.VirtualLayout(shape=shape, dtype=dtype)
            vsource = h5py.VirtualSource(src_filename, prefix + name, shape=shape)
            layout[...] = vsource
            dst_obj.create_virtual_dataset(name, layout)
            copy_attrs(item, dst_obj[name])
        elif isinstance(item, h5py.Group):
            new_group = dst_obj.create_group(name)
            copy_object(
                item,
                new_group,
                src_filename,
                prefix + name + "/",
            )


def remove_fof_groups(
    snapshot_base_folder=None,
    output_base_folder=None,
    snapshot_basename=None,
    snapshot_is_distributed=None,
    fof_base_folder=None,
    fof_is_distributed=None,
    snap_nr=None,
    size_threshold=None,
):
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

    # Make a format string for the input/output snapshot filenames
    snap_filepath = (
        f"{snapshot_base_folder}/{snapshot_basename}_{snap_nr:04d}/{snapshot_basename}_{snap_nr:04d}"
        + (".{file_nr}.hdf5" if snapshot_is_distributed else ".hdf5")
    )
    output_filepath = (
        f"{output_base_folder}/{snapshot_basename}_{snap_nr:04d}/{snapshot_basename}_{snap_nr:04d}"
        + (".{file_nr}.hdf5" if snapshot_is_distributed else ".hdf5")
    )
    if comm_rank == 0:
        with h5py.File(snap_filepath.format(file_nr=0), "r") as file:
            nr_files = file["Header"].attrs["NumFilesPerSnapshot"][0]
        output_dirname = os.path.dirname(output_filepath)
        if not os.path.exists(output_dirname):
            os.makedirs(output_dirname, exist_ok=True)
        if os.path.exists(output_filepath.format(file_nr=0)):
            print('Output file already exist')
            comm.Abort()
        if nr_files < comm_size:
            print("Can't run with more ranks than snapshot chunk files")
            comm.Abort()
    else:
        nr_files = None
    nr_files = comm.bcast(nr_files)

    # The FOF virtual files sometimes have NumFilesPerSnapshot != 1, so we
    # explicity create a list of the filenames to pass to MultiFile
    if fof_is_distributed:
        catalogue_basename = f"{fof_base_folder}/fof_output_{snap_nr:04d}/fof_output_{snap_nr:04d}.{{file_nr}}.hdf5"
        if comm_rank == 0:
            with h5py.File(catalogue_basename.format(file_nr=0), "r") as file:
                nr_fof_files = file["Header"].attrs["NumFilesPerSnapshot"][0]
        else:
            nr_fof_files = None
        nr_fof_files = comm.bcast(nr_fof_files)
        catalogue_filenames = [
            catalogue_basename.format(file_nr=i) for i in range(nr_fof_files)
        ]
    else:
        catalogue_filenames = [
            f"{fof_base_folder}/fof_output_{snap_nr:04d}/fof_output_{snap_nr:04d}.hdf5"
        ]
        if comm_rank == 0:
            with h5py.File(catalogue_filenames[0], "r") as file:
                if (file["Header"].attrs["Virtual"][0] == 1) and (comm_size > 1):
                    print("Can't read a virtual FOF file with multiple ranks")
                    comm.Abort()
        comm.barrier()

    if comm_rank == 0:
        print("Creating output files")
    # Assign files to ranks
    files_per_rank = np.zeros(comm_size, dtype=int)
    files_per_rank[:] = nr_files // comm_size
    remainder = nr_files % comm_size
    if remainder > 0:
        step = max(nr_files // (remainder + 1), 1)
        for i in range(remainder):
            files_per_rank[(i * step) % comm_size] += 1
    first_file = np.cumsum(files_per_rank) - files_per_rank
    assert sum(files_per_rank) == nr_files

    # Create output files
    for i_file in range(
        first_file[comm_rank], first_file[comm_rank] + files_per_rank[comm_rank]
    ):
        src_filename = snap_filepath.format(file_nr=i_file)
        dst_filename = output_filepath.format(file_nr=i_file)
        with (
            h5py.File(src_filename, "r") as src_file,
            h5py.File(dst_filename, "w") as dst_file,
        ):
            rel_filename = os.path.relpath(src_filename, os.path.dirname(dst_filename))
            copy_object(src_file, dst_file, rel_filename)

    if comm_rank == 0:
        print(f"Reading FoF catalogue for snapshot number {snap_nr}")
        print(f"Opening file: {catalogue_filenames}")
        with h5py.File(catalogue_filenames[0], "r") as file:
            header_total_nr_groups = file["Header"].attrs["NumGroups_Total"][0]
    comm.barrier()

    # Read FOF group ids to remove from catalogue
    local_group_ids, local_group_sizes = load_group_catalogues(catalogue_filenames)

    # Get the group null id used in this simulation
    if comm_rank == 0:
        with h5py.File(output_filepath.format(snap_nr=snap_nr, file_nr=0)) as file:
            null_group_id = int(file["Parameters"].attrs["FOF:group_id_default"])
    else:
        null_group_id = None
    null_group_id = comm.bcast(null_group_id, root=0)
    assert null_group_id is not None
    comm.barrier()

    # Find total number of haloes that require fixing and the total
    total_nr_groups_to_remove = comm.allreduce(
        (local_group_sizes < size_threshold).sum()
    )
    total_nr_groups = comm.allreduce(len(local_group_ids))
    if total_nr_groups_to_remove == 0:
        if comm_rank == 0:
            print("No FOF groups need to be removed. Exiting now.")
        return
    if comm_rank == 0:
        print(
            f"We will remove {total_nr_groups_to_remove} FOF groups from the snapshots. This is {total_nr_groups_to_remove / total_nr_groups * 100:.3f}% of the total."
        )
    comm.barrier()

    # Create the array "remove_group"
    # This can indexed using the FOFGroupIDs value from a particle
    # and indicates whether that FOF should be kept
    global_group_ids = gather_array(local_group_ids)
    global_group_sizes = gather_array(local_group_sizes)
    if comm_rank == 0:
        max_group_id = np.max(global_group_ids)
        assert header_total_nr_groups == total_nr_groups
        assert max_group_id == global_group_ids.shape[0]
        assert np.unique(global_group_ids).shape[0] == global_group_ids.shape[0]
        remove_group = np.zeros(max_group_id + 1, dtype=bool)
        remove_group[global_group_ids] = global_group_sizes < size_threshold
    else:
        remove_group = None
    remove_group = comm.bcast(remove_group)

    # Before changing anything, make sure that we have renamed the original FOF dataset
    if comm_rank == 0:
        print("Renaming original FOF dataset to FOFGroupIDs_old.")
    for i_file in range(
        first_file[comm_rank], first_file[comm_rank] + files_per_rank[comm_rank]
    ):
        with h5py.File(output_filepath.format(file_nr=i_file), "a") as file:
            for particle_type in [0, 1, 4, 5]:
                # Check we will not overwrite anything by accident
                assert f"PartType{particle_type}/FOFGroupIDs_old" not in file

                group = file[f"PartType{particle_type}"]
                group.move("FOFGroupIDs", "FOFGroupIDs_old")
    comm.barrier()

    # Load the fof group ids of each particle type.
    file = phdf5.MultiFile(
        output_filepath,
        file_nr_attr=("Header", "NumFilesPerSnapshot"),
        comm=comm,
    )
    if comm_rank == 0:
        print(f"Reading snapshot catalogue for snapshot number {snap_nr}")
    for particle_type in [0, 1, 4, 5]:
        if comm_rank == 0:
            print(f"Doing PartType{particle_type}")

        # Read memberships
        data = {}
        data[f"FOFGroupIDs"] = read_particle_groups(output_filepath, particle_type)
        comm.barrier()

        # Find those we need to remove
        mask = data[f"FOFGroupIDs"] != null_group_id
        mask[mask] = remove_group[data["FOFGroupIDs"][mask]]
        data[f"FOFGroupIDs"][mask] = null_group_id

        # Check how many particles of this type we have per file
        number_particles_per_file = file.get_elements_per_file(
            f"PartType{particle_type}/FOFGroupIDs_old"
        )
        comm.barrier()

        # Save
        file.write(
            data,
            number_particles_per_file,
            output_filepath,
            "a",
            group=f"PartType{particle_type}",
        )

        comm.barrier()
    comm.barrier()

    if comm_rank == 0:
        print("Copying attributes from old to new FOF dataset.")
    for i_file in range(
        first_file[comm_rank], first_file[comm_rank] + files_per_rank[comm_rank]
    ):
        with h5py.File(output_filepath.format(file_nr=i_file), "a") as file:
            for particle_type in [0, 1, 4, 5]:
                attributes = dict(
                    file[f"PartType{particle_type}/FOFGroupIDs_old"].attrs
                )
                for attr, value in attributes.items():
                    file[f"PartType{particle_type}/FOFGroupIDs"].attrs[attr] = value

                # Update the description of the old FOF group ids, to prevent confusion
                file[f"PartType{particle_type}/FOFGroupIDs_old"].attrs[
                    "Description"
                ] = "Old values of the Friends-Of-Friends ID membership of particles. Use the new dataset instead."
    comm.barrier()

    if comm_rank == 0:
        print("Done.")


if __name__ == "__main__":

    from virgo.mpi.util import MPIArgumentParser

    parser = MPIArgumentParser(
        comm, description="Remove FoF groups whose size is below the chosen threshold"
    )
    parser.add_argument(
        "--snapshot_base_folder",
        type=str,
        required=True,
        help="Base path of the input particle snapshots.",
    )
    parser.add_argument(
        "--snapshot_basename",
        type=str,
        required=True,
        help="Base name used for snapshots.",
    )
    parser.add_argument(
        "--snapshot_is_distributed",
        action="store_true",
        help="Whether multiple subfiles exist per snapshot.",
    )
    parser.add_argument(
        "--fof_base_folder",
        type=str,
        required=True,
        help="Base path where FoF catalogues are saved.",
    )
    parser.add_argument(
        "--fof_is_distributed",
        action="store_true",
        help="Whether multiple subfiles exist per FOF group catalogue.",
    )
    parser.add_argument(
        "--snap_nr", type=int, required=True, help="Snapshot number to processs."
    )
    parser.add_argument(
        "--size_threshold",
        type=int,
        required=True,
        help="Minimum size threshold for a FoF group to not be removed.",
    )
    parser.add_argument(
        "--output_base_folder",
        type=str,
        required=True,
        help="Base path to save the edited particle snapshots.",
    )

    args = parser.parse_args()
    if comm_rank == 0:
        print("Running with the following arguments:")
        for k, v in vars(args).items():
            print(f"  {k}: {v}")

    remove_fof_groups(**vars(args))
