#!/bin/env python

from mpi4py import MPI
comm = MPI.COMM_WORLD
comm_rank = comm.Get_rank()
comm_size = comm.Get_size()

import os
import h5py
import numpy as np

import virgo.mpi.util
import virgo.mpi.parallel_hdf5 as phdf5
import virgo.mpi.parallel_sort as psort

def log(message, quiet=False):
    if (not quiet) and (comm_rank == 0):
        print(message)

def read_snapshot(snapshot_file, snap_nr, particle_ids):
    """
    Read particle properties for the specified particle IDs.
    Returns a dict of arrays.
    """

    # Datasets to pass through from the snapshot
    passthrough_datasets = ("Coordinates", "Masses", "FOFGroupIDs")

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
        log(f"Reading snapshot particle IDs for type {ptype}", quiet=quiet)
        snapshot_ids = Mf.read(f"PartType{ptype}/ParticleIDs")
        # For each subhalo particle ID, find matching index in the snapshot (if any)
        ptr = psort.parallel_match(particle_ids, snapshot_ids, comm=comm)
        matched = (ptr>=0)
        # Loop over particle properties to pass through
        for name in passthrough_datasets:
            # Read this property from the snapshot
            if ptype == 5 and name == "Masses":
                snapshot_data = mf.read(f"PartType{ptype}/DynamicalMasses")
            else:
                snapshot_data = mf.read(f"PartType{ptype}/{name}")
            # Allocate output array, if we didn't already
            if name not in particle_data:
                shape = (len(particle_ids),)+snapshot_data.shape[1:]
                dtype = snapshot_data.dtype
                particle_data[name] = -np.ones(shape, dtype=dtype) # initialize to -1 = not found
            # Look up the value for each subhalo particle
            log(f"Looking up particle type {ptype} property {name} from snapshot", quiet=quiet)
            particle_data[name][matched,...] = psort.fetch_elements(snapshot_data, ptr[matched], comm=comm)
        # Also store the type of each matched particle
        particle_data["Type"][matched] = ptype

    # Should have matched all particles
    assert np.all(particle_data["Type"] >= 0)

    return particle_data

def read_hbt_particles(filenames, nr_local_subhalos, prop_name='SubhaloParticles'):
    """
    Read in the ParticleIDs/PotentialEnergies/BindingEnergies belonging to the subhalos 
    on this MPI rank from the specified SubSnap files. Returns a single array with
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
                nr_files = int(infile["NumberOfFiles"][0])
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
    prop_values = []
    for file_nr in range(nr_files):

        # Find range of subhalos this rank read from this file
        i1 = first_local_subhalo - first_subhalo_in_file[file_nr]
        i2 = i1 + nr_local_subhalos
        i1 = max(0, i1)
        i2 = min(subhalos_per_file[file_nr], i2)

        # Read subhalo particle IDs, if there are any in this file for this rank
        if i2 > i1:
            with h5py.File(filenames.format(file_nr=file_nr), "r") as infile:
                prop_values.append(infile[prop_name][i1:i2])

    if len(prop_values) > 0:
        # Combine arrays from different files
        prop_values = np.concatenate(prop_values)
        # Combine arrays from different subhalos
        prop_values = np.concatenate(prop_values)        
    else:
        # Some ranks may have read zero files
        prop_values = None
    prop_values = virgo.mpi.util.replace_none_with_zero_size(prop_values, comm=comm)

    # Handle case of no subhalos on any rank
    if prop_values is None:
        prop_values = np.zeros(0, dtype=int)

    return prop_values

def read_swift_cell_grid(swift_cell_file, snap_nr):
    """
    Read the SWIFT top level cell structure from a snapshot file.

    Returns the number of cells along each axis and the comoving cell size
    converted to Mpc (to match the units of the halo centres below).
    """
    from virgo.util.partial_formatter import PartialFormatter
    formatter = PartialFormatter()
    filename = formatter.format(swift_cell_file, snap_nr=snap_nr, file_nr=None)
    filename = filename.format(file_nr=0)

    mpc_in_cgs = 3.0856775814913673e24
    with h5py.File(filename, "r") as infile:
        cell_dimension = infile["Cells/Meta-data"].attrs["dimension"][:].astype(np.int64)
        cell_size = np.asarray(infile["Cells/Meta-data"].attrs["size"], dtype=np.float64)
        length_in_cgs = float(infile["Units"].attrs["Unit length in cgs (U_L)"][0])
        h = float(infile["Cosmology"].attrs["h"][0])

    # SWIFT cell sizes are comoving and expressed in SWIFT internal length units
    cell_size = cell_size * length_in_cgs / mpc_in_cgs
    return cell_dimension, cell_size, h

def compute_soap_index(cofp, nbound, length_in_mpch, swift_cell_file, snap_nr, quiet):
    """
    For each subhalo, return the index of the corresponding halo in the SOAP
    catalogue, or -1 for subhalos which SOAP does not process (unresolved
    "orphan" subhalos with Nbound == 0).

    SOAP orders halos by (SWIFT top level cell containing the halo centre, index
    of the halo in this TrackId-sorted catalogue); see the spatial_sort function
    in SOAP/SOAP/core/combine_chunks.py. This routine reproduces that ordering.

    The subhalos must already be sorted by TrackId and distributed over comm in
    the usual way. cofp is the (n_local, 3) ComovingMostBoundPosition array in HBT
    length units and length_in_mpch is the HBT Units/LengthInMpch value.
    """
    log("Computing SOAP catalogue index", quiet=quiet)

    nr_local = len(nbound)
    if comm.allreduce(nr_local) == 0:
        return -np.ones(nr_local, dtype=np.int64)

    # Read the SWIFT cell structure on rank 0 and broadcast
    if comm_rank == 0:
        cell_dimension, cell_size, h = read_swift_cell_grid(swift_cell_file, snap_nr)
    else:
        cell_dimension = cell_size = h = None
    cell_dimension, cell_size, h = comm.bcast((cell_dimension, cell_size, h))

    # SOAP's cell index hash assumes an equal number of cells along each axis
    assert cell_dimension[0] == cell_dimension[1] == cell_dimension[2], \
        "SOAP spatial sort assumes the same number of SWIFT cells on each axis"

    # Convert halo centres to comoving Mpc, matching SOAP's InputHalos/HaloCentre
    # (see cofp in SOAP/SOAP/catalogue_readers/read_hbtplus.py)
    halo_centre = np.asarray(cofp, dtype=np.float64) * (length_in_mpch / h)

    # SOAP only processes resolved subhalos
    keep = nbound > 0

    # Cell index of each halo centre (matches spatial_sort)
    cell_indices = (halo_centre // cell_size).astype(np.int64)
    if np.any(keep):
        assert np.min(cell_indices[keep]) >= 0
        for i_cell in range(3):
            assert np.max(cell_indices[keep, i_cell]) < cell_dimension[i_cell]
    cell_index = (
        cell_indices[:, 0] * cell_dimension[0] ** 2
        + cell_indices[:, 1] * cell_dimension[1]
        + cell_indices[:, 2]
    )

    # Global index of each subhalo in this TrackId-sorted catalogue. This is the
    # quantity SOAP stores as InputHalos/HaloCatalogueIndex and uses to break
    # ties between halos which fall in the same cell.
    first_local = comm.scan(nr_local) - nr_local
    catalogue_index = np.arange(nr_local, dtype=np.int64) + first_local

    # Global index among the resolved subhalos only
    nr_local_kept = int(np.sum(keep))
    first_local_kept = comm.scan(nr_local_kept) - nr_local_kept
    kept_global_index = np.arange(nr_local_kept, dtype=np.int64) + first_local_kept

    # Establish SOAP's ordering of the resolved subhalos by sorting on
    # (cell_index, catalogue_index).
    sort_key = np.zeros(
        nr_local_kept, dtype=[("cell_index", np.int64), ("catalogue_index", np.int64)]
    )
    sort_key["cell_index"] = cell_index[keep]
    sort_key["catalogue_index"] = catalogue_index[keep]
    soap_order = psort.parallel_sort(sort_key, return_index=True, comm=comm)

    # soap_order[j] is the resolved global index of the subhalo at SOAP position
    # j, so matching each resolved subhalo against soap_order gives its position.
    soap_index_kept = psort.parallel_match(kept_global_index, soap_order, comm=comm)
    assert np.all(soap_index_kept >= 0)

    # Scatter back into an array covering all subhalos, with -1 for orphans
    soap_index = -np.ones(nr_local, dtype=np.int64)
    soap_index[keep] = soap_index_kept
    return soap_index

def sort_hbt_output(basedir, snap_nr, outdir, with_particles, with_potential_energy, with_binding_energy, snapshot_file, swift_cell_file, quiet):
    """
    This reorganizes a set of HBT SubSnap files into a single file which
    contains one HDF5 dataset for each subhalo property. Subhalos are written
    in order of TrackId.

    Particle IDs in groups can be optionally copied to the output.
    """

    # Make a format string for the filenames
    filenames = f"{basedir}/{snap_nr:03d}/SubSnap_{snap_nr:03d}" + ".{file_nr}.hdf5"

    # If we're adding the SOAP index, read the HBT length unit needed to convert
    # halo centres into the same units SOAP uses
    if swift_cell_file is not None:
        if comm_rank == 0:
            with h5py.File(filenames.format(file_nr=0), "r") as infile:
                if "Units" in infile:
                    length_in_mpch = float(infile["Units/LengthInMpch"][0])
                else:
                    length_in_mpch = None
        else:
            length_in_mpch = None
        length_in_mpch = comm.bcast(length_in_mpch)
        if length_in_mpch is None:
            raise RuntimeError(
                "--swift-cell-file was given but the HBT input has no Units group, "
                "so halo centres cannot be converted to SOAP units"
            )

    # Read in the input subhalos
    log(f"Reading HBT-HERONS output for snapshot {snap_nr}", quiet=quiet)
    mf = phdf5.MultiFile(filenames, file_nr_dataset="NumberOfFiles", comm=comm)
    subhalos = mf.read("Subhalos")
    field_names = list(subhalos.dtype.fields)

    if with_particles:

        log(f"Reading particle IDs", quiet=quiet)

        # Read the particle IDs in our local subhalos
        particle_ids = read_hbt_particles(filenames, len(subhalos))
        nbound = subhalos["Nbound"]
        nr_local_particles = len(particle_ids)
        assert nr_local_particles == np.sum(nbound)

        if with_potential_energy:
            log(f"Reading potential energy", quiet=quiet)
            potential_energies = read_hbt_particles(filenames, len(subhalos), prop_name='PotentialEnergies')

        if with_binding_energy:
            log(f"Reading binding energy", quiet=quiet)
            binding_energies = read_hbt_particles(filenames, len(subhalos), prop_name='BindingEnergies')

        # Assign TrackIds to the particles
        particle_sort_key = np.repeat(subhalos["TrackId"], subhalos["Nbound"]).astype(np.int64)

        # Find maximum size of any subhalo
        if len(subhalos) > 0:
            max_subhalo_size = np.amax(subhalos["Nbound"])
        else:
            max_subhalo_size = 0
        max_subhalo_size = comm.allreduce(max_subhalo_size, op=MPI.MAX)

        # Convert trackid to sort key which includes ordering by energy
        particle_sort_key *= max_subhalo_size
        offset = 0
        for n in subhalos["Nbound"]:
            particle_sort_key[offset:offset+n] += np.arange(n, dtype=int)
            offset += n

    # Find total number of subhalos
    total_nr_subhalos = comm.allreduce(len(subhalos))

    # Convert array of structs to dict of arrays
    data = {}
    for name in field_names:
        data[name] = np.ascontiguousarray(subhalos[name])
    del subhalos

    # Establish TrackId ordering for the subhalos
    log("Sorting subhalos by TrackId", quiet=quiet)
    order = psort.parallel_sort(data["TrackId"], return_index=True, comm=comm)

    # Sort the subhalo properties by TrackId
    for name in field_names:
        if name != "TrackId":
            log(f"Reordering subhalo property: {name}", quiet=quiet)
            data[name] = psort.fetch_elements(data[name], order, comm=comm)
    del order

    # Optionally add the index of each subhalo in the corresponding SOAP catalogue
    if swift_cell_file is not None:
        data["SOAPIndex"] = compute_soap_index(
            data["ComovingMostBoundPosition"], data["Nbound"],
            length_in_mpch, swift_cell_file, snap_nr, quiet,
        )
        field_names.append("SOAPIndex")

    if with_particles:

        # Sort particle IDs too
        log("Sorting particles by TrackId", quiet=quiet)
        order = psort.parallel_sort(particle_sort_key, return_index=True, comm=comm)
        log(f"Reordering particle IDs by TrackId and energy", quiet=quiet)
        particle_ids = psort.fetch_elements(particle_ids, order, comm=comm)
        if with_potential_energy:
            log(f"Reordering potential energies by TrackId and energy", quiet=quiet)
            potential_energies = psort.fetch_elements(potential_energies, order, comm=comm)
        if with_binding_energy:
            log(f"Reordering binding energies by TrackId and energy", quiet=quiet)
            binding_energies = psort.fetch_elements(binding_energies, order, comm=comm)
        del order
        del particle_sort_key

        # Compute offset to each subhalo after sorting by TrackId
        nbound = data["Nbound"]
        nr_local_particles = sum(nbound)
        particle_offset = np.cumsum(nbound) - nbound # offset on this rank
        particle_offset += (comm.scan(nr_local_particles) - nr_local_particles) # convert to global offset

    if snapshot_file is not None:
        particle_data = read_snapshot(snapshot_file, snap_nr, particle_ids)

    # Write subhalo properties to the output file
    output_filename = f"{outdir}/OrderedSubSnap_{snap_nr:03d}.hdf5"
    log(f"Writing file: {output_filename}", quiet=quiet)
    if comm_rank == 0:
        os.makedirs(outdir, exist_ok=True)
    with h5py.File(output_filename, "w", driver="mpio", comm=comm) as outfile:

        # Create groups
        subhalo_group = outfile.create_group("Subhalos")
        if with_particles:
            particle_group = outfile.create_group("Particles")

        # Write HBT subhalo property fields
        for name in field_names:
            phdf5.collective_write(subhalo_group, name, data[name], comm)

        # Write out particle info
        if with_particles:
            phdf5.collective_write(particle_group, "ParticleIDs", particle_ids, comm)            
            phdf5.collective_write(subhalo_group, "ParticleOffset", particle_offset, comm)            
            if with_potential_energy:
                phdf5.collective_write(particle_group, "PotentialEnergies", potential_energies, comm)            
            if with_binding_energy:
                phdf5.collective_write(particle_group, "BindingEnergies", binding_energies, comm)            
            if snapshot_file is not None:
                for name, data in particle_data.items():
                    phdf5.collective_write(particle_group, name, data, comm)

    # Copy metadata from the first file
    comm.barrier()
    log("Copying metadata groups", quiet=quiet)
    if comm_rank == 0:
        input_filename = f"{basedir}/{snap_nr:03d}/SubSnap_{snap_nr:03d}" + ".0.hdf5"
        with h5py.File(input_filename, "r") as input_file, h5py.File(output_filename, "r+") as output_file:        
            for name in ("Cosmology", "Header", "Units"):
                input_file.copy(name, output_file)
            # Add some other datasets usually found in HBT SubSnap files
            output_file["NumberOfFiles"] = (1,)
            output_file["SnapshotId"] = input_file["SnapshotId"][...]
            output_file["NumberOfSubhalosInAllFiles"] = (total_nr_subhalos,)
            if swift_cell_file is not None:
                output_file["Subhalos/SOAPIndex"].attrs["Description"] = (
                    "Index of this subhalo in the corresponding SOAP catalogue, or "
                    "-1 for subhalos not included in SOAP (orphan subhalos with "
                    "Nbound=0)."
                )

    comm.barrier()

if __name__ == "__main__":

    from virgo.mpi.util import MPIArgumentParser

    parser = MPIArgumentParser(comm, description="Reorganize HBT-HERONS SubSnap outputs")
    parser.add_argument("basedir", type=str, help="Location of the HBT-HERONS output")
    parser.add_argument("snap_nr", type=int, help="Index of the snapshot to process")
    parser.add_argument("outdir",  type=str, help="Directory in which to write the output")
    parser.add_argument("--with-particles", action="store_true", help="Also copy the particle IDs to the output")
    parser.add_argument("--with-potential-energy", action="store_true", help="Also copy the particle potential energies to the output")
    parser.add_argument("--with-binding-energy", action="store_true", help="Also copy the particle binding energies to the output")
    parser.add_argument("--snapshot-file", type=str, help="Format string for snapshot files (f-string using {snap_nr}, {file_nr})")
    parser.add_argument("--swift-cell-file", type=str, help="Format string for SWIFT snapshot files (f-string using {snap_nr}, {file_nr}). If given, add a Subhalos/SOAPIndex dataset pointing to each subhalo's position in the SOAP catalogue")
    parser.add_argument("--quiet", action="store_true", help="Suppress logging")

    args = parser.parse_args()

    log(f'Running on {comm_size} ranks with the following arguments:')
    for key, value in vars(args).items():
        log(f'  {key}: {value}')

    if args.snapshot_file or args.with_potential_energy or args.with_binding_energy:
        assert args.with_particles

    sort_hbt_output(**vars(args))

    log("Done.")
