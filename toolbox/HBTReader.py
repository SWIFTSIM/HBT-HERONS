"""
Class to read HBT-HERONS outputs.
NOTE: only use for small simulations, as this method is very slow to
retrieve the properties of all subhaloes. Suggested alternative is to
instead run the script "./SortCatalogues.py", which will sort catalogues in
ascending TrackId order.

To use it, initialize the reader by pointing it to the directory containing the
subhalo catalogues, e.g.

  from HBTReader import HBTReader
  catalogues = HBTReader("<PATH_TO_SUBHALO_CATALOGUES>")

  # Load all subhaloes and their properties in the last analysed snapshot.
  subhaloes = reader.LoadSubhalos()

  # Loads the Nbound of all subhaloes in the last analysed snapshot.
  Nbound = reader.LoadSubhalos("Nbound")

  # Loads the Nbound of the eight subhalo entry in the last analysed snapshot.
  subhalo_Nbound = reader.LoadSubhalos(subindex = 8)

  # Gets the property evolution of TrackId = 2
  subhalo_evolution = reader.GetTrackEvolution(2)
"""

import h5py
import numpy as np
from glob import glob
from numpy.lib.recfunctions import append_fields

def get_hbt_snapnum(snapname, is_sorted):
    """
    Get the snapshot number from the name of a HBT-HERONS catalogue file.

    Parameters
    ==========
    snapname: str
        Name of the HBT-HERONS catalogue file.
    is_sorted: bool
        Whether the catalogue is sorted.

    Returns
    =======
    int
        Snapshot number.
    """

    if is_sorted:
        return int(snapname.rsplit('OrderedSubSnap_')[1].split('.')[0])
    else:
        return int(snapname.format(filetype="Sub", subfile_nr=0).rsplit('SubSnap_')[1].split('.')[0])

class HBTReader:
    """
    Class to read HBT-HERONS catalogues.
    """
    mass_units = None
    length_units = None
    velocity_units = None

    def __init__(self, base_dir, sorted_catalogues = False):
        """
        Initialize HBTReader to read data from base_dir where all the outputs
        are saved.

        Parameters
        ==========
        base_dir: str
            Base directory of where the HBT-HERONS catalogues are saved.
        sorted_catalogues: bool, opt
            Whether the catalogues have been merged into a single file and
            sorted in ascending TrackId (see ./catalogue_cleanup/SortCatalogues.py).
            It defaults to False.
        """

        self.base_dir = base_dir
        self.__sorted_catalogues = sorted_catalogues
        self.__file_list = self.GetFileList()

    def GetFileList(self):
        """
        Get the paths to the catalogue files, sorted in ascending output number.

        Returns
        =======
        list of str
            List of paths to the catalogue files, sorted in ascending output
            number.
        """

        if self.__sorted_catalogues:
            file_list = sorted(glob(self.base_dir + "/OrderedSubSnap_*.hdf5"), key=lambda x: get_hbt_snapnum(x, self.__sorted_catalogues))
        else:
            # Create f-formatted strings because we may have several subfiles per snapshot.
            file_list = sorted(glob(self.base_dir + "/*/SubSnap_*.0.hdf5"), key=lambda x: get_hbt_snapnum(x, self.__sorted_catalogues))
            file_list = [path.replace(".0.hdf5",".{subfile_nr}.hdf5").replace("SubSnap","{filetype}Snap") for path in file_list]

        # Did we find any files?
        assert len(file_list) > 0, "No catalogue files found in the provided directory."
        self.SnapshotIdList = np.array([get_hbt_snapnum(path, self.__sorted_catalogues) for path in file_list])

        return file_list

    def GetFileName(self, snap_nr, subfile_nr=0, filetype='Sub'):
        """
        Returns the path to a catalogue file.

        Parameters
        ==========
        snap_nr: int
            The snapshot number of the catalogue we are interested in.
        subfile_nr: int, opt
            The subfile_nr we are interested in. Defaults to 0.
        filetype: str, opt
            Whether to load the bound subhalo information ('Sub') or the source
            subhalo information ('Src'). Defaults to 'Sub'.

        Returns
        =======
        str
            Path to the subfile.
        """

        index = np.where(self.SnapshotIdList == snap_nr)[0]
        if len(index) == 0:
            raise ValueError(f"The requested snapshot number ({snap_nr}) is not available. Possible values:\n {self.SnapshotIdList}")

        return self.__file_list[index[0]].format(subfile_nr = subfile_nr, filetype = filetype)

    def LoadNestedSubhalos(self, snap_nr=None):
        """
        Returns the list of nested subhalo indices for each subhalo.

        Parameters
        ==========
        snap_nr : int, opt
            Snapshot number of the catalogue we are interested in. It defaults
            to the last snapshot with currently available catalogues.
        """
        nests = []
        for i in range(max(self.nfiles, 1)):
            with h5py.File(self.GetFileName(snap_nr, i), 'r') as subfile:
                nests.extend(subfile['NestedSubhalos'][...])
        return np.array(nests)

    def LoadSubhalos(self, snap_nr=None, subhalo_index=None, property_selection=None, show_progress=False):
        """
        Load subhalos from a single snapshot, with the option to load a subset
        of properties and subhaloes.

        Parameters
        ==========
        snap_nr: int, opt
            The snapshot number we are interested in. It defaults to the
            latest snapshot with available catalogues.
        subhalo_index: int, opt
            If specified, only load the subhalo in the specified entry. Note
            that this does NOT correspond to the TrackId of the subhalo, as
            the catalogues are not sorted in TrackId. If not provided, all
            subhaloes will be loaded.
        property_selection: tuple or list of strings, opt
            If specified, only load the specified properties. If not provided,
            all subhalo properties will be loaded.
        show_progess: bool, opt
            Prints progress bar indicating how many subfiles have been loaded.
            Defaults to false.

        Returns
        =======
        subhalos: np.ndarray
            Specified properties for the requested subhaloes at the snapshot
            of interest. The array contains multiple dtypes, with each
            corresponding to a different subhalo property. They can be accessed
            via indexing, e.g. subhalos["PROPERTY_NAME"]
        """

        if snap_nr is None:
            snap_nr = self.SnapshotIdList.max()

        load_single_subhalo = subhalo_index is not None
        if (load_single_subhalo):
            if not isinstance(subhalo_index, (np.integer, int)):
                raise TypeError("Parameter subhalo_index is not of the required type (int).")
            if subhalo_index >= self.GetNumberOfSubhalos(snap_nr):
                raise ValueError(f"Selected subhalo entry ({subhalo_index}) is larger than the number of existing subhaloes ({self.GetNumberOfSubhalos(snap_nr)})")

        # Handle defaults, and list inputs
        if property_selection is None:
            property_selection = np.s_[:]
        else:
            if type(property_selection) is list:
                property_selection = tuple(property_selection)

        # Determine number of files for requested output
        with h5py.File(self.GetFileName(snap_nr), 'r') as subfile:
            number_subfiles = subfile['NumberOfFiles'][0]

        offset = 0
        subhalos = []
        for subfile_nr in range(number_subfiles):
            if show_progress:
                print(".", end="")

            with h5py.File(self.GetFileName(snap_nr, subfile_nr), 'r') as subfile:
                nsub = subfile['Subhalos'].shape[0]

                # Nothing to load
                if nsub == 0:
                    continue

                # Keep iterating over subfiles until we find our target entry.
                if load_single_subhalo:
                    if (offset+nsub > subhalo_index):
                        subhalos.append(subfile['Subhalos'][property_selection][subhalo_index-offset])
                        break
                    else:
                        offset += nsub
                        continue

                # Load everything
                subhalos.append(subfile['Subhalos'][property_selection])

        if len(subhalos):
            subhalos = np.hstack(subhalos)
        else:
            subhalos = np.array(subhalos)

        if show_progress:
            print()

        return subhalos

    def GetNumberOfSubhalos(self, snap_nr=None):
        """
        Returns the total number of subhaloes in a given catalogue output.

        Parameters
        ==========
        snap_nr : int, opt
            Snapshot number of the catalogue we are interested in. It defaults
            to the last snapshot with currently available catalogues.

        Returns
        =======
        int
            The total number of subhaloes (resolved + unresolved) that exist
            at the current snapshot.
        """

        if snap_nr is None:
            snap_nr = self.SnapshotIdList.max()

        with h5py.File(self.GetFileName(snap_nr, 0), 'r') as file:
            return file["NumberOfSubhalosInAllFiles"][0]

    def LoadParticleIDs(self, TrackId=None, snap_nr=None, filetype='Sub'):
        """
        Load the particles that are bound or are part of the source subhalo
        for the specified subhalo.

        Parameters
        ==========
        TrackId: int, opt
            The TrackId of the subhalo whose particles we are interested in. If
            not specified, this function will return the particle IDs of all
            TrackIds in the catalogue.
        snap_nr : int, opt
            Snapshot number of the catalogue we are interested in. It defaults
            to the last snapshot with currently available catalogues.
        file_type: str, opt
            Whether to load particles bound ('Sub') or associated ('Src') of
            the subhalo.

        Returns
        =======
        subhalo_particles: np.ndarray
            IDs of the particles that are either bound or part of the source
            subhalo of all subhaloes or a specific subhalo. They are ordered
            in the same manner as the subhalo arrays, and bound particles are
            ordered in binding energy.
        """

        if snap_nr is None:
            snap_nr = self.SnapshotIdList.max()

        if filetype not in ["Sub","Src"]:
            raise ValueError(f"Requested filetype ({filetype}) is not valid. Only \"Sub\" or \"Src\" are accepted.")
        key_to_load = "SubhaloParticles" if filetype == "Sub" else "SrchaloParticles"

        # Find the index location of the subhalo entry if we requested a given
        # TrackId. If we do not find it, it does not exist at this point.
        load_single_subhalo = TrackId is not None
        if load_single_subhalo:

            if not isinstance(TrackId, (np.integer, int)):
                raise TypeError("Parameter TrackId is not of the required type (int).")

            subhalo_index = np.flatnonzero(self.LoadSubhalos(snap_nr, property_selection='TrackId') == TrackId)
            if len(subhalo_index) == 0:
                raise LookupError("The requested TrackId does not exist in the snapshot of interest.")
            subhalo_index = subhalo_index[0]

        else:
            subhalo_index = None

        # Determine number of files for requested output, only present in SubSnap
        # files.
        with h5py.File(self.GetFileName(snap_nr), 'r') as subfile:
            number_subfiles = subfile['NumberOfFiles'][0]

        offset = 0
        subhalo_particles = []
        for subfile_nr in range(number_subfiles):
            with h5py.File(self.GetFileName(snap_nr, subfile_nr, filetype), 'r') as subfile:
                nsub = subfile['Subhalos'].shape[0] if filetype == "Sub" else subfile[key_to_load].shape[0]

                # Nothing to load
                if nsub == 0:
                    continue

                # Keep iterating over subfiles until we find our target entry.
                if load_single_subhalo:
                    if (offset+nsub > subhalo_index):
                        subhalo_particles.append(subfile[key_to_load][subhalo_index-offset])
                        break
                    else:
                        offset += nsub
                        continue

                # Load everything
                subhalo_particles.append(subfile['SubhaloParticles'][:])

        subhalo_particles = np.hstack(subhalo_particles)

        return subhalo_particles

    def GetTrackSnapshot(self, TrackId, snap_nr=None, fields=None):
        """
        Load the properties of a given TrackId in the specified snapshot.

        Parameters
        ==========
        TrackId: int
            The TrackId of the subhalo whose properties we are interested in.
        snap_nr: int
            The snapshot we are interested in.
        fields: tuple or list of properties, opt
            The properties we are interested in loading. If not defined, we load
            everything.
        """

        # Find the  index location of the requested TrackId. If we do not find
        # it, it does not exist at this point.
        subhalo_index = np.flatnonzero(self.LoadSubhalos(snap_nr, property_selection='TrackId') == TrackId)
        if len(subhalo_index) == 0:
            raise LookupError("The requested TrackId does not exist in the snapshot of interest.")
        subhalo_index = subhalo_index[0]

        subhalo = self.LoadSubhalos(snap_nr, subhalo_index=subhalo_index, property_selection=fields)

        return subhalo

    def GetTrackEvolution(self, TrackId, fields=None):
        """
        Load the entire evolution of a given TrackId, from when it was first
        resolved until the latest snapshot with available catalogues (even if it
        the subhalo is unresolved by that time).

        Parameters
        ==========
        TrackId: int
            The TrackId of the subhalo whose evolution we are interested in.
        fields: tuple or list of properties, opt
            The properties we are interested in loading. If not defined, we load
            everything.
        """

        # Only load information from when the subhalo was resolved. We try using
        # the old and new dataset name for backwards compatibility.
        try:
            SnapshotOfBirth = self.GetTrackSnapshot(TrackId)['SnapshotOfBirth']
        except:
            SnapshotOfBirth = self.GetTrackSnapshot(TrackId)['SnapshotIndexOfBirth']
        snapshots_to_load = self.SnapshotIdList[self.SnapshotIdList >= SnapshotOfBirth]

        track = np.array(
            [self.GetTrackSnapshot(TrackId, snap_nr, fields=fields)
             for snap_nr in snapshots_to_load])

        scales = np.array([self.GetScaleFactor(snap_nr) for snap_nr in snapshots_to_load])

        return append_fields(
            track, ['Snapshot', 'ScaleFactor'], [snapshots_to_load, scales], usemask=False)

    def GetScaleFactor(self, snap_nr):
        """
        Returns the scale factor of a given catalogue output.

        Parameters
        ==========
        snap_nr : int, opt
            Snapshot number of the catalogue we are interested in.

        Returns
        =======
        float
            The scale factor of the snapshot.
        """
        return h5py.File(self.GetFileName(snap_nr), 'r')['Cosmology/ScaleFactor'][0]

    def GetScaleFactorDict(self):
        """
        Returns a dictionary mapping each snapshot number to its corresponding
        scale factor.

        Returns
        =======
        dict of (int, float)
            Mapping between snapshot number and scale factor.
        """
        return dict([(snap_nr, self.GetScaleFactor(snap_nr)) for snap_nr in self.SnapshotIdList])

    def GetMassUnits_Msunh(self):
        """
        Returns the mass units of the catalogue outputs in Msun/h.

        Returns
        =======
        float
            Mass units of the catalogue outputs, in Msun/h.
        """
        if self.mass_units is None:
            self.mass_units = h5py.File(self.GetFileName(self.SnapshotIdList[0]), 'r')['Units/MassInMsunh'][0]
        return self.mass_units

    def GetLengthUnits_Mpch(self):
        """
        Returns the length units of the catalogue outputs.

        Returns
        =======
        float
            Length units of the catalogue outputs, in Mpc/h.
        """
        if self.length_units is None:
            self.length_units = h5py.File(self.GetFileName(self.SnapshotIdList[0]), 'r')['Units/LengthInMpch'][0]
        return self.length_units

    def GetVelocityUnits_kms(self):
        """
        Returns the velocity units of the catalogue outputs in km/s.

        Returns
        =======
        float
            Velocity units of the catalogue outputs, in km/s.
        """
        if self.velocity_units is None:
            self.velocity_units = h5py.File(self.GetFileName(self.SnapshotIdList[0]), 'r')['Units/VelInKmS'][0]
        return self.velocity_units