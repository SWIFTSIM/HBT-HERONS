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

def generate_custom_array_dtypes(h5py_group, requested_properties):
    """
    Returns a list of tuples that is used to create custom dtype arrays.

    Returns
    =======
    list of tuples
        A list where each tuple specifies the name, dtype and shape of a custom
        numpy dtype.
    """

    # Build an array with custom dtypes for the requested properties.
    array_dtypes = []
    for property in requested_properties:
        if h5py_group[f"{property}"].ndim == 1:
            array_dtypes.append((property, h5py_group[f"{property}"].dtype))
        else:
            array_dtypes.append((property, h5py_group[f"{property}"].dtype, h5py_group[f"{property}"].shape[1]))

    return array_dtypes

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

    #===========================================================================
    # File path functions.
    #===========================================================================
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

        # For raw catalogues we will have a certain number of files per snapshot.
        if not self.__sorted_catalogues:
            with h5py.File(file_list[0].format(subfile_nr=0, filetype="Sub"), 'r') as subfile:
                self.nfiles = subfile["NumberOfFiles"][0]

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

        index = np.flatnonzero(self.SnapshotIdList == snap_nr)

        if len(index) == 0:
            raise ValueError(f"The requested snapshot number ({snap_nr}) is not available. Possible values:\n {self.SnapshotIdList}")

        return self.__file_list[index[0]].format(subfile_nr = subfile_nr, filetype = filetype)

    #===========================================================================
    # Data loading functions.
    #===========================================================================
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

    def LoadSubhalos(self, snap_nr=None, subhalo_index=None, property_selection=None, show_progress=False):
        """
        Load subhalos from a single snapshot, with the option to load a subset
        of properties and subhaloes.

        Parameters
        ==========
        snap_nr: int, opt
            The snapshot number we are interested in. It defaults to the
            latest snapshot with available catalogues.
        subhalo_index: int or np.ndarray, opt
            If specified, only load the subhalo(es) in the specified entries,
            which should be in ascending value order. For unsorted catalogues,
            this does not correspond to the TrackId of the subhalo, but it does
            so for the sorted version. If not provided, all subhaloes will be loaded.
        property_selection: tuple or list of strings, opt
            If specified, only load the specified properties. If not provided,
            all subhalo properties will be loaded.
        show_progess: bool, opt
            For unsorted catalogues, enable the printing of a progress bar
            indicating how many subfiles have been loaded. Defaults to False.

        Returns
        =======
        subhalos: np.ndarray
            Specified properties for the requested subhaloes at the snapshot
            of interest. The array contains multiple dtypes, with each
            corresponding to a different subhalo property. They can be accessed
            via indexing, e.g. subhalos["PROPERTY_NAME"]
        """

        if self.__sorted_catalogues:
            return self.__LoadSubhalos_SortedCatalogues(snap_nr, subhalo_index, property_selection)
        else:
            return self.__LoadSubhalos_UnsortedCatalogues(snap_nr, subhalo_index, property_selection, show_progress)

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

        if self.__sorted_catalogues and filetype == "Src":
            raise ValueError("Source subhalo information is not available in the sorted catalogues.")

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

            subhalo_index = np.flatnonzero(self.LoadSubhalos(snap_nr, property_selection=['TrackId'])["TrackId"] == TrackId)
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

    #===========================================================================
    # Functions to load properties of specified TrackIds.
    #===========================================================================
    def GetTrackSnapshot(self, TrackId, snap_nr, fields=None):
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

        if self.__sorted_catalogues:
            return self.LoadSubhalos(snap_nr, TrackId, fields)
        else:
            # Find the index location of the requested TrackId. If we do not find
            # it, it does not exist at this point.
            subhalo_index = np.flatnonzero(self.LoadSubhalos(snap_nr, property_selection=['TrackId'])["TrackId"] == TrackId)
            if len(subhalo_index) == 0:
                raise LookupError("The requested TrackId does not exist in the snapshot of interest.")
            subhalo_index = subhalo_index[0]

            return self.LoadSubhalos(snap_nr, subhalo_index=subhalo_index, property_selection=fields)

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
            SnapshotOfBirth = self.GetTrackSnapshot(TrackId, self.SnapshotIdList.max())['SnapshotOfBirth']
        except:
            SnapshotOfBirth = self.GetTrackSnapshot(TrackId, self.SnapshotIdList.max())['SnapshotIndexOfBirth']
        snapshots_to_load = self.SnapshotIdList[self.SnapshotIdList >= SnapshotOfBirth]

        track = np.array(
            [self.GetTrackSnapshot(TrackId, snap_nr, fields=fields)
             for snap_nr in snapshots_to_load])

        scales = np.array([self.GetScaleFactor(snap_nr) for snap_nr in snapshots_to_load])

        return append_fields(
            track, ['Snapshot', 'ScaleFactor'], [snapshots_to_load, scales], usemask=False)

    #===========================================================================
    # Functions to identify progenitor subhaloes.
    #===========================================================================
    def GetSinkProgenitors(self, TrackId):
        """
        Returns the TrackId of subhaloes that sunk into the subhalo of interest.

        Parameters
        ==========
        TrackId: int
            The TrackId of the subhalo whose sunk progenitors we are interested in.

        Returns
        =======
        np.ndarray of int
            The TrackId of subhaloes that sunk into the subhalo of interest.
        """

        # If we are calling this function, we might need SnapshotOfDeath
        # and SnapshotOfSink information more than once. Since these datasets are
        # unique and based on the last available output. Let's ensure they remain
        # loaded if not already.
        self.__load_merger_snapshot_information()

        if not hasattr(self, "SinkTrackId"):
            self.__SinkTrackId = self.LoadSubhalos(self.SnapshotIdList.max(), property_selection=['SinkTrackId'])['SinkTrackId']
        mask_progenitors = self.__SinkTrackId == TrackId

        if self.__sorted_catalogues:
            return np.flatnonzero(mask_progenitors & self.__mask_sunk_subhaloes)
        else:
            return np.sort(self.LoadSubhalos(self.SnapshotIdList.max(), property_selection=["TrackId"])["TrackId"][mask_progenitors & self.__mask_sunk_subhaloes])

    def GetDisruptionProgenitors(self, TrackId):
        """
        Returns the TrackId of subhaloes that disrupted into the subhalo of interest.

        Parameters
        ==========
        TrackId: int
            The TrackId of the subhalo whose disrupted progenitors we are interested in.

        Returns
        =======
        np.ndarray of int
            The TrackId of subhaloes that disrupted into the subhalo of interest.
        """

        # If we are calling this function, we might need SnapshotOfDeath
        # and SnapshotOfSink information more than once. Since these datasets are
        # unique and based on the last available output. Let's ensure they remain
        # loaded if not already.
        self.__load_merger_snapshot_information()

        if not hasattr(self, "__DescendantTrackId"):
            self.__DescendantTrackId = self.LoadSubhalos(self.SnapshotIdList.max(), property_selection=['DescendantTrackId'])['DescendantTrackId']
        mask_progenitors = self.__DescendantTrackId == TrackId

        if self.__sorted_catalogues:
            return np.flatnonzero(mask_progenitors & self.__mask_disrupted_subhaloes)
        else:
            return np.sort(self.LoadSubhalos(self.SnapshotIdList.max(), property_selection=["TrackId"])["TrackId"][mask_progenitors & self.__mask_disrupted_subhaloes])

    def GetAllProgenitors(self, TrackId, only_direct_progenitors=False):
        """
        Returns the TrackId of all subhaloes that merged into the subhalo of interest.
        If direct_progenitors is set to True, only returns those that merged
        directly onto the subhalo of interest.

        Parameters
        ==========
        TrackId: int
            The TrackId of the subhalo whose progenitors we are interested in.
        only_direct_progenitors: bool, opt
            Whether to only return the direct progenitors of the subhalo of interest.
            Defaults to False.
        """

        if only_direct_progenitors:
            return np.hstack([self.GetDisruptionProgenitors(TrackId), self.GetSinkProgenitors(TrackId)])
        else:
            all_progenitors = [] # To keep track of all progenitors we found

            # We start a queue with the request TrackId, and keep on adding and
            # removing subhaloes until we do not need to find the progenitors of
            # any other.
            queue = [TrackId]
            while len(queue) > 0:

                # Grab the next TrackId in the queue.
                current_TrackId = queue.pop(0)

                # Get the progenitors of the current TrackId.
                current_TrackId_progenitors = self.GetAllProgenitors(current_TrackId, only_direct_progenitors=True).tolist()
                all_progenitors.extend(current_TrackId_progenitors)
                queue.extend(current_TrackId_progenitors)

            return np.sort(np.asarray(all_progenitors))

    #===========================================================================
    # Scale factor handling.
    #===========================================================================
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

    #===========================================================================
    # Unit handling.
    #===========================================================================
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

    #===========================================================================
    # Internal helper functions.
    #===========================================================================
    def __LoadSubhalos_UnsortedCatalogues(self, snap_nr=None, subhalo_index=None, property_selection=None, show_progress=False):
        """
        Load subhalos from the unsorted HBT-HERONS catalogues, from a single
        simulation output.

        Parameters
        ==========
        snap_nr: int, opt
            The snapshot number we are interested in. It defaults to the
            latest snapshot with available catalogues.
        subhalo_index: int or np.ndarray, opt
            Only load the subhalo(es) in the specified entries, which must be
            in ascending order if multiple are requested. If not provided, all
            subhaloes will be loaded. Note that this value does NOT correspond
            to the TrackId of the subhalo.
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
        total_number_subhaloes = self.GetNumberOfSubhalos(snap_nr)

        if subhalo_index is not None:
            if not isinstance(subhalo_index, (np.integer, int, np.ndarray)):
                raise TypeError("Parameter TrackId is not of the required type (int or np.ndarray).")

            if isinstance(subhalo_index, (np.integer, int)):
                subhalo_index = np.array([subhalo_index])

            if subhalo_index.max() >= total_number_subhaloes:
                raise ValueError(f"Largest requested subhalo index ({subhalo_index.max()}) is larger than the number of existing subhaloes ({total_number_subhaloes})")

        # Handle defaults, and list inputs
        if property_selection is None:
            with h5py.File(self.GetFileName(snap_nr), 'r') as subfile:
                property_selection = subfile['Subhalos'].dtype.names
        else:
            if type(property_selection) is list:
                property_selection = tuple(property_selection)

        with h5py.File(self.GetFileName(snap_nr), 'r') as subfile:
            # Number of subfiles we may need to iterate over.
            number_subfiles = subfile['NumberOfFiles'][0]

            # Create subhalo array with custom dtypes for the requested properties.
            subhaloes_dtype = generate_custom_array_dtypes(subfile['Subhalos'], property_selection)
            subhaloes_data  = np.empty(len(subhalo_index) if subhalo_index is not None else number_subhaloes, dtype=subhaloes_dtype)

        # Create counters to keep track where to position loaded subhalo data. array for the subhaloes, and keep on filling it as we load each
        subfile_offset = 0
        for subfile_nr in range(number_subfiles):

            if show_progress:
                print(".", end="")

            # This ensures we only iterate over the subfiles that contain the
            # requested subhalo entries. If we are loading everything, this will
            # never trigger.
            if subfile_offset > subhalo_index.max() if subhalo_index is not None else False:
                break

            with h5py.File(self.GetFileName(snap_nr, subfile_nr), 'r') as subfile:
                number_subhaloes = subfile['Subhalos'].shape[0]

                # Nothing to load
                if number_subhaloes == 0:
                    continue

                # Iterate over subfiles until we find our target(s).
                if subhalo_index is not None:

                    # We find which of the requested entries are located in this
                    # subfile.
                    present_subhalo_index_mask = (subhalo_index >= subfile_offset) & (subhalo_index < subfile_offset + number_subhaloes)

                    # Load into the entries we want
                    for property in property_selection:
                        subhaloes_data[property][present_subhalo_index_mask] = subfile['Subhalos'][property][subhalo_index - subfile_offset]

                    subfile_offset += number_subhaloes

                else: # Load everything
                    for property in property_selection:
                        subhaloes_data[property] = subfile['Subhalos'][property]

        if show_progress:
            print()

        return subhaloes_data

    def __LoadSubhalos_SortedCatalogues(self, snap_nr=None, TrackId=None, property_selection=None):
        """
        Load subhalos from the sorted HBT-HERONS catalogues, from a single
        simulation output.

        Parameters
        ==========
        snap_nr: int, opt
            The snapshot number we are interested in. It defaults to the
            latest snapshot with available catalogues.
        TrackId: int or np.ndarray, opt
            Only load the subhalo(es) with the specified TrackId, which must be
            in ascending order if multiple are requested. If not provided, all
            subhaloes will be loaded.
        property_selection: tuple or list of strings, opt
            If specified, only load the specified properties. If not provided,
            all subhalo properties will be loaded.

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
        number_subhaloes = self.GetNumberOfSubhalos(snap_nr)

        if TrackId is not None:
            if not isinstance(TrackId, (np.integer, int, np.ndarray)):
                raise TypeError("Parameter TrackId is not of the required type (int or np.ndarray).")

            if isinstance(TrackId, (np.integer, int)):
                TrackId = np.array([TrackId])

            if TrackId.max() >= number_subhaloes:
                raise ValueError(f"Largest requested TrackId ({TrackId.max()}) is larger than the number of existing subhaloes ({number_subhaloes})")

        subhaloes_data  = []
        with h5py.File(self.GetFileName(snap_nr), 'r') as catalogue_file:

            # If we have not specified specific properties to load, we load
            # everything.
            if property_selection is None:
                property_selection = list(catalogue_file['Subhalos'].keys())

            subhaloes_dtype = generate_custom_array_dtypes(catalogue_file['Subhalos'], property_selection)
            subhaloes_data  = np.empty(len(TrackId) if TrackId is not None else number_subhaloes, dtype=subhaloes_dtype)

            for property in property_selection:
                if TrackId is None:
                    subhaloes_data[property] = catalogue_file[f"Subhalos/{property}"][()]
                else:
                    subhaloes_data[property] = catalogue_file[f"Subhalos/{property}"][TrackId]

        return subhaloes_data

    def __load_merger_snapshot_information(self):
        """
        Identify how each subhalo in the simulation merged with others and store
        for later use.
        """
        if not hasattr(self, "__mask_disrupted_subhaloes") or not hasattr(self, "__mask_sunk_subhaloes"):

            subhalo_data = self.LoadSubhalos(self.SnapshotIdList.max(), property_selection=['SnapshotOfDeath', 'SnapshotOfSink'])

            # Identify subhaloes that disrupted or underwent unresolved sinking.
            self.__mask_disrupted_subhaloes = ( subhalo_data['SnapshotOfDeath'] != -1) & \
                                              ((subhalo_data['SnapshotOfSink']  == -1) | (subhalo_data['SnapshotOfSink'] > subhalo_data['SnapshotOfDeath']))
            # Identify subhaloes that sunk
            self.__mask_sunk_subhaloes = (subhalo_data['SnapshotOfDeath'] != -1) & \
                                         (subhalo_data['SnapshotOfSink']  == subhalo_data['SnapshotOfDeath'])