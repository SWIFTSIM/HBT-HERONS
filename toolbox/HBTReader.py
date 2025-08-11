"""
Class to read the HBT outputs.
NOTE: only use for small simulations, as this method is very slow to
retrieve the properties of all subhaloes. Suggested alternative is to
instead run the script "./SortCatalogues.py", which will sort catalogues in
ascending TrackId order.

To use it, initialize the reader with the parameter file under the subhalo directory, e.g.,

  from HBTReader import HBTReader

  reader=HBTReader('subcat')

  snapshotnumber=-1 # or 0~MaxSnap. -1 means last snapshot

  subs=reader.LoadSubhalos(snapshotnumber) #load all

  nbound=reader.LoadSubhalos(snapshotnumber, 'Nbound') #only Nbound

  sub2=reader.LoadSubhalos(snapshotnumber, subindex=2) #only subhalo 2

  track2=reader.GetTrackEvolution(2) #track 2

"""
import h5py
import os.path
import numpy as np
from glob import glob
# TODO: check whether they are required
from numpy.lib.recfunctions import append_fields

class ConfigReader:
    """
    Class to read the HBT-HERONS configuration files.
    """
    Options = {}

    def __init__(self, path_to_config_file):
        """
        Loads which options are present in the HBT-HERONS parameter file, and
        their values.

        Parameters
        ==========
        path_to_config_file: str
            Path to the HBT-HERONS configuration file.
        """

        with open(path_to_config_file, 'r') as file:
            for line in file:
                self.__parse_line(line)

    def __parse_line(self, line):
        """
        Extracts the name and value of a HBT-HERONS parameter, if it is present
        in the provided line.

        Parameters
        ==========
        line: str
            A single line of a parameter file.
        """
        pair = line.lstrip().split("#", 1)[0].split("[", 1)[0].split()
        if len(pair) == 2:
            self.Options[pair[0]] = pair[1]
        elif len(pair) > 2:
            self.Options[pair[0]] = pair[1:]

    def __getitem__(self, ParameterName):
        return self.Options[ParameterName]

def get_hbt_snapnum(snapname):
    return int(snapname.rsplit('SubSnap_')[1].split('.')[0])

class HBTReader:
    """
    Class to read HBT-HERONS catalogues.
    """

    def __init__(self, base_dir):
        """
        Initialize HBTReader to read data from base_dir where all the outputs
        are saved. A parameter file must exist there (Parameters.log dumped by
        HBT-HERONS during runtime).
        """

        self.Options = ConfigReader(base_dir +'/Parameters.log').Options
        self.base_dir = base_dir

        # To know which files to open.
        self.MinimumSnapshotIndex = int(self.Options['MinSnapshotIndex'])
        self.MaximumSnapshotIndex = int(self.Options['MaxSnapshotIndex'])
        if "SnapshotIdList" in self.Options:
            self.SnapshotIdList = np.array(self.Options["SnapshotIdList"], int)
        else:
            self.SnapshotIdList = np.arange(self.MinimumSnapshotIndex, self.MaximumSnapshotIndex)

        # Generate an f-formated list of files
        self._file_list = sorted(glob(self.base_dir+'/*/SubSnap_*.0.hdf5'),key=get_hbt_snapnum)
        self._file_list = [path.replace(".0.hdf5",".{subfile_nr}.hdf5") for path in self._file_list]



        # Do we have the same number of files as we expect? If not, remove the
        # missing catalogues.
        if len(self._file_list) != len(self.SnapshotIdList):
            print(f"HBT-HERONS run not finished yet. Only found {len(self._file_list)} outputs out of {len(self.SnapshotIdList)} total.")
            self.SnapshotIdList = self.SnapshotIdList[:len(self._file_list)]

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

        return self._file_list[index[0]].format(subfile_nr = subfile_nr)

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

    def LoadParticleIDs(self, snap_nr=None, subhalo_index=None, filetype='Sub'):
        """
        Load the particles that are bound or are part of the source subhalo
        for the specified subhalo.

        Parameters
        ==========
        snap_nr : int, opt
            Snapshot number of the catalogue we are interested in. It defaults
            to the last snapshot with currently available catalogues.
        subhalo_index: int
            The array entry to load from the subhalo catalogues. Note
            that this does NOT correspond to the TrackId of the subhalo, as
            the catalogues are not sorted in TrackId.

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

        load_single_subhalo = subhalo_index is not None
        if (load_single_subhalo):
            if not isinstance(subhalo_index, (np.integer, int)):
                raise TypeError("Parameter subhalo_index is not of the required type (int).")
            if subhalo_index >= self.GetNumberOfSubhalos(snap_nr):
                raise ValueError(f"Selected subhalo entry ({subhalo_index}) is larger than the number of existing subhaloes ({self.GetNumberOfSubhalos(snap_nr)})")

        # Determine number of files for requested output
        with h5py.File(self.GetFileName(snap_nr), 'r') as subfile:
            number_subfiles = subfile['NumberOfFiles'][0]

        # Determine number of files for requested output
        with h5py.File(self.GetFileName(snap_nr), 'r') as subfile:
            number_subfiles = subfile['NumberOfFiles'][0]

        offset = 0
        subhalo_particles = []
        for subfile_nr in range(number_subfiles):
            with h5py.File(self.GetFileName(snap_nr, subfile_nr), 'r') as subfile:
                nsub = subfile['Subhalos'].shape[0]

                # Nothing to load
                if nsub == 0:
                    continue

                # Keep iterating over subfiles until we find our target entry.
                if load_single_subhalo:
                    if (offset+nsub > subhalo_index):
                        subhalo_particles.append(subfile['SubhaloParticles'][subhalo_index-offset])
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