"""
Class to read the HBT outputs.

To use it, initialize the reader with the parameter file under the subhalo directory, e.g.,

  from HBTReader import HBTReader

  reader=HBTReader('subcat')

  snapshotnumber=-1 # or 0~MaxSnap. -1 means last snapshot

  subs=reader.LoadSubhalos(snapshotnumber) #load all

  nbound=reader.LoadSubhalos(snapshotnumber, 'Nbound') #only Nbound

  sub2=reader.LoadSubhalos(snapshotnumber, subindex=2) #only subhalo 2

  track2=reader.GetTrack(2) #track 2

"""
# TODO: check whether they are required (possibly Python2 leftover)
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import sys
import h5py
import os.path
import numpy as np
from glob import glob
# TODO: check whether they are required
import numbers
from numpy.lib.recfunctions import append_fields

def PeriodicDistance(x, y, BoxSize, axis=-1):
    d = x-y
    d[d > BoxSize/2] = d[d > BoxSize/2]-BoxSize
    d[d < -BoxSize/2] = d[d < -BoxSize/2]+BoxSize
    return np.sqrt(np.sum(d**2, axis=axis))

def distance(x, y, axis=-1):
    return np.sqrt(np.sum((x-y)**2, axis=axis))

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
    NOTE: only use for small simulations, as this method is very slow to 
    retrieve the properties of all subhaloes. Suggested alternative is to
    instead run the script "./SortCatalogues.py"
    """

    def __init__(self, base_path):
        """ 
        Initialize HBTReader to read data from base_path where all the outputs
        are saved. A parameter file must exist there (Parameters.log dumped by 
        HBT-HERONS during runtime).
        """

        self.Options = ConfigReader(base_path +'/Parameters.log').Options
        self.rootdir = self.Options['SubhaloPath']
        self.MaxSnap = int(self.Options['MaxSnapshotIndex'])
        self.BoxSize = float(self.Options['BoxSize'])
        self.Softening = float(self.Options['SofteningHalo'])

        # To know which files to open.
        self.MinimumSnapshotIndex = int(self.Options['MinSnapshotIndex'])
        self.MaximumSnapshotIndex = int(self.Options['MaxSnapshotIndex'])
        if "SnapshotIdList" in self.Options:
            self.SnapshotIdList = np.array(self.Options["SnapshotIdList"], int)
        else:
            self.SnapshotIdList = np.arange(self.MinimumSnapshotIndex, self.MaximumSnapshotIndex)

        # Generate an f-formated list of files
        file_list = sorted(glob(self.rootdir+'/*/SubSnap_*.0.hdf5'),key=get_hbt_snapnum)
        file_list = [path.replace(".0.hdf5",".{subfile_nr}.hdf5") for path in file_list]

        # Do we have the same number of files as we expect? 
        if len(file_list) != len(self.SnapshotIdList):
            print(f"HBT-HERONS run not finished yet. Only found {len(file_list)} outputs out of {len(self.SnapshotIdList)} total.")
            # Remove the expect (and not-yet-done) catalogues.
            self.SnapshotIdList = self.SnapshotIdList[:len(file_list)]

        if 'MinSnapshotIndex' in self.Options:
            self.MinSnap = int(self.Options['MinSnapshotIndex'])
        else:
            self.MinSnap = 0

    def GetFileName(self, snap_nr, subfile=0, filetype='Sub'):
        """
        Returns the path to a catalogue file.

        Parameters
        ==========
        snap_nr: int
            The snapshot number of the catalogue we are interested in.
        subfile: int, opt
            The subfile we are interested in.
        filetype: str, opt
            Whether to load the bound subhalo information ('Sub') or the source
            subhalo information ('Src').
        """
        if snap_nr < 0:
            snap_nr = self.MaxSnap+1+snap_nr
        if self.nfiles:
            return self.rootdir+'/%03d/' % isnap+filetype+'Snap_%03d.%d.hdf5' \
                % (snap_nr, subfile)
        else:
            return self.rootdir+'/'+filetype+'Snap_%03d.hdf5' % (snap_nr)

    def LoadNestedSubhalos(self, snap_nr=-1):
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

    def LoadSubhalos(self, isnap=-1, selection=None, show_progress=False):
        """
        load subhalos from snapshot isnap (default =-1, means final snapshot;
        isnap<0 will count backward from final snapshot)

        `selection` can be a single field, a list of the field names or a single
        subhalo index. e.g., selection=('Rank', 'Nbound') will load only the Rank
        and Nbound fields of subhaloes. selection=3 will only load subhalo with
        subindex 3. Default will load all fields of all subhaloes.

        ...Note: subindex specifies the order of the subhalo in the file at the
        current snapshot, i.e., subhalo=AllSubhalo[subindex].    subindex==trackId
        for single file output, but subindex!=trackId for mpi multiple-file outputs.

        You can also use numpy slice for selection, e.g., selection=np.s_[:10,
        'Rank','HostHaloId'] will select the 'Rank' and 'HostHaloId' of the first 10
        subhaloes. You can also specify multiple subhaloes by passing a list of
        (ordered) subindex, e.g., selection=((1,2,3),). However, currently only a
        single subhalo can be specified for multiple-file hbt data (not restricted
        for single-file data).
        """

        subhalos = []
        offset = 0
        trans_index = False
        if selection is None:
            selection = np.s_[:]
        else:
            trans_index = isinstance(selection, numbers.Integral)

        if type(selection) is list:
            selection = tuple(selection)

        for i in range(max(self.nfiles, 1)):
            if show_progress:
                sys.stdout.write(".")
                sys.stdout.flush()
            with h5py.File(self.GetFileName(isnap, i), 'r') as subfile:
                nsub = subfile['Subhalos'].shape[0]
                if nsub == 0:
                    continue
                if trans_index:
                    if offset+nsub > selection:
                        subhalos.append(subfile['Subhalos'][selection-offset])
                        break
                    offset += nsub
                else:
                    subhalos.append(subfile['Subhalos'][selection])
        if len(subhalos):
            subhalos = np.hstack(subhalos)
        else:
            subhalos = np.array(subhalos)
        if show_progress:
            print()
        return subhalos

    def GetNumberOfSubhalos(self, snap_nr=-1):
        """
        Returns the total number of subhaloes in a given catalogue output.

        Parameters
        ==========
        snap_nr : int, opt
            Snapshot number of the catalogue we are interested in. It defaults
            to the last snapshot with currently available catalogues.
        """
        with h5py.File(self.GetFileName(snap_nr, 0), 'r') as f:
            if self.nfiles:
                return f['TotalNumberOfSubhalosInAllFiles'][...]
            else:
                return f['Subhalos'].shape[0]

    def LoadParticles(self, isnap=-1, subindex=None, filetype='Sub'):
        """ 
        load subhalo particle list at snapshot isnap.

        if subindex is given, only load subhalo of the given index (the order it
        appears in the file, subindex==trackId for single file output, but not for
        mpi multiple-file outputs). otherwise load all the subhaloes.

        default filetype='Sub' will load subhalo particles. set filetype='Src' to
        load source subhalo particles instead (for debugging purpose only).
        """

        subhalos = []
        offset = 0
        for i in range(max(self.nfiles, 1)):
            with h5py.File(self.GetFileName(isnap,  i, filetype), 'r') as subfile:
                if subindex is None:
                    subhalos.append(subfile[filetype+'haloParticles'][...])
                else:
                    nsub = subfile[filetype+'haloParticles'].shape[0]
                    if offset+nsub > subindex:
                        subhalos.append(
                            subfile[filetype+'haloParticles'][subindex-offset])
                        break
                    offset += nsub
        subhalos = np.hstack(subhalos)
        return subhalos

    def GetParticleProperties(self, subindex, isnap=-1):
        """
        load subhalo particle properties for subhalo with index subindex (the
        order it appears in the file, subindex==trackId for single file output, but
        not for mpi multiple-file outputs)
        """

        offset = 0
        for i in range(max(self.nfiles, 1)):
            with h5py.File(self.GetFileName(isnap,  i), 'r') as subfile:
                nsub = subfile['Subhalos'].shape[0]
                if offset+nsub > subindex:
                    # for compatibility with old data
                    try:
                        return subfile['ParticleProperties/Sub%d' % (subindex-offset)][...]
                    except:
                        return subfile['ParticleProperties'][subindex-offset]
                offset += nsub
        raise RuntimeError("subhalo %d not found" % subindex)

    def GetSub(self, trackId, isnap=-1):
        """
        load a subhalo with the given trackId at snapshot isnap
        """
        if self.nfiles:
            subid = np.flatnonzero(self.LoadSubhalos(
                isnap, 'TrackId') == trackId)[0]
        else:
            subid = trackId
        return self.LoadSubhalos(isnap, subid)

    def GetTrackSnapshot(self, trackId, isnap, fields=None):
        """
        Get track information for a single snapshot
        """
        s = self.GetSub(trackId, isnap)
        if fields is not None:
            return s[fields]
        return s

    def GetTrack(self, trackId, fields=None):
        """ 
        load an entire track of the given trackId 
        """
        track = []
        snaps = []
        scales = []
        snapbirth = self.GetSub(trackId)['SnapshotOfBirth']
        if hasattr(snapbirth, '__iter__'):
            snapbirth = snapbirth[0]
        snaps = np.arange(snapbirth, self.MaxSnap+1, dtype=int)
        track = np.array(
            [self.GetTrackSnapshot(trackId, isnap, fields=fields)
             for isnap in snaps])
        scales = np.array([self.GetScaleFactor(isnap) for isnap in snaps])
        return append_fields(
            track, ['Snapshot', 'ScaleFactor'], [snaps, scales], usemask=False)

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
        return dict([(i, self.GetScaleFactor(i))
                     for i in range(self.MinSnap, self.MaxSnap+1)])