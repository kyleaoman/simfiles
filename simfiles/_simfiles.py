import warnings
from kyleaoman_utilities.hdf5_io import hdf5_get
from importlib.util import spec_from_file_location, module_from_spec
from os.path import expanduser

# SimFiles is a dict with added features, notably __getattr__ and __setattr__,
# and automatic loading of data from simulation files as defined using a config
# file.


class SimFiles(dict):
    """
    Provides a generic interface to simulation hdf5 files.

    SimFiles is a dict with added features, notably __getattr__ and
    __setattr__, and automatic loading of data from simulation files based on a
    configuration file.

    Parameters
    ----------
    snap_id : index
        An identifier for a specific simulation snapshot. The exact format is
        defined in the configuration for the simulation in question.

    configfile : str
        Path to the configuration file to use (default: None).

    ncpu : int
        Number of processors on which to run (default: 2).

    share_mode : bool
        Setting 'True' disables the __delitem__ method (default: False) and
        suppresses warnings for repeated loading of the same keys.

    single_file : int
        Specify to load from only a specific hdf5 file 'piece' of the snapshot.
        Assumes 'pieces' end in '.X.hdf5' where X is an integer.

    Returns
    -------
    out : SimFiles
        A SimFiles object configured using the file provided.

    Examples
    --------
    The following example sets up a SimFiles instance, loads a few keys, and
    accesses the loaded data, for APOSTLE simulation data on the cavi system::

        from simfiles import SimFiles
        from simfiles.configs.APOSTLE_cavi import __file__ as configfile
        import namedtuple
        snap_id = namedtuple('snap_id', ['res', 'phys', 'vol', 'snap'])
        mysnap = snap_id(res=3, phys='hydro', vol=1, snap=127)
        SF = SimFiles(mysnap, configfile=configfile)
        SF.load(keys=('m_s', 'xyz_s', 'vxyz_s'))
        # print mass and coordinates of one particle
        # both dict-like and attribute-like access are supported
        # this config file supports units via astropy.units
        print(SF.m_s[0], SF['xyz_s'][0])
    """

    def __init__(self, snap_id, configfile=None, ncpu=2, share_mode=False,
                 single_file=None):

        self.snap_id = snap_id
        self.configfile = expanduser(configfile)
        self.ncpu = ncpu
        self.share_mode = share_mode
        self.single_file = single_file

        self._read_config()

        return

    def _read_config(self):

        try:
            spec = spec_from_file_location('config', self.configfile)
            config = module_from_spec(spec)
            spec.loader.exec_module(config)
        except FileNotFoundError:
            raise FileNotFoundError("SimFiles: configfile '{:s}' not found."
                                    .format(self.configfile))
        try:
            snapshots = config.snapshots
        except AttributeError:
            raise ValueError("SimFiles: configfile missing 'snapshots' "
                             "definition.")

        try:
            self._snapshot = snapshots[self.snap_id]
        except KeyError:
            raise ValueError("SimFiles: unknown snapshot (not defined in "
                             "configfile).")

        try:
            self._extractors = config.extractors
        except AttributeError:
            raise ValueError("Simfiles: configfile missing 'extractors' "
                             "definition.")

        return

    def __setattr__(self, key, value):
        return self.__setitem__(key, value)

    def __getattr__(self, key):
        try:
            return self.__getitem__(key)
        except KeyError:
            raise AttributeError("'SimFiles' object has no attribute '{:s}'."
                                 .format(key))

    def __delitem__(self, key):
        if not self.share_mode:
            return super().__delitem__(key)
        else:
            return

    def load(self, keys=tuple(), filetype=None, intervals=None, verbose=False):
        """
        Load data for a set of keys.

        Parameters
        ----------
        keys : iterable
            List of keys to load (default: tuple()).

        filetype : str
            Advanced use only, override filetype defined in config file
            (default: None).

        intervals : iterable
            List containing lists of 2-tuples, one for each key. Each 2-tuple
            represents an interval of indices from the underlying data table
            to load (default: None).

        verbose : bool
            Setting 'True' prints messages upon loading each key (default:
            False).
        """
        loaded_keys = set()

        try:
            keys = tuple(keys)
        except TypeError:
            raise TypeError("SimFiles.load: keys must interpretable as tuple.")

        if (keys != tuple()) and (intervals is None):
            intervals = (None, ) * len(keys)
        for key, interval in zip(keys, intervals):
            loaded_keys.update(self._load_key(
                key,
                filetype=filetype,
                interval=interval,
                verbose=verbose
            ))

        return loaded_keys

    def fields(self, keytype='all'):
        """
        Return a list of available keys, optionally for a specific keytype.

        Parameters
        ----------
        keytype : str
            Specify which type of keys to include (default: 'all').
        """
        if keytype == 'all':
            return [k for k in self._extractors.keys()]
        else:
            return [k for k, E in self._extractors.items()
                    if E.keytype == keytype]

    def _dependencies(self, _dependencies_list, filetype=None, interval=None,
                      verbose=False):

        loaded_keys = set()

        for k in _dependencies_list:
            if k not in self:
                loaded_keys.update(self._load_key(
                    k,
                    filetype=filetype,
                    interval=interval,
                    verbose=verbose
                ))

        return loaded_keys

    def _load_key(self, key, filetype=None, interval=None, verbose=False):

        loaded_keys = set()

        if key in self:
            if self.share_mode:
                return tuple()
            else:
                warnings.warn("SimFiles._load_key: overwriting key '{:s}', may"
                              " be possible to suppress by changing load "
                              "order.".format(key), RuntimeWarning)

        loaded_keys.update(self._dependencies(
            self._extractors[key].dependencies,
            filetype=filetype,
            interval=interval,
            verbose=verbose
        ))

        E = self._extractors[key]
        path, fname = None, None
        use_filetype = E.filetype if filetype is None else filetype
        try:
            path, fname = self._snapshot[use_filetype]
        except KeyError:
            raise ValueError("SimFiles: filetype '{:s}' unknown."
                             .format(use_filetype))
        if (self.single_file is not None) and \
           ((use_filetype == 'particle') or (use_filetype == 'snapshot')):
            # will force loading only one file for particles
            fname = fname + '.{0:.0f}'.format(self.single_file)
        self[key] = E.convert(
            self,
            hdf5_get(
                path,
                fname,
                E.hpath,
                attr=E.attr,
                ncpu=self.ncpu,
                interval=interval
            ),
            path,
            fname,
            E.hpath
        )
        if E.units is not None:
            self[key] = self[key] * E.units
        if E.unit_convert is not None:
            self[key] = self[key].to(E.unit_convert)

        loaded_keys.update((key, ))

        return loaded_keys
