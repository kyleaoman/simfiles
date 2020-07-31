import warnings
from importlib.util import spec_from_file_location, module_from_spec
from os.path import expanduser, dirname, join
from ._hdf5_io import hdf5_get

# SimFiles is a dict with added features, notably __getattr__ and __setattr__,
# and automatic loading of data from simulation files as defined using a config
# file.


def dealias(func):
    def dealias_wrapper(self, key, *args, **kwargs):
        if not key.startswith('_') and hasattr(self, '_aliases'):
            key = self._aliases.get(key, key)
        return func(self, key, *args, **kwargs)
    return dealias_wrapper


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

        self._aliases = dict()
        self._dealiases = dict()
        try:
            aliaspath = dirname(config.__file__)
            aliasfile = join(aliaspath, config.aliasfile)
        except AttributeError:
            pass
        else:
            with open(aliasfile) as f:
                lines = f.readlines()
            for line in lines:
                v, k = line.strip().split()
                if k.startswith('_'):
                    raise ValueError("Aliases may not start with '_'.")
                self._aliases[k] = v
                self._dealiases[v] = k
            if not set(self._aliases.values()).issubset(
                    self._extractors.keys()
            ):
                unknown = set(self._aliases.values()) - \
                    set(self._extractors.keys())
                warnings.warn(
                    'Aliases exist for unknown keys:\n {:s}.'.format(
                        '\n '.join(unknown)
                    ),
                    RuntimeWarning
                )

        return

    @dealias
    def __setattr__(self, key, value):
        return self.__setitem__(key, value)

    @dealias
    def __getattr__(self, key):
        try:
            return self.__getitem__(key)
        except KeyError:
            raise AttributeError("'SimFiles' object has no attribute '{:s}'."
                                 .format(key))

    __getitem__ = dealias(dict.__getitem__)
    __setitem__ = dealias(dict.__setitem__)

    @dealias
    def __delitem__(self, key):
        if not self.share_mode:
            return super().__delitem__(key)
        else:
            return

    @dealias
    def __delattr__(self, key):
        if key in self.keys():
            if not self.share_mode:
                return super().__delitem__(key)
        else:
            return super().__delattr__(key)

    def load(self, keys=tuple(), filetype=None, intervals=None, verbose=True):
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
            True).
        """

        keys = [self._aliases.get(key, key) for key in keys]

        loaded_keys = set()

        try:
            keys = tuple(keys)
        except TypeError:
            raise TypeError("SimFiles.load: keys must be iterable.")

        if intervals is None:
            intervals = (None, ) * len(keys)
        for key, interval in zip(keys, intervals):
            loaded_keys.update(self._load_key(
                key,
                filetype=filetype,
                interval=interval,
                verbose=verbose
            ))

        return loaded_keys

    def fields(self, keytype='all', aliases=True):
        """
        Return a list of available keys, optionally for a specific keytype.

        Parameters
        ----------
        keytype : str
            Specify which type of keys to include. This can be one of the
            keytypes defined in the extractors, or 'all', or 'aliased'
            (default: 'all').

        aliases : bool
            If True, the keys will be replaced by their aliases in the list
            (default: True).
        """
        if keytype == 'all':
            retval = [k for k in self._extractors.keys()]
        elif keytype == 'aliased':
            retval = list(self._aliases.values())
        else:
            retval = [k for k, E in self._extractors.items()
                      if E.keytype == keytype]
        if aliases:
            retval = [self._dealiases.get(k, k) for k in retval]
        return retval

    def _dependencies(self, _dependencies_list, filetype=None, interval=None,
                      verbose=True):

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

    def _load_key(self, key, filetype=None, interval=None, verbose=True):

        loaded_keys = set()

        if key in self:
            if self.share_mode:
                return tuple()
            else:
                if verbose:
                    warnings.warn(
                        "SimFiles._load_key: overwriting key '{:s}' (alias: "
                        "{:s}), may be possible to suppress by changing load "
                        "order.".format(key, self._dealiases.get(key, 'None')),
                        RuntimeWarning
                    )

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
        if verbose:
            alias = self._dealiases.get(key, None)
            astr = ' (alias: {:s})'.format(alias) if alias else ''
            print('SimFiles: loaded {:s}{:s}.'.format(key, astr))

        loaded_keys.update((key, ))

        return loaded_keys
