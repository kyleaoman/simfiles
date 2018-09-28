import numpy as np
import warnings
from collections import namedtuple
from kyleaoman_utilities.hdf5_io import hdf5_get
from importlib.util import spec_from_file_location, module_from_spec
from os.path import expanduser

# SimFiles is a dict with added features, notably __getattr__ and __setattr__, and automatic loading
# of data from simulation files as defined using a config file.

class SimFiles(dict):
    def __init__(self, snap_id, configfile=None, ncpu=2, share_mode=False):

        self.snap_id = snap_id
        self.configfile = configfile
        self.ncpu = ncpu
        self.share_mode = share_mode

        self._read_config()

        return
        
    def _read_config(self):

        try:
            spec = spec_from_file_location('config', expanduser(self.configfile))
            config = module_from_spec(spec)
            spec.loader.exec_module(config)
        except FileNotFoundError:
            raise FileNotFoundError("SimFiles: configfile '" + self.configfile + "' not found.")
        try:
            snapshots = config.snapshots
        except AttributeError:
            raise ValueError("SimFiles: configfile missing 'snapshots' definition.")            

        try:
            self._snapshot = snapshots[self.snap_id]
        except KeyError:
            raise ValueError("SimFiles: unknown snapshot (not defined in configfile).")

        del snapshots
        
        try:
            self._extractors = config.extractors
        except AttributeError:
            raise ValueError("Simfiles: configfile missing 'extractors' definition.")
        
        return        

    def __setattr__(self, key, value):
        return self.__setitem__(key, value)

    def __getattr__(self, key):
        try:
            return self.__getitem__(key)
        except KeyError:
            raise AttributeError("'SimFiles' object has no attribute '"+str(key)+"'.")

    def __delitem__(self, key):
        if self.share_mode == False:
            return super(SimFiles, self).__delitem__(key)
        else:
            return
    
    def load(self, keys=None, filetype=None, intervals=None):

        loaded_keys = set()
        
        if type(keys) is not tuple:
            raise ValueError('SimFiles.load: keys must be tuple.')

        if (keys != None) and (intervals == None):
            intervals = (None, ) * len(keys)
        for key, interval in zip(keys, intervals):
            loaded_keys.update(self._load_key(key, filetype=filetype, interval=interval))

        return loaded_keys

    def fields(self, keytype = 'all'):
        if keytype == 'all':
            return [k for k in self._extractors.keys()]
        else:
            return [k for k, E in self._extractors.items() if E.keytype == keytype]
    
    def _dependencies(self, _dependencies_list, filetype=None, interval=None):

        loaded_keys = set()
        
        for k in _dependencies_list:
            if k not in self:
                loaded_keys.update(self._load_key(k, filetype=filetype, interval=interval))

        return loaded_keys

    def _load_key(self, key, filetype=None, interval=None):
    
        loaded_keys = set()
    
        if key in self:
            if self.share_mode:
                return tuple()
            else:
                warnings.warn("SimFiles._load_key: overwriting key '"+key+"', may be possible to suppress by changing load order.", RuntimeWarning)

        loaded_keys.update(self._dependencies(self._extractors[key].dependencies, filetype=filetype, interval=interval))

        E = self._extractors[key]
        path, fname = None, None
        try:
            path, fname = self._snapshot[E.filetype if filetype is None else filetype]
        except KeyError:
            raise ValueError("SimFiles: filetype '" + E.filetype if filetype is None else filetype + "' unknown.")
        self[key] = E.convert(self, hdf5_get(path, fname, E.hpath, attr=E.attr, ncpu=self.ncpu, interval=interval), path, fname, E.hpath) 
        if E.units is not None:
            self[key] = self[key] * E.units
        if E.unit_convert is not None:
            self[key] = self[key].to(E.unit_convert)
            
        loaded_keys.update((key, ))
        
        return loaded_keys
