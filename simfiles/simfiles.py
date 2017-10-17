import numpy as np
import warnings
from collections import namedtuple
from hdf5_io import get as hdf5_get

class SimFiles(dict):
    def __init__(self, snap_id, configfile=None):

        config = dict()
        try:
            execfile(configfile, config)
        except IOError:
            raise IOError("SimFiles: configfile '" + configfile + "' not found.")
        try:
            self._snapshots = config['snapshots']
        except KeyError:
            raise ValueError("SimFiles: configfile missing 'snapshots' definition.")            

        try:
            self._snapshot = self._snapshots[snap_id]
        except KeyError:
            raise ValueError("SimFiles: unknown snapshot (not defined in configfile).")
        
        try:
            self._Extractors = config['extractors']
        except KeyError:
            raise ValueError("Simfiles: configfile missing 'extractors' definition.")
        
        return

    def __setattr__(self, key, value):
        return self.__setitem__(key, value)

    def __getattr__(self, key):
        try:
            return self.__getitem__(key)
        except KeyError:
            raise AttributeError("'SimFiles' object has no attribute '"+str(key)+"'.")
    
    def load(self, filetype, keys=tuple()):

        loaded_keys = set()

        if (type(keys) != tuple):
            raise ValueError('SimFiles.load(file_id, keys=tuple()): keys must be tuple.')
        
        for key in keys:
            try:
                path, filename = self._snapshot[filetype]
            except KeyError:
                raise ValueError("SimFiles: filetype '" + filetype + "' unknown.")
            loaded_keys.update(self._load_key(path, filename, key))

        return loaded_keys

    def fields(self, keytype = 'all'):
        if keytype == 'all':
            return [k for k in self._Extractors.keys()]
        else:
            return [k for k, E in self._Extractors.items() if E.keytype == keytype]
    
    def _dependencies(self, _dependencies_list, path, fname):

        loaded_keys = set()
        
        for k in _dependencies_list:
            if k not in self:
                loaded_keys.update(self._load_key(path, fname, k))

        return loaded_keys

    def _load_key(self, path, fname, key):
    
        loaded_keys = set()
    
        if key in self:
            warnings.warn("ApostleFileset._load_key: overwriting key '"+key+"', may be possible to suppress by changing load order", RuntimeWarning)

        loaded_keys.update(self._dependencies(self._Extractors[key].dependencies, path, fname))

        E = self._Extractors[key]
        self[key] = E.convert(self, hdf5_get(path, fname, E.hpath, attr=E.attr), path, fname, E.hpath) * E.units
            
        loaded_keys.update((key, ))
        
        return loaded_keys

        
