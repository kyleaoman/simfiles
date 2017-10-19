import numpy as np
import warnings
from collections import namedtuple
from utilities.hdf5_io import hdf5_get

class SimFiles(dict):
    def __init__(self, snap_id, configfile=None):

        self.snap_id = snap_id
        self.configfile = configfile

        self._read_config()

        return
        
    def _read_config(self):

        config = dict()
        try:
            execfile(self.configfile, config)
        except IOError:
            raise IOError("SimFiles: configfile '" + self.configfile + "' not found.")
        try:
            snapshots = config['snapshots']
        except KeyError:
            raise ValueError("SimFiles: configfile missing 'snapshots' definition.")            

        try:
            self._snapshot = snapshots[self.snap_id]
        except KeyError:
            raise ValueError("SimFiles: unknown snapshot (not defined in configfile).")

        del snapshots
        
        try:
            self._extractors = config['extractors']
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
    
    def load(self, keys=None, filetype=None):

        loaded_keys = set()
        
        if type(keys) is not tuple:
            raise ValueError('SimFiles.load: keys must be tuple.')
        
        for key in keys:
            loaded_keys.update(self._load_key(key, filetype=filetype))

        return loaded_keys

    def fields(self, keytype = 'all'):
        if keytype == 'all':
            return [k for k in self._extractors.keys()]
        else:
            return [k for k, E in self._extractors.items() if E.keytype == keytype]
    
    def _dependencies(self, _dependencies_list, filetype=None):

        loaded_keys = set()
        
        for k in _dependencies_list:
            if k not in self:
                loaded_keys.update(self._load_key(k, filetype=filetype))

        return loaded_keys

    def _load_key(self, key, filetype=None):
    
        loaded_keys = set()
    
        if key in self:
            warnings.warn("SimFiles._load_key: overwriting key '"+key+"', may be possible to suppress by changing load order.", RuntimeWarning)

        loaded_keys.update(self._dependencies(self._extractors[key].dependencies, filetype=filetype))

        E = self._extractors[key]
        path, fname = None, None
        try:
            path, fname = self._snapshot[E.filetype if filetype is None else filetype]
        except KeyError:
            raise ValueError("SimFiles: filetype '" + E.filetype if filetype is None else filetype + "' unknown.")
        self[key] = E.convert(self, hdf5_get(path, fname, E.hpath, attr=E.attr), path, fname, E.hpath) * E.units
            
        loaded_keys.update((key, ))
        
        return loaded_keys

        
