# simfiles
Framework to simplify loading data from a set of files corresponding to a simulation snapshot, setup via a configuration file. Examine [`example.py`](https://github.com/kyleaoman/simfiles/blob/master/simfiles/configs/example.py) for simple and advanced configuration examples.

## Installation:
 - Download via web UI, or `git clone https://github.com/kyleaoman/simfiles.git`
 - Install dependencies if necessary (see [`setup.py`](https://github.com/kyleaoman/simfiles/blob/master/setup.py)), some may be found in [other repositories by kyleaoman](https://github.com/kyleaoman?tab=repositories).
 - Global install (Linux): 
   - cd to directory with [`setup.py`](https://github.com/kyleaoman/simfiles/blob/master/setup.py)
   - run `sudo pip install -e .` (`-e` installs via symlink, so pulling repository will do a 'live' update of the installation)
 - User install (Linux):
   - cd to directory with [`setup.py`](https://github.com/kyleaoman/simfiles/blob/master/setup.py)
   - ensure `~/lib/python2.7/site-packages` or similar is on your `PYTHONPATH` (e.g. `echo $PYTHONPATH`), if not, add it (perhaps in `.bash_profile` or similar)
   - run `pip install --prefix ~ -e .` (`-e` installs via symlink, so pulling repository will do a 'live' update of the installation)
 - cd to a directory outside the module and launch `python`; you should be able to do `from simfiles import SimFiles`

## Usage:

```python
from simfiles import SimFiles
F = SimFiles(snap_id, configfile=None)
```

Initializes a `SimFiles` object and reads configuration.
 - `snap_id`: a unique identifier for a snapshot in the simulation dataset defined in the configfile; see [example configuration](https://github.com/kyleaoman/simfiles/blob/master/simfiles/configs/example.py) for details.
 - `configfile`: path to a configfile; see [example configuration](https://github.com/kyleaoman/simfiles/blob/master/simfiles/configs/example.py) for a description of what should be defined in this file.
 
```python
F.load(keys=None, filetype=None)
```

Loads the keys specified according to the configuration. Loaded keys can be accessed as `dict` entries (`F['keyname']`) or as attributes (`F.keyname`). Loading the same key twice should be avoided, doing so will result in a warning.
 - `keys`: a tuple of keys, as defined via the `extractors` in the configfile.
 - `filetype`: this will force the keys to be read from a specific filetype; use `None` to use the defaults as defined in the corresponding `extractors`. Prefer the default unless absolutely necessary.

```python
F.fields(keytype='all')
```

Returns a list of keys known to the object, which may be filtered to a specific key type.
 - `keytype`: keytype (as defined in `extractors`) to return, default is all known keys.
