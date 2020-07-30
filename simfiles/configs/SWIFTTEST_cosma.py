from simfiles._setup_cfg import snapshots, extractor, extractors
from collections import namedtuple
from astropy import units as U
from itertools import product
import numpy as np
from simfiles._hdf5_io import hdf5_get

# define snapshot unique id tuple format
snap_id = namedtuple('snap_id', [])

path_base = '/snap7/scratch/dp004/dc-borr1/softening-runs/AN_REF'

# next line defines a snapshot by its id and specifies where to find its
# files
snapshots[snap_id()] = {
    'snapshot': (path_base, 'eagle_0028'),  # omit .X.hdf5
}

# define a mnemonic suffix for each particle type in EAGLE
T = {
    'g': '0',
    'dm': '1',
    's': '4',
}

# SWIFT gives the conversion to cgs as an attribute, use this concisely
def to_cgs(path, fname, hpath):
    return hdf5_get(
        path,
        fname,
        hpath,
        attr='Conversion factor to physical CGS '
        '(including cosmological corrections)'
    )
 
# a
extractors['a'] = extractor(
    keytype='header',
    filetype='snapshot',
    dependencies=tuple(),
    hpath='/Header',
    attr='Scale-factor',
    convert=lambda vals, raw, path, fname, hpath:
        raw,
    units=U.dimensionless_unscaled,
    unit_convert=None
)

# h
extractors['h'] = extractor(
    keytype='header',
    filetype='snapshot',
    dependencies=tuple(),
    hpath='/Cosmology',
    attr='h',
    convert=lambda vals, raw, path, fname, hpath:
        raw,
    units=U.dimensionless_unscaled,
    unit_convert=None
)

# Lbox
extractors['Lbox'] = extractor(
    keytype='header',
    filetype='snapshot',
    dependencies=tuple(),
    hpath='/Header',
    attr='BoxSize',
    convert=lambda vals, raw, path, fname, hpath:
        raw[0],  # assume equal length in each dimension
    units=U.Mpc,
    unit_convert=None
)

# xyz_*
for ptype in T.keys():
    extractors['xyz_' + ptype] = extractor(
        keytype='particle_'+ptype,
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/PartType' + T[ptype] + '/Coordinates',
        attr=None,
        convert=lambda vals, raw, path, fname, hpath:
            raw * to_cgs(path, fname, hpath),
        units=U.cm,
        unit_convert=U.Mpc
    )

# vxyz_*
for ptype in T.keys():
    extractors['vxyz_' + ptype] = extractor(
        keytype='particle_'+ptype,
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/PartType' + T[ptype] + '/Velocities',
        attr=None,
        convert=lambda vals, raw, path, fname, hpath:
            raw * to_cgs(path, fname, hpath),
        units=U.cm * U.s ** -1,
        unit_convert=U.km * U.s ** -1
    )
