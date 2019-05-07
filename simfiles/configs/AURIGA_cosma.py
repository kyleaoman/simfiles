from os.path import join
from simfiles._setup_cfg import snapshots, extractor, extractors
from collections import namedtuple
from astropy import units as U
from itertools import product
import numpy as np

# define snapshot unique id tuple format
snap_id = namedtuple('snap_id', ['level', 'phys', 'halo', 'snap'])

path_base = '/cosma5/data/Gigagalaxy/pakmorr/'

for level, halo, phys, snapnum in product(
        range(3, 5),
        ['{:d}'.format(i) for i in range(1, 30)]
        + ['L{:d}'.format(i) for i in range(1, 11)],
        ('DM', 'MHD'),
        range(128)
):

    uscore = '' if ((level == 3) and (halo in {'16', '24', '6'})) else '_'
    Lnew = '_new' if 'L' in halo else ''
    path_prefix = join(
        path_base,
        'level{:d}_{:s}{:s}'.format(level, phys, Lnew),
        'halo{:s}{:s}'.format(uscore, halo),
        'output'
    )

    group_path = join(path_prefix, 'groups_{:d}'.format(snapnum))
    group_file = 'fof_subhalo_tab_{:d}'.format(snapnum)
    snapshot_path = join(path_prefix, 'snapdir_{:d}'.format(snapnum))
    snapshot_file = 'snapshot_{:d}'.format(snapnum)

    if (level == 3) and not (halo in {'6', '16', '21', '23', '24', '27'}):
        continue

    # next line defines a snapshot by its id and specifies where to find its
    # files
    snapshots[snap_id(level=level, phys=phys, halo=halo, snap=snapnum)] = {
        'group': (group_path, group_file),  # omit .X.hdf5
        'snapshot': (snapshot_path, snapshot_file),  # omit .X.hdf5
    }

# define a mnemonic suffix for each particle type in AURIGA
T = {
    'g': '0',
    'dm': '1',
    'b2': '2',
    'b3': '3',
    's': '4',
    'bh': '5'
}

# a
extractors['a'] = extractor(
    keytype='header',
    filetype='snapshot',
    dependencies=tuple(),
    hpath='/Header',
    attr='Time',
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
    hpath='/Header',
    attr='HubbleParam',
    convert=lambda vals, raw, path, fname, hpath:
    raw,
    units=U.dimensionless_unscaled,
    unit_convert=None
)

# Lbox
extractors['Lbox'] = extractor(
    keytype='header',
    filetype='snapshot',
    dependencies=('h'),
    hpath='/Header',
    attr='BoxSize',
    convert=lambda vals, raw, path, fname, hpath:
    raw / vals.h,
    units=U.Mpc,
    unit_convert=U.kpc
)

# redshift
extractors['redshift'] = extractor(
    keytype='header',
    filetype='snapshot',
    dependencies=tuple(),
    hpath='/Header',
    attr='Redshift',
    convert=lambda vals, raw, path, fname, hpath:
    raw,
    units=U.dimensionless_unscaled,
    unit_convert=None
)

# Omega0
extractors['Omega0'] = extractor(
    keytype='header',
    filetype='snapshot',
    dependencies=tuple(),
    hpath='/Header',
    attr='Omega0',
    convert=lambda vals, raw, path, fname, hpath:
    raw,
    units=U.dimensionless_unscaled,
    unit_convert=None
)

# p_mass
extractors['p_mass'] = extractor(
    keytype='header',
    filetype='snapshot',
    dependencies=('h', ),
    hpath='/Header',
    attr='MassTable',
    convert=lambda vals, raw, path, fname, hpath:
    raw[1] * 1E10 / vals.h,
    units=U.solMass,
    unit_convert=None
)

# nsubhalos
extractors['nsubhalos'] = extractor(
    keytype='fofgroup',
    filetype='group',
    dependencies=tuple(),
    hpath='/Group/GroupNsubs',
    attr=None,
    convert=lambda vals, raw, path, fname, hpath:
    raw,
    units=U.dimensionless_unscaled,
    unit_convert=None
)

# firstsub
extractors['firstsub'] = extractor(
    keytype='fofgroup',
    filetype='group',
    dependencies=tuple(),
    hpath='/Group/GroupFirstSub',
    attr=None,
    convert=lambda vals, raw, path, fname, hpath:
    raw,
    units=U.dimensionless_unscaled,
    unit_convert=None
)

# gns
extractors['gns'] = extractor(
    keytype='group',
    filetype='group',
    dependencies=tuple(),
    hpath='/Subhalo/GrNr',
    attr=None,
    convert=lambda vals, raw, path, fname, hpath:
    raw,
    units=U.dimensionless_unscaled,
    unit_convert=None
)


def subval(s):
    return lambda vals, raw, path, fname, hpath: \
        raw[:, int(T[s])]


# l_*
for ptype in T.keys():
    extractors['l_'+ptype] = extractor(
        keytype='fofgroup',
        filetype='group',
        dependencies=tuple(),
        hpath='/Group/GroupLenType',
        attr=None,
        convert=subval(ptype),
        units=U.dimensionless_unscaled,
        unit_convert=None
    )

# sl_*
for ptype in T.keys():
    extractors['sl_'+ptype] = extractor(
        keytype='group',
        filetype='group',
        dependencies=tuple(),
        hpath='/Subhalo/SubhaloLenType',
        attr=None,
        convert=subval(ptype),
        units=U.dimensionless_unscaled,
        unit_convert=None
    )

# cops
extractors['cops'] = extractor(
    keytype='group',
    filetype='group',
    dependencies=('h', 'a'),
    hpath='/Subhalo/CentreOfPotential',
    attr=None,
    convert=lambda vals, raw, path, fname, hpath:
    raw * vals.a / vals.h,
    units=U.Mpc,
    unit_convert=U.kpc
)

# vcents
extractors['vcents'] = extractor(
    keytype='group',
    filetype='group',
    dependencies=('h', 'a'),
    hpath='/Subhalo/SubhaloVel',
    attr=None,
    convert=lambda vals, raw, path, fname, hpath:
    raw * np.sqrt(vals.a),
    units=U.km * U.s ** -1,
    unit_convert=None
)

# nID
extractors['nID'] = extractor(
    keytype='group',
    filetype='group',
    dependencies=tuple(),
    hpath='/Subhalo/SubhaloLen',
    attr=None,
    convert=lambda vals, raw, path, fname, hpath:
    raw,
    units=U.dimensionless_unscaled,
    unit_convert=None
)


def subval(s):
    return lambda vals, raw, path, fname, hpath: \
        raw[:, int(T[s])] * 1E10 / vals.h


# msubfind_*
for ptype in T.keys():
    extractors['msubfind_' + ptype] = extractor(
        keytype='group',
        filetype='group',
        dependencies=('h', 'a'),
        hpath='/Subhalo/SubhaloMassType',
        attr=None,
        convert=subval(ptype),
        units=U.solMass,
        unit_convert=None
    )

# nfof
extractors['nfof'] = extractor(
    keytype='header',
    filetype='group',
    dependencies=tuple(),
    hpath='/Header',
    attr='Ngroups_Total',
    convert=lambda vals, raw, path, fname, hpath:
    raw,
    units=U.dimensionless_unscaled,
    unit_convert=None
)

# nsub
extractors['nsub'] = extractor(
    keytype='header',
    filetype='group',
    dependencies=tuple(),
    hpath='/Header',
    attr='Nsubgroups_Total',
    convert=lambda vals, raw, path, fname, hpath:
    raw,
    units=U.dimensionless_unscaled,
    unit_convert=None
)

# M200
extractors['M200'] = extractor(
    keytype='fofgroup',
    filetype='group',
    dependencies=('h', ),
    hpath='/Group/Group_M_Crit200',
    attr=None,
    convert=lambda vals, raw, path, fname, hpath:
    raw * 1E10 / vals.h,
    units=U.solMass,
    unit_convert=None
)

# R200
extractors['R200'] = extractor(
    keytype='fofgroup',
    filetype='group',
    dependencies=('h', 'a'),
    hpath='/Group/Group_R_Crit200',
    attr=None,
    convert=lambda vals, raw, path, fname, hpath:
    raw * vals.a / vals.h,
    units=U.Mpc,
    unit_convert=U.kpc
)

# ids_*
for ptype in T.keys():
    extractors['ids_' + ptype] = extractor(
        keytype='particle_'+ptype,
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/PartType' + T[ptype] + '/ParticleIDs',
        attr=None,
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.dimensionless_unscaled,
        unit_convert=None
    )

# xyz_*
for ptype in T.keys():
    extractors['xyz_' + ptype] = extractor(
        keytype='particle_'+ptype,
        filetype='snapshot',
        dependencies=('h', 'a'),
        hpath='/PartType' + T[ptype] + '/Coordinates',
        attr=None,
        convert=lambda vals, raw, path, fname, hpath:
        raw * vals.a / vals.h,
        units=U.Mpc,
        unit_convert=U.kpc
    )

# vxyz_*
for ptype in T.keys():
    extractors['vxyz_' + ptype] = extractor(
        keytype='particle_'+ptype,
        filetype='snapshot',
        dependencies=('a', ),
        hpath='/PartType' + T[ptype] + '/Velocities',
        attr=None,
        convert=lambda vals, raw, path, fname, hpath:
        raw * np.sqrt(vals.a),
        units=U.km * U.s ** -1,
        unit_convert=None
    )

# m_g, m_b2, m_b3, m_s, m_bh
for ptype in ['g', 'b2', 'b3', 's', 'bh']:
    extractors['m_' + ptype] = extractor(
        keytype='particle_'+ptype,
        filetype='snapshot',
        dependencies=('h', ),
        hpath='/PartType' + T[ptype] + '/Masses',
        attr=None,
        convert=lambda vals, raw, path, fname, hpath:
        raw * 1E10 / vals.h,
        units=U.solMass,
        unit_convert=None
    )

# m_dm
extractors['m_dm'] = extractor(
    keytype='particle_dm',
    filetype='snapshot',
    dependencies=('p_mass',),
    hpath='/PartType1/ParticleIDs',
    attr=None,
    convert=lambda vals, raw, path, fname, hpath:
    np.ones(raw.shape, dtype=np.float) * vals.p_mass.to(U.solMass) / U.solMass,
    units=U.solMass,
    unit_convert=None
)

# rho_g
extractors['rho_g'] = extractor(
    keytype='particle_g',
    filetype='snapshot',
    dependencies=('code_to_g', 'code_to_cm', 'h', 'a'),
    hpath='/PartType0/Density',
    attr=None,
    convert=lambda vals, raw, path, fname, hpath:
    raw * 1E10 * np.power(vals.h, 2) / vals.a,
    units=U.solMass / U.Mpc ** 3,
    unit_convert=U.solMass / U.kpc ** 3
)

# -----------------------------------------------------------------------------------------------------
