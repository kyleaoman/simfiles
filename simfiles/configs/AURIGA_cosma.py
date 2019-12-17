from os.path import join
from simfiles._setup_cfg import snapshots, extractor, extractors
from collections import namedtuple
from astropy import units as U, constants as C
from itertools import product
from Hdecompose.BlitzRosolowsky2006 import molecular_frac
from Hdecompose.SpringelHernquist2003 import auriga_correct_neutral_frac \
    as correct_neutral_frac
import numpy as np

# define snapshot unique id tuple format
snap_id = namedtuple('snap_id', ['level', 'phys', 'halo', 'snap'])

path_base = '/cosma5/data/Gigagalaxy/pakmorr/'

for level, halo, phys, snapnum in product(
        range(3, 5),
        ['{:d}'.format(i) for i in range(1, 31)]
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

# column order for abundance tables
elements = ['H', 'He', 'C', 'N', 'O', 'Ne', 'Mg', 'Si', 'Fe']

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
    hpath='/Subhalo/SubhaloGrNr',
    attr=None,
    convert=lambda vals, raw, path, fname, hpath:
    raw,
    units=U.dimensionless_unscaled,
    unit_convert=None
)

# sgns
extractors['sgns'] = extractor(
    keytype='group',
    filetype='group',
    dependencies=('firstsub', ),
    hpath='/Subhalo/SubhaloGrNr',
    attr=None,
    convert=lambda vals, raw, path, fname, hpath:
    np.arange(raw.size) - vals.firstsub[raw],
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
    hpath='/Subhalo/SubhaloPos',
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

# msubfind
extractors['msubfind'] = extractor(
    keytype='group',
    filetype='group',
    dependencies=('h'),
    hpath='/Subhalo/SubhaloMass',
    attr=None,
    convert=lambda vals, raw, path, fname, hpath:
    raw * 1E10 / vals.h,
    units=U.solMass,
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
        dependencies=('h', ),
        hpath='/Subhalo/SubhaloMassType',
        attr=None,
        convert=subval(ptype),
        units=U.solMass,
        unit_convert=None
    )

# vmax
extractors['vmax'] = extractor(
    keytype='group',
    filetype='group',
    dependencies=('a', ),
    hpath='/Subhalo/SubhaloVmax',
    attr=None,
    convert=lambda vals, raw, path, fname, hpath:
    raw * np.sqrt(vals.a),
    units=U.km * U.s ** -1,
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

# sft_s
extractors['sft_s'] = extractor(
    keytype='particle_s',
    filetype='snapshot',
    dependencies=tuple(),
    hpath='/PartType4/GFM_StellarFormationTime',
    attr=None,
    convert=lambda vals, raw, path, fname, hpath:
    raw,
    units=U.dimensionless_unscaled,
    unit_convert=None
)

# sfr_g
extractors['sfr_g'] = extractor(
    keytype='particle_g',
    filetype='snapshot',
    dependencies=tuple(),
    hpath='/PartType0/StarFormationRate',
    attr=None,
    convert=lambda vals, raw, path, fname, hpath:
    raw,
    units=U.solMass * U.yr ** -1,
    unit_convert=None
)

# rho_g
extractors['rho_g'] = extractor(
    keytype='particle_g',
    filetype='snapshot',
    dependencies=('h', 'a'),
    hpath='/PartType0/Density',
    attr=None,
    convert=lambda vals, raw, path, fname, hpath:
    raw * 1E10 * np.power(vals.h, 2) / vals.a,
    units=U.solMass / U.Mpc ** 3,
    unit_convert=U.solMass / U.kpc ** 3
)

# xe_g
extractors['xe_g'] = extractor(
    keytype='particle_g',
    filetype='snapshot',
    dependencies=tuple(),
    hpath='/PartType0/ElectronAbundance',
    attr=None,
    convert=lambda vals, raw, path, fname, hpath:
    raw,
    units=U.dimensionless_unscaled,
    unit_convert=None
)

# u_g
extractors['u_g'] = extractor(
    keytype='particle_g',
    filetype='snapshot',
    dependencies=tuple(),
    hpath='/PartType0/InternalEnergy',
    attr=None,
    convert=lambda vals, raw, path, fname, hpath:
    raw,
    units=U.km ** 2 / U.s ** 2,
    unit_convert=None
)

# T_g
extractors['T_g'] = extractor(
    keytype='particle_g',
    filetype='snapshot',
    dependencies=('xe_g', 'fH_g'),
    hpath='/PartType0/InternalEnergy',
    attr=None,
    convert=lambda vals, raw, path, fname, hpath:
    (5 / 3 - 1) * raw  * 1E10 * 4 * C.m_p.to(U.g).value
    / C.k_B.to(U.erg / U.K).value
    / (1 + 3 * vals.fH_g + 4 * vals.fH_g * vals.xe_g),
    units=U.K,
    unit_convert=None
)


def subval(ie):
    return lambda vals, raw, path, fname, hpath: \
        raw[:, ie]


# f*_g
for ie, e in enumerate(elements):
    extractors['f{:s}_g'.format(e)] = extractor(
        keytype='particle_g',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/PartType0/GFM_Metals',
        attr=None,
        convert=subval(ie),
        units=U.dimensionless_unscaled,
        unit_convert=None
    )

# mHneutral_g
extractors['mHneutral_g'] = extractor(
    keytype='particle_g',
    filetype='snapshot',
    dependencies=('m_g', 'fH_g', 'sfr_g', 'u_g'),
    hpath='/PartType0/NeutralHydrogenAbundance',
    attr=None,
    convert=lambda vals, raw, path, fname, hpath:
    vals.m_g.to(U.solMass).value * vals.fH_g * correct_neutral_frac(
        raw,
        vals.sfr_g,
        vals.u_g
    ),
    units=U.solMass,
    unit_convert=None
)

# mHI_g
extractors['mHI_g'] = extractor(
    keytype='particle_g',
    filetype='snapshot',
    dependencies=('m_g', 'fH_g', 'sfr_g', 'u_g', 'T_g', 'rho_g'),
    hpath='/PartType0/NeutralHydrogenAbundance',
    attr=None,
    convert=lambda vals, raw, path, fname, hpath:
    vals.m_g.to(U.solMass).value * vals.fH_g *
    correct_neutral_frac(
        raw,
        vals.sfr_g,
        vals.u_g
    ) * (1 - molecular_frac(
        vals.T_g,
        vals.rho_g,
        mu=1.22,
        Auriga_corrections=True,
        SFR=vals.sfr_g,
        fNeutral=correct_neutral_frac(
            raw,
            vals.sfr_g,
            vals.u_g
        )
    )),
    units=U.solMass,
    unit_convert=None
)

# -----------------------------------------------------------------------------------------------------
