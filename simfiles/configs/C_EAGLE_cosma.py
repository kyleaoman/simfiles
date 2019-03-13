from simfiles._setup_cfg import snapshots, extractor, extractors
from collections import namedtuple
from astropy import units as U
from itertools import product
import numpy as np
from astropy.cosmology import FlatLambdaCDM
from ._hdf5_io import hdf5_get
from Hdecompose.atomic_frac import atomic_frac
from Hdecompose.RahmatiEtal2013 import neutral_frac

# annoying redshift text suffixes for C-EAGLE hdf5 files
suffix = [
    '000_z014p003', '001_z006p772', '002_z004p614',
    '003_z003p512', '004_z002p825', '005_z002p348',
    '006_z001p993', '007_z001p716', '008_z001p493',
    '009_z001p308', '010_z001p151', '011_z001p017',
    '012_z000p899', '013_z000p795', '014_z000p703',
    '015_z000p619', '016_z000p543', '017_z000p474',
    '018_z000p411', '019_z000p366', '020_z000p352',
    '021_z000p297', '022_z000p247', '023_z000p199',
    '024_z000p155', '025_z000p113', '026_z000p101',
    '027_z000p073', '028_z000p036', '029_z000p000',
]

# define snapshot unique id tuple format
snap_id = namedtuple('snap_id', ['sample', 'ce', 'model', 'snap'])

path_base = '/cosma5/data/dp004/C-EAGLE/'

volumes = \
    {
        'Complete_Sample': {
            'CE_{0:02d}'.format(k): ['HYDRO'] for k in range(30)
        },
        'Hydrangea': {
            'CE-{0:.0f}'.format(k): ['DM', 'HYDRO'] for k in [
                0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 18,
                21, 22, 24, 25, 28, 29
            ]
        }
    }

# CompleteSample partial outputs: 6 7 12 13 14 15 16 18 22 24 25 28 29 (these
#                                 symlinked to hydrangea)
# partial output = groups 006 008 011 017 026 029; particles 026 029;
#                  snaps 017 026 029
# Hydrangea particles 026 029; groups 006 008 011 017 026 029;
#           snaps 011 026 029
# Hydrangea CE-0, CE-1 have extra data
# Hydrangea DMO exists for: CE >= 6 (029 data);
#                           6, 7 and 12 have group data for 008 011 017 029
# DMO simulations don't seem to match resolution (?)


box_list = [(sample, ce, model)
            for sample, v in volumes.items()
            for ce, vv in v.items()
            for model in vv]

for (sample, ce, model), snap in product(box_list, range(30)):

    if sample == 'Complete_Sample':
        path_prefix = '/'.join((path_base, sample, ce, 'data'))
    elif sample == 'Hydrangea':
        path_prefix = '/'.join((path_base, sample, ce, model, 'data'))

    group_path = path_prefix + '/groups_' + suffix[snap]
    group_file = 'eagle_subfind_tab_' + suffix[snap]
    particle_path = path_prefix + '/particledata_' + suffix[snap]
    particle_file = 'eagle_subfind_particles_' + suffix[snap]
    snapshot_path = path_prefix + '/snapshot_' + suffix[snap]
    snapshot_file = 'snap_' + suffix[snap]

    # next line defines a snapshot by its id and specifies where to find files
    snapshots[snap_id(sample=sample, ce=ce, model=model, snap=snap)] = {
        'group': (group_path, group_file),  # omit .X.hdf5
        'particle': (particle_path, particle_file),  # omit .X.hdf5
        'snapshot': (snapshot_path, snapshot_file),  # omit .X.hdf5
    }

# define a mnemonic suffix for each particle type in EAGLE
T = {
    'g': '0',
    'dm': '1',
    'b2': '2',
    'b3': '3',
    's': '4',
    'bh': '5'
}

# define abbreviations for elemental abundance fields
elements = {
    'H': 'ElementAbundance/Hydrogen',
    'He': 'ElementAbundance/Helium',
    'C': 'ElementAbundance/Carbon',
    'Mg': 'ElementAbundance/Magnesium',
    'Fe': 'ElementAbundance/Iron',
    'Ne': 'ElementAbundance/Neon',
    'N': 'ElementAbundance/Nitrogen',
    'O': 'ElementAbundance/Oxygen',
    'Si': 'ElementAbundance/Silicon',
    'Z': 'Metallicity'
}

# define abbreviations for annoying softening length suffixes in EAGLE
softstrings = {
    'g': 'Gas',
    'dm': 'Halo',
    'b2': 'Disk',
    'b3': 'Bulge',
    's': 'Stars',
    'bh': 'Bndry'
}


# many fields specify exponents for h and a; use this to apply them concisely
def h_a_powers(vals, path, fname, hpath):
    args = (path, fname, hpath)
    return np.power(vals.h, hdf5_get(*args, attr='h-scale-exponent')) * \
        np.power(vals.a, hdf5_get(*args, attr='aexp-scale-exponent'))


# function to get molecular weight once other parameters are defined
def mu(vals):
    return 1. / (vals.fH + .25 * vals.fHe)


# convenience function to get cosmology utility
def cosmo(vals):
    return FlatLambdaCDM(
        H0=vals.h * 100. * U.km * U.s ** -1 * U.Mpc ** -1,
        Om0=vals.Omega0,
        Ob0=vals.Omegab
    )


# a
extractors['a'] = extractor(
    keytype='header',
    filetype='particle',
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
    filetype='particle',
    dependencies=tuple(),
    hpath='/Header',
    attr='HubbleParam',
    convert=lambda vals, raw, path, fname, hpath:
    raw,
    units=U.dimensionless_unscaled,
    unit_convert=None
)

# code_to_g
extractors['code_to_g'] = extractor(
    keytype='header',
    filetype='particle',
    dependencies=tuple(),
    hpath='/Units',
    attr='UnitMass_in_g',
    convert=lambda vals, raw, path, fname, hpath:
        raw,
    units=U.dimensionless_unscaled,
    unit_convert=None
)

# code_to_cm
extractors['code_to_cm'] = extractor(
    keytype='header',
    filetype='particle',
    dependencies=tuple(),
    hpath='/Units',
    attr='UnitLength_in_cm',
    convert=lambda vals, raw, path, fname, hpath:
        raw,
    units=U.dimensionless_unscaled,
    unit_convert=None
)

# code_to_cm_s
extractors['code_to_cm_s'] = extractor(
    keytype='header',
    filetype='particle',
    dependencies=tuple(),
    hpath='/Units',
    attr='UnitVelocity_in_cm_per_s',
    convert=lambda vals, raw, path, fname, hpath:
        raw,
    units=U.dimensionless_unscaled,
    unit_convert=None
)

# Lbox
extractors['Lbox'] = extractor(
    keytype='header',
    filetype='particle',
    dependencies=('code_to_cm', 'h'),
    hpath='/Header',
    attr='BoxSize',
    convert=lambda vals, raw, path, fname, hpath:
        raw / vals.h * vals.code_to_cm,
    units=U.cm,
    unit_convert=U.kpc
)

# proton_mass
extractors['proton_mass'] = extractor(
    keytype='header',
    filetype='particle',
    dependencies=tuple(),
    hpath='/Constants',
    attr='PROTONMASS',
    convert=lambda vals, raw, path, fname, hpath:
        raw,
    units=U.g,
    unit_convert=None
)

# redshift
extractors['redshift'] = extractor(
    keytype='header',
    filetype='particle',
    dependencies=tuple(),
    hpath='/Header',
    attr='Redshift',
    convert=lambda vals, raw, path, fname, hpath:
        raw,
    units=U.dimensionless_unscaled,
    unit_convert=None
)

# fH
extractors['fH'] = extractor(
    keytype='header',
    filetype='particle',
    dependencies=tuple(),
    hpath='/RuntimePars',
    attr='InitAbundance_Hydrogen',
    convert=lambda vals, raw, path, fname, hpath:
        raw,
    units=U.dimensionless_unscaled,
    unit_convert=None
)

# fHe
extractors['fHe'] = extractor(
    keytype='header',
    filetype='particle',
    dependencies=tuple(),
    hpath='/RuntimePars',
    attr='InitAbundance_Helium',
    convert=lambda vals, raw, path, fname, hpath:
        raw,
    units=U.dimensionless_unscaled,
    unit_convert=None
)

# Omega0
extractors['Omega0'] = extractor(
    keytype='header',
    filetype='particle',
    dependencies=tuple(),
    hpath='/Header',
    attr='Omega0',
    convert=lambda vals, raw, path, fname, hpath:
        raw,
    units=U.dimensionless_unscaled,
    unit_convert=None
)

# Omegab
extractors['Omegab'] = extractor(
    keytype='header',
    filetype='particle',
    dependencies=tuple(),
    hpath='/Header',
    attr='OmegaBaryon',
    convert=lambda vals, raw, path, fname, hpath:
        raw,
    units=U.dimensionless_unscaled,
    unit_convert=None
)

# gamma
extractors['gamma'] = extractor(
    keytype='header',
    filetype='particle',
    dependencies=tuple(),
    hpath='/RuntimePars',
    attr='EOS_Jeans_GammaEffective',
    convert=lambda vals, raw, path, fname, hpath:
        raw,
    units=U.dimensionless_unscaled,
    unit_convert=None
)

# T0
extractors['T0'] = extractor(
    keytype='header',
    filetype='particle',
    dependencies=tuple(),
    hpath='/RuntimePars',
    attr='EOS_Jeans_TempNorm_K',
    convert=lambda vals, raw, path, fname, hpath:
        raw,
    units=U.K,
    unit_convert=None
)

# p_mass
extractors['p_mass'] = extractor(
    keytype='header',
    filetype='particle',
    dependencies=('code_to_g', 'h'),
    hpath='/Header',
    attr='MassTable',
    convert=lambda vals, raw, path, fname, hpath:
        raw[1] / vals.h * vals.code_to_g,
    units=U.g,
    unit_convert=None
)


def subval(s):
    return lambda vals, raw, path, fname, hpath: \
        min(
            raw * vals.a / vals.h,
            vals[s].to(U.cm) / U.cm / vals.code_to_cm
        ) * vals.code_to_cm


# eps_*
for ptype in T.keys():
    extractors['eps_' + ptype] = extractor(
        keytype='header',
        filetype='particle',
        dependencies=('code_to_cm', 'h', 'a', 'eps_maxphys_' + ptype),
        hpath='/RuntimePars',
        attr='Softening' + softstrings[ptype],
        convert=subval('eps_maxphys_' + ptype),
        units=U.cm,
        unit_convert=U.kpc
    )

# eps_maxphys_*
for ptype in T.keys():
    extractors['eps_maxphys_' + ptype] = extractor(
        keytype='header',
        filetype='particle',
        dependencies=('code_to_cm', 'h', 'a'),
        hpath='/RuntimePars',
        attr='Softening' + softstrings[ptype] + 'MaxPhys',
        convert=lambda vals, raw, path, fname, hpath:
            raw / vals.h * vals.code_to_cm,
        units=U.cm,
        unit_convert=U.kpc
    )

# nsubhalos
extractors['nsubhalos'] = extractor(
    keytype='fofgroup',
    filetype='group',
    dependencies=tuple(),
    hpath='/FOF/NumOfSubhalos',
    attr=None,
    convert=lambda vals, raw, path, fname, hpath:
        raw,
    units=U.dimensionless_unscaled,
    unit_convert=None
)

# fshid
extractors['fshid'] = extractor(
        keytype='fofgroup',
        filetype='group',
        dependencies=tuple(),
        hpath='/FOF/FirstSubhaloID',
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
    hpath='/Subhalo/GroupNumber',
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
    dependencies=tuple(),
    hpath='/Subhalo/SubGroupNumber',
    attr=None,
    convert=lambda vals, raw, path, fname, hpath:
        raw,
    units=U.dimensionless_unscaled,
    unit_convert=None
)

# cops
extractors['cops'] = extractor(
    keytype='group',
    filetype='group',
    dependencies=('code_to_cm', 'h', 'a'),
    hpath='/Subhalo/CentreOfPotential',
    attr=None,
    convert=lambda vals, raw, path, fname, hpath:
        raw * h_a_powers(vals, path, fname, hpath) * vals.code_to_cm,
    units=U.cm,
    unit_convert=U.kpc
)

# vcents
extractors['vcents'] = extractor(
    keytype='group',
    filetype='group',
    dependencies=('code_to_cm_s', 'h', 'a'),
    hpath='/Subhalo/Velocity',
    attr=None,
    convert=lambda vals, raw, path, fname, hpath:
        raw * h_a_powers(vals, path, fname, hpath) * vals.code_to_cm_s,
    units=U.cm * U.s ** -1,
    unit_convert=None
)

# nID
extractors['nID'] = extractor(
    keytype='group',
    filetype='group',
    dependencies=tuple(),
    hpath='/Subhalo/SubLength',
    attr=None,
    convert=lambda vals, raw, path, fname, hpath:
        raw,
    units=None,
    unit_convert=None
)

# offID
extractors['offID'] = extractor(
    keytype='group',
    filetype='group',
    dependencies=tuple(),
    hpath='/Subhalo/SubOffset',
    attr=None,
    convert=lambda vals, raw, path, fname, hpath:
        raw,
    units=None,
    unit_convert=None
)

# msubfind
extractors['msubfind'] = extractor(
    keytype='group',
    filetype='group',
    dependencies=('code_to_g', 'h', 'a'),
    hpath='/Subhalo/MassType',
    attr=None,
    convert=lambda vals, raw, path, fname, hpath:
        raw * h_a_powers(vals, path, fname, hpath) * vals.code_to_g,
    units=U.g,
    unit_convert=U.solMass
)


def subval(s):
    return lambda vals, raw, path, fname, hpath: \
        raw[:, int(T[s])] * h_a_powers(vals, path, fname, hpath) \
        * vals.code_to_g


# msubfind_*
for ptype in T.keys():
    extractors['msubfind_' + ptype] = extractor(
        keytype='group',
        filetype='group',
        dependencies=('code_to_g', 'h', 'a'),
        hpath='/Subhalo/MassType',
        attr=None,
        convert=subval(ptype),
        units=U.g,
        unit_convert=U.solMass
    )

# nfof
extractors['nfof'] = extractor(
    keytype='header',
    filetype='group',
    dependencies=tuple(),
    hpath='/FOF',
    attr='TotNgroups',
    convert=lambda vals, raw, path, fname, hpath:
        raw,
    units=U.dimensionless_unscaled,
    unit_convert=None
)

# M200
extractors['M200'] = extractor(
    keytype='fofgroup',
    filetype='group',
    dependencies=('code_to_g', 'h', 'a'),
    hpath='/FOF/Group_M_Crit200',
    attr=None,
    convert=lambda vals, raw, path, fname, hpath:
        raw * h_a_powers(vals, path, fname, hpath) * vals.code_to_g,
    units=U.g,
    unit_convert=U.solMass
)

# R200
extractors['R200'] = extractor(
    keytype='fofgroup',
    filetype='group',
    dependencies=('code_to_cm', 'h', 'a'),
    hpath='/FOF/Group_R_Crit200',
    attr=None,
    convert=lambda vals, raw, path, fname, hpath:
        raw * h_a_powers(vals, path, fname, hpath) * vals.code_to_cm,
    units=U.cm,
    unit_convert=U.kpc
)

# ids
extractors['ids'] = extractor(
    keytype='idgroup',
    filetype='group',
    dependencies=tuple(),
    hpath='/IDs/ParticleID',
    attr=None,
    convert=lambda vals, raw, path, fname, hpath:
        raw,
    units=U.dimensionless_unscaled,
    unit_convert=None
)

# ids_*
for ptype in T.keys():
    extractors['ids_' + ptype] = extractor(
        keytype='particle_'+ptype,
        filetype='particle',
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
        filetype='particle',
        dependencies=('code_to_cm', 'h', 'a'),
        hpath='/PartType' + T[ptype] + '/Coordinates',
        attr=None,
        convert=lambda vals, raw, path, fname, hpath:
            raw * h_a_powers(vals, path, fname, hpath) * vals.code_to_cm,
        units=U.cm,
        unit_convert=U.kpc
    )

# vxyz_*
for ptype in T.keys():
    extractors['vxyz_' + ptype] = extractor(
        keytype='particle_'+ptype,
        filetype='particle',
        dependencies=('code_to_cm_s', 'h', 'a'),
        hpath='/PartType' + T[ptype] + '/Velocity',
        attr=None,
        convert=lambda vals, raw, path, fname, hpath:
            raw * h_a_powers(vals, path, fname, hpath) * vals.code_to_cm_s,
        units=U.cm * U.s ** -1,
        unit_convert=U.km * U.s ** -1
    )

# ng_*
for ptype in T.keys():
    extractors['ng_' + ptype] = extractor(
        keytype='particle_'+ptype,
        filetype='particle',
        dependencies=tuple(),
        hpath='/PartType' + T[ptype] + '/GroupNumber',
        attr=None,
        convert=lambda vals, raw, path, fname, hpath:
            raw,
        units=U.dimensionless_unscaled,
        unit_convert=None
    )

# nsg_*
for ptype in T.keys():
    extractors['nsg_'+ptype] = extractor(
        keytype='particle_'+ptype,
        filetype='particle',
        dependencies=tuple(),
        hpath='/PartType' + T[ptype] + '/SubGroupNumber',
        attr=None,
        convert=lambda vals, raw, path, fname, hpath:
            raw,
        units=U.dimensionless_unscaled,
        unit_convert=None
    )

# m_g, m_s, m_bh
for ptype in ['g', 's', 'bh']:
    extractors['m_' + ptype] = extractor(
        keytype='particle_'+ptype,
        filetype='particle',
        dependencies=('code_to_g', 'h', 'a'),
        hpath='/PartType' + T[ptype] + '/Mass',
        attr=None,
        convert=lambda vals, raw, path, fname, hpath:
            raw * h_a_powers(vals, path, fname, hpath) * vals.code_to_g,
        units=U.g,
        unit_convert=U.solMass
    )

# m_dm
extractors['m_dm'] = extractor(
    keytype='particle_dm',
    filetype='particle',
    dependencies=('p_mass',),
    hpath='/PartType1/ParticleIDs',
    attr=None,
    convert=lambda vals, raw, path, fname, hpath:
        np.ones(raw.shape, dtype=np.float) * vals.p_mass.to(U.g) / U.g,
    units=U.g,
    unit_convert=U.solMass
)

# T_g
extractors['T_g'] = extractor(
    keytype='particle_g',
    filetype='particle',
    dependencies=('h', 'a'),
    hpath='/PartType0/Temperature',
    attr=None,
    convert=lambda vals, raw, path, fname, hpath:
        raw * h_a_powers(vals, path, fname, hpath),
    units=U.K,
    unit_convert=None
)

# rho_g
extractors['rho_g'] = extractor(
    keytype='particle_g',
    filetype='particle',
    dependencies=('code_to_g', 'code_to_cm', 'h', 'a'),
    hpath='/PartType0/Density',
    attr=None,
    convert=lambda vals, raw, path, fname, hpath:
    raw * h_a_powers(vals, path, fname, hpath) * vals.code_to_g *
    np.power(vals.code_to_cm, -3),
    units=U.g * U.cm ** -3,
    unit_convert=None
)

# *abundance_g, *abundance_s, sm*abundance_g, sm*abundance_s
for ptype in ['g', 's']:
    for element in elements.keys():
        for prefix, smooth in {'sm': 'Smoothed', '': ''}.items():
            extractors[prefix+element+'abundance_'+ptype] = extractor(
                keytype='particle_'+ptype,
                filetype='particle',
                dependencies=('h', 'a'),
                hpath='/PartType'+T[ptype]+'/'+smooth+elements[element],
                attr=None,
                convert=lambda vals, raw, path, fname, hpath:
                    raw * h_a_powers(vals, path, fname, hpath),
                units=U.dimensionless_unscaled,
                unit_convert=None
            )
# SFR_g
extractors['SFR_g'] = extractor(
    keytype='particle_g',
    filetype='particle',
    dependencies=('h', 'a'),
    hpath='/PartType0/StarFormationRate',
    attr=None,
    convert=lambda vals, raw, path, fname, hpath:
        raw * h_a_powers(vals, path, fname, hpath),
    units=U.solMass / U.yr,
    unit_convert=None
)

# hsm_g, hsm_s
for ptype in ['g', 's']:
    extractors['hsm_' + ptype] = extractor(
        keytype='particle_'+ptype,
        filetype='particle',
        dependencies=('code_to_cm', 'h', 'a'),
        hpath='/PartType'+T[ptype]+'/SmoothingLength',
        attr=None,
        convert=lambda vals, raw, path, fname, hpath:
            raw * h_a_powers(vals, path, fname, hpath) * vals.code_to_cm,
        units=U.cm,
        unit_convert=U.kpc
    )

# age_s
extractors['age_s'] = extractor(
    keytype='particle_s',
    filetype='particle',
    dependencies=('h', 'Omega0', 'Omegab', 'redshift'),
    hpath='/PartType4/StellarFormationTime',
    attr=None,
    convert=lambda vals, raw, path, fname, hpath:
    (cosmo(vals).age(vals.redshift) -
     cosmo(vals).age(1. / raw - 1.)).to(U.yr) / U.yr,
    units=U.yr,
    unit_convert=U.Gyr
)

# mHI_g
extractors['mHI_g'] = extractor(
    keytype='particle_g',
    filetype='particle',
    dependencies=(
        'redshift',
        'rho_g',
        'Habundance_g',
        'proton_mass',
        'SFR_g',
        'fH',
        'fHe',
        'T_g',
        'code_to_g',
        'T0',
        'gamma'
    ),
    hpath='/PartType0/Mass',
    attr=None,
    convert=lambda vals, raw, path, fname, hpath: raw *
    vals.Habundance_g *
    h_a_powers(vals, path, fname, hpath) *
    atomic_frac(
        vals.redshift,
        vals.rho_g * vals.Habundance_g / (mu(vals) * vals.proton_mass),
        vals.T_g,
        vals.rho_g,
        vals.Habundance_g,
        onlyA1=True,
        EAGLE_corrections=True,
        SFR=vals.SFR_g,
        mu=mu(vals),
        gamma=vals.gamma,
        fH=vals.fH,
        T0=vals.T0,
    ) *
    vals.code_to_g,
    units=U.g,
    unit_convert=U.solMass
)

# mHneutral_g
extractors['mHneutral_g'] = extractor(
    keytype='particle_g',
    filetype='particle',
    dependencies=(
        'redshift',
        'rho_g',
        'Habundance_g',
        'proton_mass',
        'SFR_g',
        'fH',
        'fHe',
        'T_g',
        'code_to_g',
        'T0',
        'gamma'
    ),
    hpath='/PartType0/Mass',
    attr=None,
    convert=lambda vals, raw, path, fname, hpath: raw *
    vals.Habundance_g *
    h_a_powers(vals, path, fname, hpath) *
    neutral_frac(
        vals.redshift,
        vals.rho_g * vals.Habundance_g / (mu(vals) * vals.proton_mass),
        vals.T_g,
        onlyA1=True,
        EAGLE_corrections=True,
        SFR=vals.SFR_g,
        mu=mu(vals),
        gamma=vals.gamma,
        fH=vals.fH,
        Habundance=vals.Habundance_g,
        T0=vals.T0,
        rho=vals.rho_g,
    ) *
    vals.code_to_g,
    units=U.g,
    unit_convert=U.solMass
)

# ------------------------------------------------------------------------------
