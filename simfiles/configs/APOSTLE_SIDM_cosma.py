from simfiles._setup_cfg import snapshots, extractor, extractors
from collections import namedtuple
from astropy import units as U
from itertools import product
import numpy as np
from astropy.cosmology import FlatLambdaCDM
from .._hdf5_io import hdf5_get
from Hdecompose.atomic_frac import atomic_frac

# annoying redshift text suffixes for EAGLE/APOSTLE hdf5 files
suffix = [
    '000_z020p000', '001_z019p503', '002_z019p017', '003_z018p543',
    '004_z018p080', '005_z017p628', '006_z017p187', '007_z016p756',
    '008_z016p335', '009_z015p925', '010_z015p524', '011_z015p132',
    '012_z014p750', '013_z014p377', '014_z014p013', '015_z013p657',
    '016_z013p310', '017_z012p971', '018_z012p640', '019_z012p317',
    '020_z012p002', '021_z011p694', '022_z011p393', '023_z011p099',
    '024_z010p813', '025_z010p533', '026_z010p260', '027_z009p993',
    '028_z009p733', '029_z009p478', '030_z009p230', '031_z008p988',
    '032_z008p751', '033_z008p520', '034_z008p295', '035_z008p075',
    '036_z007p860', '037_z007p650', '038_z007p445', '039_z007p245',
    '040_z007p050', '041_z006p859', '042_z006p673', '043_z006p491',
    '044_z006p314', '045_z006p140', '046_z005p971', '047_z005p806',
    '048_z005p645', '049_z005p487', '050_z005p334', '051_z005p184',
    '052_z005p037', '053_z004p894', '054_z004p755', '055_z004p618',
    '056_z004p485', '057_z004p355', '058_z004p228', '059_z004p105',
    '060_z003p984', '061_z003p866', '062_z003p750', '063_z003p638',
    '064_z003p528', '065_z003p421', '066_z003p316', '067_z003p214',
    '068_z003p114', '069_z003p017', '070_z002p921', '071_z002p828',
    '072_z002p738', '073_z002p649', '074_z002p563', '075_z002p478',
    '076_z002p396', '077_z002p316', '078_z002p237', '079_z002p160',
    '080_z002p085', '081_z002p012', '082_z001p941', '083_z001p871',
    '084_z001p803', '085_z001p737', '086_z001p672', '087_z001p609',
    '088_z001p547', '089_z001p487', '090_z001p428', '091_z001p370',
    '092_z001p314', '093_z001p259', '094_z001p206', '095_z001p154',
    '096_z001p103', '097_z001p053', '098_z001p004', '099_z000p957',
    '100_z000p910', '101_z000p865', '102_z000p821', '103_z000p778',
    '104_z000p736', '105_z000p695', '106_z000p654', '107_z000p615',
    '108_z000p577', '109_z000p540', '110_z000p503', '111_z000p468',
    '112_z000p433', '113_z000p399', '114_z000p366', '115_z000p333',
    '116_z000p302', '117_z000p271', '118_z000p241', '119_z000p211',
    '120_z000p183', '121_z000p155', '122_z000p127', '123_z000p101',
    '124_z000p075', '125_z000p049', '126_z000p024', '127_z000p000'
]

# define snapshot unique id tuple format
snap_id = namedtuple('snap_id', ['res', 'phys', 'vol', 'snap'])

path_base = '/cosma5/data/Eagle/Apostle_SIDM/'
res_str = {1: 'HR', 2: 'MR', 3: 'LR'}
vol_str = {1: 'V1', 2: 'V2', 3: 'V3', 4: 'V4', 5: 'V5', 6: 'V6',
           7: 'S1', 8: 'S2', 9: 'S3', 10: 'S4', 11: 'S5', 12: 'S6'}
phys_str = {'hydro': 'REF'}

for res, vol, phys, snapnum, crosssec in product(
        range(1, 4), range(1, 13), ('hydro', 'DMO'), range(128), ('10gcm2', )):

    path_prefix = path_base + vol_str[vol] + '_' + res_str[res] + '_' + \
        phys_str[phys] + '_XS_' + crosssec + '/Apos/data/'

    group_path = path_prefix + '/groups_' + suffix[snapnum]
    group_file = 'eagle_subfind_tab_' + suffix[snapnum]
    particle_path = path_prefix + '/particledata_' + suffix[snapnum]
    particle_file = 'eagle_subfind_particles_' + suffix[snapnum]
    snapshot_path = path_prefix + '/snapshot_' + suffix[snapnum]
    snapshot_file = 'snap_' + suffix[snapnum]

    if (res == 1) and (vol == 1) and (phys == 'hydro') and \
       (crosssec == '10gcm2'):
        pass
    else:
        continue

    # next line defines a snapshot by its id and specifies where to find its
    # files
    snapshots[snap_id(res=res, phys=phys, vol=vol, snap=snapnum)] = {
        'group': (group_path, group_file),  # omit .X.hdf5
        'particle': (particle_path, particle_file),  # omit .X.hdf5
        'snapshot': (snapshot_path, snapshot_file),  # omit .X.hdf5
    }

# define a mnemonic suffix for each particle type in EAGLE/APOSTLE
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

# define abbreviations for annoying softening length suffixes in EAGLE/APOSTLE
softstrings = {
    'g': 'Gas',
    'dm': 'Halo',
    'b2': 'Disk',
    'b3': 'Bulge',
    's': 'Stars',
    'bh': 'Bndry'
}


# many EAGLE/APOSTLE fields specify exponents for h and a; use this function to
# apply them concisely
def h_a_powers(vals, path, fname, hpath):
    return np.power(
        vals.h, hdf5_get(path, fname, hpath, attr='h-scale-exponent')
    ) * \
        np.power(
            vals.a, hdf5_get(path, fname, hpath, attr='aexp-scale-exponent')
        )


# convenience function to get molecular weight once other parameters are
# defined
def mu(vals):
    return 1. / (vals.fH + .25 * vals.fHe)


# convenience function to get an astropy cosmology utility
def cosmo(vals): \
    return FlatLambdaCDM(
        H0=vals.h * 100. * U.km * U.s ** -1 * U.Mpc ** -1,
        Om0=vals.Omega0,
        Ob0=vals.Omegab
    )


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

# code_to_g
extractors['code_to_g'] = extractor(
    keytype='header',
    filetype='snapshot',
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
    filetype='snapshot',
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
    filetype='snapshot',
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
    filetype='snapshot',
    dependencies=('code_to_cm', 'h'),
    hpath='/Header',
    attr='BoxSize',
    convert=lambda vals, raw, path, fname, hpath:
    raw / vals.h * vals.code_to_cm,
    units=U.cm,
    unit_convert=U.kpc
)

# G
extractors['G'] = extractor(
    keytype='header',
    filetype='snapshot',
    dependencies=tuple(),
    hpath='/Constants',
    attr='GRAVITY',
    convert=lambda vals, raw, path, fname, hpath:
    raw,
    units=U.cm ** 3 * U.s ** -2 * U.g ** -1,
    unit_convert=U.km ** 2 * U.s ** -2 * U.Msun ** -1 * U.kpc
)

# proton_mass
extractors['proton_mass'] = extractor(
    keytype='header',
    filetype='snapshot',
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
    filetype='snapshot',
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
    filetype='snapshot',
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
    filetype='snapshot',
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
    filetype='snapshot',
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
    filetype='snapshot',
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
    filetype='snapshot',
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
    filetype='snapshot',
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
    filetype='snapshot',
    dependencies=('code_to_g', 'h'),
    hpath='/Header',
    attr='MassTable',
    convert=lambda vals, raw, path, fname, hpath:
    raw[1] / vals.h * vals.code_to_g,
    units=U.g,
    unit_convert=None
)


# this helper function allows correct definition of lambda when looping over
# ptype in defining eps_*
def subval(s):
    return lambda vals, raw, path, fname, hpath:\
        min(
            raw * vals.a / vals.h,
            vals[s].to(U.cm) / U.cm / vals.code_to_cm,
        ) * vals.code_to_cm


# eps_*
for ptype in T.keys():
    extractors['eps_' + ptype] = extractor(
        keytype='header',
        filetype='snapshot',
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
        filetype='snapshot',
        dependencies=('code_to_cm', 'h', 'a'),
        hpath='/RuntimePars',
        attr='Softening' + softstrings[ptype] + 'MaxPhys',
        convert=lambda vals, raw, path, fname, hpath:
        raw / vals.h * vals.code_to_cm,
        units=U.cm,
        unit_convert=U.kpc
    )

# contamination
extractors['contamination'] = extractor(
    keytype='fofgroup',
    filetype='group',
    dependencies=tuple(),
    hpath='/FOF/ContaminationCount',
    attr=None,
    convert=lambda vals, raw, path, fname, hpath:
    raw,
    units=U.dimensionless_unscaled,
    unit_convert=None
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


def subval(s):
    return lambda vals, raw, path, fname, hpath: \
        raw[:, int(T[s])]


# sl_*
for ptype in T.keys():
    extractors['sl_'+ptype] = extractor(
        keytype='group',
        filetype='group',
        dependencies=tuple(),
        hpath='/Subhalo/SubLengthType',
        attr=None,
        convert=subval(ptype),
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
    units=U.dimensionless_unscaled,
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
    units=U.dimensionless_unscaled,
    unit_convert=None
)


def subval(s):
    return lambda vals, raw, path, fname, hpath: \
        raw[:, int(T[s])] * h_a_powers(vals, path, fname, hpath) * \
        vals.code_to_g


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
        hpath='/PartType' + T[ptype] + '/Velocities',
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

# m_g, m_b2, m_b3, m_s, m_bh
for ptype in ['g', 'b2', 'b3', 's', 'bh']:
    extractors['m_' + ptype] = extractor(
        keytype='particle_'+ptype,
        filetype='particle',
        dependencies=('code_to_g', 'h', 'a'),
        hpath='/PartType' + T[ptype] + '/Masses',
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
    (cosmo(vals).age(vals.redshift)
     - cosmo(vals).age(1. / raw - 1.)).to(U.yr) / U.yr,
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
    hpath='/PartType0/Masses',
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

# -----------------------------------------------------------------------------------------------------
