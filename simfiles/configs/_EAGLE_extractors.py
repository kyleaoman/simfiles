import numpy as np
from astropy import units as U
from simfiles._hdf5_io import hdf5_get
from simfiles._setup_cfg import extractor

elements = ('Hydrogen', 'Helium', 'Carbon', 'Magnesium', 'Iron', 'Neon',
            'Nitrogen', 'Oxygen', 'Silicon')

softstrings = {
    0: 'Gas',
    1: 'Halo',
    2: 'Disk',
    3: 'Bulge',
    4: 'Stars',
    5: 'Bndry',
}

overdensities = ('Crit200', 'Crit2500', 'Crit500', 'Mean200', 'Mean2500',
                 'Mean500', 'TopHat200')


def h_a_powers(vals, path, fname, hpath):
    args = (path, fname, hpath)
    return np.power(
        vals.Header_attr_HubbleParam,
        hdf5_get(*args, attr='h-scale-exponent')
    ) * np.power(
        vals.Header_attr_Time,
        hdf5_get(*args, attr='aexp-scale-exponent')
    )


def generate_eagle_extractors(
        T=(0, 1, 4, 5),
        Mstring='Mass',
        Vstring='Velocity'
):

    extractors = dict()

    extractors['Header_attr_Time'] = extractor(
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

    extractors['Header_attr_HubbleParam'] = extractor(
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

    extractors['Units_attr_UnitMass_in_g'] = extractor(
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

    extractors['Units_attr_UnitLength_in_cm'] = extractor(
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

    extractors['Units_attr_UnitVelocity_in_cm_per_s'] = extractor(
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

    extractors['Header_attr_BoxSize'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=(
            'Units_attr_UnitLength_in_cm',
            'Header_attr_HubbleParam'
        ),
        hpath='/Header',
        attr='BoxSize',
        convert=lambda vals, raw, path, fname, hpath:
        raw / vals.Header_attr_HubbleParam * vals.Units_attr_UnitLength_in_cm,
        units=U.cm,
        unit_convert=U.kpc
    )

    # extractors['Constants_attr_AVOGADRO'] = None

    # extractors['Constants_attr_BOLTZMANN'] = None

    # extractors['Constants_attr_C'] = None

    # extractors['Constants_attr_CM_PER_MPC'] = None

    # extractors['Constants_attr_ELECTRONCHARGE'] = None

    # extractors['Constants_attr_ELECTRONMASS'] = None

    # extractors['Constants_attr_EV_TO_ERG'] = None

    # extractors['Constants_attr_GAMMA'] = None

    # extractors['Constants_attr_GAS_CONST'] = None

    extractors['Constants_attr_GRAVITY'] = extractor(
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

    # extractors['Constants_attr_HUBBLE'] = None

    # extractors['Constants_attr_PI'] = None

    # extractors['Constants_attr_PLANCK'] = None

    extractors['Constants_attr_PROTONMASS'] = extractor(
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

    # extractors['Constants_attr_RAD_CONST'] = None

    # extractors['Constants_attr_SEC_PER_MEGAYEAR'] = None

    # extractors['Constants_attr_SEC_PER_YEAR'] = None

    # extractors['Constants_attr_SOLAR_LUM'] = None

    # extractors['Constants_attr_SOLAR_MASS'] = None

    # extractors['Constants_attr_STEFAN'] = None

    # extractors['Constants_attr_THOMPSON'] = None

    # extractors['Constants_attr_T_CMB0'] = None

    # extractors['Constants_attr_Z_Solar'] = None

    extractors['Header_attr_Redshift'] = extractor(
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

    extractors['RuntimePars_attr_InitAbundance_Hydrogen'] = extractor(
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

    extractors['RuntimePars_attr_InitAbundance_Helium'] = extractor(
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

    extractors['Header_attr_Omega0'] = extractor(
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

    extractors['Header_attr_OmegaBaryon'] = extractor(
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

    extractors['RuntimePars_attr_EOS_Jeans_GammaEffective'] = extractor(
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

    extractors['RuntimePars_attr_EOS_Jeans_TempNorm_K'] = extractor(
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

    extractors['Header_attr_MassTable'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=(
            'Units_attr_UnitMass_in_g',
            'Header_attr_HubbleParam'
        ),
        hpath='/Header',
        attr='MassTable',
        convert=lambda vals, raw, path, fname, hpath:
        raw[1] / vals.Header_attr_HubbleParam * vals.Units_attr_UnitMass_in_g,
        units=U.g,
        unit_convert=None
    )

    for Ti in T:
        extractors[
            'RuntimePars_attr_Softening{:s}'.format(softstrings[Ti])
        ] = extractor(
            keytype='header',
            filetype='snapshot',
            dependencies=(
                'Units_attr_UnitLength_in_cm',
                'Header_attr_HubbleParam',
                'Header_attr_Time'
            ),
            hpath='/RuntimePars',
            attr='Softening{:s}'.format(softstrings[Ti]),
            convert=lambda vals, raw, path, fname, hpath:
            raw
            / vals.Header_attr_HubbleParam
            * vals.Units_attr_UnitLength_in_cm,
            units=U.cm,
            unit_convert=U.kpc
        )

    for Ti in T:
        extractors[
            'RuntimePars_attrs_softening{:s}MaxPhys'.format(softstrings[Ti])
        ] = extractor(
            keytype='header',
            filetype='snapshot',
            dependencies=(
                'Units_attr_UnitLength_in_cm',
                'Header_attr_HubbleParam',
                'Header_attr_Time'
            ),
            hpath='/RuntimePars',
            attr='Softening{:s}MaxPhys'.format(softstrings[Ti]),
            convert=lambda vals, raw, path, fname, hpath:
            raw
            / vals.Header_attr_HubbleParam
            * vals.Units_attr_UnitLength_in_cm,
            units=U.cm,
            unit_convert=U.kpc
        )

    extractors['FOF_ContaminationCount'] = extractor(
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

    extractors['FOF_ContaminationMass'] = extractor(
        keytype='fofgroup',
        filetype='group',
        dependencies=(
            'Units_attr_UnitMass_in_g',
            'Header_attr_HubbleParam',
            'Header_attr_Time'
        ),
        hpath='/FOF/ContaminationMass',
        attr=None,
        convert=lambda vals, raw, path, fname, hpath:
        raw * h_a_powers(vals, path, fname, hpath)
        * vals.Units_attr_UnitMass_in_g,
        units=U.g,
        unit_convert=U.Msun
    )

    extractors['FOF_First_SubhaloID'] = extractor(
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

    extractors['FOF_GroupCentreOfPotential'] = extractor(
        keytype='fofgroup',
        filetype='group',
        dependencies=(
            'Units_attr_UnitLength_in_cm',
            'Header_attr_HubbleParam',
            'Header_attr_Time'
        ),
        hpath='/FOF/GroupCentreOfPotential',
        attr=None,
        convert=lambda vals, raw, path, fname, hpath:
        raw * h_a_powers(vals, path, fname, hpath)
        * vals.Units_attr_UnitLength_in_cm,
        units=U.cm,
        unit_convert=U.kpc
    )

    extractors['FOF_GroupLength'] = extractor(
        keytype='fofgroup',
        filetype='group',
        dependencies=tuple(),
        hpath='/FOF/GroupLength',
        attr=None,
        convert=lambda vals, raw, path, fname, hapth:
        raw,
        units=U.dimensionless_unscaled,
        unit_convert=None
    )

    extractors['FOF_GroupMass'] = extractor(
        keytype='fofgroup',
        filetype='group',
        dependencies=(
            'Units_attr_UnitMass_in_g',
            'Header_attr_HubbleParam',
            'Header_attr_Time'
        ),
        hpath='/FOF/GroupMass',
        attr=None,
        convert=lambda vals, raw, path, fname, hpath:
        raw * h_a_powers(vals, path, fname, hpath)
        * vals.Units_attr_UnitMass_in_g,
        units=U.g,
        unit_convert=U.Msun
    )

    extractors['FOF_GroupOffset'] = extractor(
        keytype='fofgroup',
        filetype='group',
        dependencies=tuple(),
        hpath='/FOF/GroupOffset',
        attr=None,
        convert=lambda vals, raw, path, fname, hapth:
        raw,
        units=U.dimensionless_unscaled,
        unit_convert=None
    )

    for overdensity in overdensities:
        extractors['FOF_Group_M_{:s}'.format(overdensity)] = extractor(
            keytype='fofgroup',
            filetype='group',
            dependencies=(
                'Units_attr_UnitMass_in_g',
                'Header_attr_HubbleParam',
                'Header_attr_Time'
            ),
            hpath='/FOF/Group_M_Crit200',
            attr=None,
            convert=lambda vals, raw, path, fname, hpath:
            raw * h_a_powers(vals, path, fname, hpath)
            * vals.Units_attr_UnitMass_in_g,
            units=U.g,
            unit_convert=U.solMass
        )

        extractors['FOF_Group_R_{:s}'.format(overdensity)] = extractor(
            keytype='fofgroup',
            filetype='group',
            dependencies=(
                'Units_attr_UnitLength_in_cm',
                'Header_attr_HubbleParam',
                'Header_attr_Time'
            ),
            hpath='/FOF/Group_R_Crit200',
            attr=None,
            convert=lambda vals, raw, path, fname, hpath:
            raw * h_a_powers(vals, path, fname, hpath)
            * vals.Units_attr_UnitLength_in_cm,
            units=U.cm,
            unit_convert=U.kpc
        )

    extractors['FOF_NumOfSubhalos'] = extractor(
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

    extractors['Subhalo_GroupNumber'] = extractor(
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

    extractors['Subhalo_SubGroupNumber'] = extractor(
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

    extractors['Subhalo_CentreOfPotential'] = extractor(
        keytype='group',
        filetype='group',
        dependencies=(
            'Units_attr_UnitLength_in_cm',
            'Header_attr_HubbleParam',
            'Header_attr_Time'
        ),
        hpath='/Subhalo/CentreOfPotential',
        attr=None,
        convert=lambda vals, raw, path, fname, hpath:
        raw * h_a_powers(vals, path, fname, hpath)
        * vals.Units_attr_UnitLength_in_cm,
        units=U.cm,
        unit_convert=U.kpc
    )

    extractors['Subhalo_Velocity'] = extractor(
        keytype='group',
        filetype='group',
        dependencies=(
            'Units_attr_UnitVelocity_in_cm_per_s',
            'Header_attr_HubbleParam',
            'Header_attr_Time'
        ),
        hpath='/Subhalo/Velocity',
        attr=None,
        convert=lambda vals, raw, path, fname, hpath:
        raw * h_a_powers(vals, path, fname, hpath)
        * vals.Units_attr_UnitVelocity_in_cm_per_s,
        units=U.cm * U.s ** -1,
        unit_convert=U.km * U.s ** -1
    )

    extractors['Subhalo_SubLength'] = extractor(
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

    extractors['Subhalo_SubOffset'] = extractor(
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

    extractors['Subhalo_MassType'] = extractor(
        keytype='group',
        filetype='group',
        dependencies=(
            'Units_attr_UnitMass_in_g',
            'Header_attr_HubbleParam',
            'Header_attr_Time'
        ),
        hpath='/Subhalo/MassType',
        attr=None,
        convert=lambda vals, raw, path, fname, hpath:
        raw * h_a_powers(vals, path, fname, hpath)
        * vals.Units_attr_UnitMass_in_g,
        units=U.g,
        unit_convert=U.solMass
    )

    extractors['FOF_TotNgroups'] = extractor(
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

    extractors['FOF_Group_R_Crit200'] = extractor(
        keytype='fofgroup',
        filetype='group',
        dependencies=(
            'Units_attr_UnitLength_in_cm',
            'Header_attr_HubbleParam',
            'Header_attr_Time'
        ),
        hpath='/FOF/Group_R_Crit200',
        attr=None,
        convert=lambda vals, raw, path, fname, hpath:
        raw * h_a_powers(vals, path, fname, hpath)
        * vals.Units_attr_UnitLength_in_cm,
        units=U.cm,
        unit_convert=U.kpc
    )

    extractors['Subhalo_Vmax'] = extractor(
        keytype='group',
        filetype='group',
        dependencies=(
            'Units_attr_UnitVelocity_in_cm_per_s',
            'Header_attr_HubbleParam',
            'Header_attr_Time'
        ),
        hpath='/Subhalo/Vmax',
        attr=None,
        convert=lambda vals, raw, path, fname, hpath:
        raw * h_a_powers(vals, path, fname, hpath)
        * vals.Units_attr_UnitVelocity_in_cm_per_s,
        units=U.cm * U.s ** -1,
        unit_convert=U.km * U.s ** -1
    )

    extractors['IDs_ParticleID'] = extractor(
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

    for Ti in T:
        extractors['PartType{:d}_ParticleIDs'.format(Ti)] = extractor(
            keytype='particle{:d}'.format(Ti),
            filetype='particle',
            dependencies=tuple(),
            hpath='/PartType{:d}/ParticleIDs'.format(Ti),
            attr=None,
            convert=lambda vals, raw, path, fname, hpath:
            raw,
            units=U.dimensionless_unscaled,
            unit_convert=None
        )

    for Ti in T:
        extractors['PartType{:d}_Coordinates'.format(Ti)] = extractor(
            keytype='particle{:d}'.format(Ti),
            filetype='particle',
            dependencies=(
                'Units_attr_UnitLength_in_cm',
                'Header_attr_HubbleParam',
                'Header_attr_Time'
            ),
            hpath='/PartType{:d}/Coordinates'.format(Ti),
            attr=None,
            convert=lambda vals, raw, path, fname, hpath:
            raw * h_a_powers(vals, path, fname, hpath)
            * vals.Units_attr_UnitLength_in_cm,
            units=U.cm,
            unit_convert=U.kpc
        )

    for Ti in T:
        extractors['PartType{:d}_{:s}'.format(Ti, Vstring)] = extractor(
            keytype='particle{:d}'.format(Ti),
            filetype='particle',
            dependencies=(
                'Units_attr_UnitVelocity_in_cm_per_s',
                'Header_attr_HubbleParam',
                'Header_attr_Time'
            ),
            hpath='/PartType{:d}/{:s}'.format(Ti, Vstring),
            attr=None,
            convert=lambda vals, raw, path, fname, hpath:
            raw * h_a_powers(vals, path, fname, hpath)
            * vals.Units_attr_UnitVelocity_in_cm_per_s,
            units=U.cm * U.s ** -1,
            unit_convert=U.km * U.s ** -1
        )

    for Ti in T:
        extractors['PartType{:d}_GroupNumber'.format(Ti)] = extractor(
            keytype='particle{:d}'.format(Ti),
            filetype='particle',
            dependencies=tuple(),
            hpath='/PartType{:d}/GroupNumber'.format(Ti),
            attr=None,
            convert=lambda vals, raw, path, fname, hpath:
            raw,
            units=U.dimensionless_unscaled,
            unit_convert=None
        )

    for Ti in T:
        extractors['PartType{:d}_SubGroupNumber'.format(Ti)] = extractor(
            keytype='particle{:d}'.format(Ti),
            filetype='particle',
            dependencies=tuple(),
            hpath='/PartType{:d}/SubGroupNumber'.format(Ti),
            attr=None,
            convert=lambda vals, raw, path, fname, hpath:
            raw,
            units=U.dimensionless_unscaled,
            unit_convert=None
        )

    for Ti in (0, 4, 5):
        extractors['PartType{:d}_{:s}'.format(Ti, Mstring)] = extractor(
            keytype='particle{:d}'.format(Ti),
            filetype='particle',
            dependencies=(
                'Units_attr_UnitMass_in_g',
                'Header_attr_HubbleParam',
                'Header_attr_Time'
            ),
            hpath='/PartType{:d}/{:s}'.format(Ti, Mstring),
            attr=None,
            convert=lambda vals, raw, path, fname, hpath:
            raw * h_a_powers(vals, path, fname, hpath)
            * vals.Units_attr_UnitMass_in_g,
            units=U.g,
            unit_convert=U.solMass
        )

    extractors['PartType0_Temperature'] = extractor(
        keytype='particle0',
        filetype='particle',
        dependencies=(
            'Header_attr_HubbleParam',
            'Header_attr_Time'
        ),
        hpath='/PartType0/Temperature',
        attr=None,
        convert=lambda vals, raw, path, fname, hpath:
        raw * h_a_powers(vals, path, fname, hpath),
        units=U.K,
        unit_convert=None
    )

    extractors['PartType0_Density'] = extractor(
        keytype='particle0',
        filetype='particle',
        dependencies=(
            'Units_attr_UnitMass_in_g',
            'Units_attr_UnitLength_in_cm',
            'Header_attr_HubbleParam',
            'Header_attr_Time'
        ),
        hpath='/PartType0/Density',
        attr=None,
        convert=lambda vals, raw, path, fname, hpath:
        raw * h_a_powers(vals, path, fname, hpath)
        * vals.Units_attr_UnitMass_in_g
        * np.power(vals.Units_attr_UnitLength_in_cm, -3),
        units=U.g * U.cm ** -3,
        unit_convert=None
    )

    for Ti in (0, 4):
        for smooth in ('Smoothed', ''):
            for element in elements:
                extractors[
                    'PartType{:d}_{:s}ElementAbundance_{:s}'.format(
                        Ti, smooth, element)
                ] = extractor(
                    keytype='particle{:d}'.format(Ti),
                    filetype='particle',
                    dependencies=(
                        'Header_attr_HubbleParam',
                        'Header_attr_Time'
                    ),
                    hpath='/PartType{:d}/{:s}ElementAbundance/{:s}'.format(
                        Ti, smooth, element),
                    attr=None,
                    convert=lambda vals, raw, path, fname, hpath:
                    raw * h_a_powers(vals, path, fname, hpath),
                    units=U.dimensionless_unscaled,
                    unit_convert=None
                )
            extractors[
                'PartType{:d}__{:s}Metallicity'.format(Ti, smooth)
            ] = extractor(
                keytype='particle{:d}'.format(Ti),
                filetype='particle',
                dependencies=(
                    'Header_attr_HubbleParam',
                    'Header_attr_Time'
                ),
                hpath='/PartType{:d}/{:s}ElementAbundance/{:s}'.format(
                    Ti, smooth, element),
                attr=None,
                convert=lambda vals, raw, path, fname, hpath:
                raw * h_a_powers(vals, path, fname, hpath),
                units=U.dimensionless_unscaled,
                unit_convert=None
            )

    extractors['PartType0_StarFormationRate'] = extractor(
        keytype='particle0',
        filetype='particle',
        dependencies=(
            'Header_attr_HubbleParam',
            'Header_attr_Time'
        ),
        hpath='/PartType0/StarFormationRate',
        attr=None,
        convert=lambda vals, raw, path, fname, hpath:
        raw * h_a_powers(vals, path, fname, hpath),
        units=U.solMass / U.yr,
        unit_convert=None
    )

    for Ti in (0, 4):
        extractors['PartType{:d}_SmoothingLength'.format(Ti)] = extractor(
            keytype='particle{:d}'.format(Ti),
            filetype='particle',
            dependencies=(
                'Units_attr_UnitLength_in_cm',
                'Header_attr_HubbleParam',
                'Header_attr_Time'
            ),
            hpath='/PartType{:d}/SmoothingLength'.format(Ti),
            attr=None,
            convert=lambda vals, raw, path, fname, hpath:
            raw * h_a_powers(vals, path, fname, hpath)
            * vals.Units_attr_UnitLength_in_cm,
            units=U.cm,
            unit_convert=U.kpc
        )

    extractors['PartType4_StellarFormationTime'] = extractor(
        keytype='particle4',
        filetype='particle',
        dependencies=tuple(),
        hpath='/PartType4/StellarFormationTime',
        attr=None,
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.dimensionless_unscaled,
        unit_convert=None
    )

    return extractors

# ~~~HEADER FIELDS~~~
# BoxSize
# E(z)
# ExpansionFactor
# Flag_Cooling
# Flag_DoublePrecision
# Flag_Feedback
# Flag_IC_Info
# Flag_Metals
# Flag_Sfr
# Flag_StellarAge
# H(z)
# HubbleParam
# MassTable
# NTask
# Ngroups
# Nids
# Nsubgroups
# NumFilesPerSnapshot
# NumPart_ThisFile
# NumPart_Total
# NumPart_Total_HighWord
# Omega0
# OmegaBaryon
# OmegaLambda
# Redshift
# RunLabel
# SendOffSetTask
# Time
# TotNgroups
# TotNids
# TotNsubgroups

# ~~~UNITS FIELDS~~~
# UnitDensity_in_cgs
# UnitEnergy_in_cgs
# UnitLength_in_cm
# UnitMass_in_g
# UnitPressure_in_cgs
# UnitTime_in_s
# UnitVelocity_in_cm_per_s
