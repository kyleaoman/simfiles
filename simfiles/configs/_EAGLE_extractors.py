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
        Vstring='Velocity',
        EOSstring='OnEquationOfState'
):

    extractors = dict()

    extractors['Constants_attr_AVOGADRO'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/Constants',
        attr='AVOGADRO',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.dimensionless_unscaled,
        unit_convert=None
    )

    extractors['Constants_attr_BOLTZMANN'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/Constants',
        attr='BOLTZMANN',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.erg * U.K ** -1,
        unit_convert=None
    )

    extractors['Constants_attr_C'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/Constants',
        attr='C',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.cm * U.s ** -1,
        unit_convert=U.km * U.s ** -1
    )

    extractors['Constants_attr_CM_PER_MPC'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/Constants',
        attr='CM_PER_MPC',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.cm * U.Mpc ** -1,
        unit_convert=None
    )

    extractors['Constants_attr_ELECTRONCHARGE'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/Constants',
        attr='ELECTRONCHARGE',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.statC,
        unit_convert=None
    )

    extractors['Constants_attr_ELECTRONMASS'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/Constants',
        attr='ELECTRONMASS',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.g,
        unit_convert=None
    )

    extractors['Constants_attr_EV_TO_ERG'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/Constants',
        attr='EV_TO_ERG',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.eV * U.erg ** -1,
        unit_convert=None
    )

    extractors['Constants_attr_GAMMA'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/Constants',
        attr='GAMMA',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.dimensionless_unscaled,
        unit_convert=None
    )

    extractors['Constants_attr_GAS_CONST'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/Constants',
        attr='GAS_CONST',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.erg * U.K ** -1 * U.mol ** -1,
        unit_convert=None
    )

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

    extractors['Constants_attr_HUBBLE'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/Constants',
        attr='HUBBLE',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.s ** -1,
        unit_convert=None
    )

    extractors['Constants_attr_PI'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/Constants',
        attr='PI',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.dimensionless_unscaled,
        unit_convert=None
    )

    extractors['Constants_attr_PLANCK'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/Constants',
        attr='PLANCK',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.erg * U.s,
        unit_convert=None
    )

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

    extractors['Constants_attr_RAD_CONST'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/Constants',
        attr='RAD_CONST',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.erg * U.cm ** -3 * U.K ** -4,
        unit_convert=None
    )

    extractors['Constants_attr_SEC_PER_MEGAYEAR'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/Constants',
        attr='SEC_PER_MEGAYEAR',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.s * U.Myr ** -1,
        unit_convert=None
    )

    extractors['Constants_attr_SEC_PER_YEAR'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/Constants',
        attr='SEC_PER_YEAR',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.s * U.yr ** -1,
        unit_convert=None
    )

    extractors['Constants_attr_SOLAR_LUM'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/Constants',
        attr='SOLAR_LUM',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.erg * U.s ** -1,
        unit_convert=None
    )

    extractors['Constants_attr_SOLAR_MASS'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/Constants',
        attr='SOLAR_MASS',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.g,
        unit_convert=None
    )

    extractors['Constants_attr_STEFAN'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/Constants',
        attr='STEFAN',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.erg * U.cm ** -3 * U.K ** -4,
        unit_convert=None
    )

    extractors['Constants_attr_THOMPSON'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/Constants',
        attr='THOMPSON',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.cm ** 2,
        unit_convert=None
    )

    extractors['Constants_attr_T_CMB0'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/Constants',
        attr='T_CMB0',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.K,
        unit_convert=None
    )

    extractors['Constants_attr_Z_Solar'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/Constants',
        attr='Z_Solar',
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

    extractors['Header_attr_E(z)'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/Header',
        attr='E(z)',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.dimensionless_unscaled,
        unit_convert=None
    )

    extractors['Header_attr_ExpansionFactor'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/Header',
        attr='ExpansionFactor',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.dimensionless_unscaled,
        unit_convert=None
    )

    extractors['Header_attr_Flag_Cooling'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/Header',
        attr='Flag_Cooling',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=None,
        unit_convert=None
    )

    extractors['Header_attr_Flag_DoublePrecision'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/Header',
        attr='Flag_DoublePrecision',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=None,
        unit_convert=None
    )

    extractors['Header_attr_Flag_Feedback'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/Header',
        attr='Flag_Feedback',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=None,
        unit_convert=None
    )

    extractors['Header_attr_Flag_IC_Info'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/Header',
        attr='Flag_IC_Info',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=None,
        unit_convert=None
    )

    extractors['Header_attr_Flag_Metals'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/Header',
        attr='Flag_Metals',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=None,
        unit_convert=None
    )

    extractors['Header_attr_Flag_Sfr'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/Header',
        attr='Flag_Sfr',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=None,
        unit_convert=None
    )

    extractors['Header_attr_Flag_StellarAge'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/Header',
        attr='Flag_StellarAge',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=None,
        unit_convert=None
    )

    extractors['Header_attr_H(z)'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/Header',
        attr='H(z)',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.s ** -1,
        unit_convert=U.km * U.s ** -1 * U.Mpc ** -1
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

    extractors['Header_attr_NumPart_Total'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/Header',
        attr='NumPart_Total',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=None,
        unit_convert=None
    )

    extractors['Header_attr_NumPart_Total_HighWord'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/Header',
        attr='NumPart_Total_HighWord',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=None,
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

    extractors['Header_attr_OmegaLambda'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/Header',
        attr='OmegaLambda',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.dimensionless_unscaled,
        unit_convert=None
    )

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

    extractors['Header_attr_TotNgroups'] = extractor(
        keytype='header',
        filetype='group',
        dependencies=tuple(),
        hpath='/Header',
        attr='TotNgroups',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=None,
        unit_convert=None
    )

    extractors['Header_attr_TotNsubgroups'] = extractor(
        keytype='header',
        filetype='group',
        dependencies=tuple(),
        hpath='/Header',
        attr='TotNsubgroups',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=None,
        unit_convert=None
    )

    extractors['FOF_attr_TotNgroups'] = extractor(
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

    extractors['Subhalo_attr_TotNgroups'] = extractor(
        keytype='header',
        filetype='group',
        dependencies=tuple(),
        hpath='/Subhalo',
        attr='TotNgroups',
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

    for Ti in set((0, 2, 3, 4, 5)).intersection(T):
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

    for Ti in set((0, 4)).intersection(T):
        extractors['PartType{:d}_AExpMaximumTemperature'.format(Ti)] = \
            extractor(
                keytype='particle{:d}'.format(Ti),
                filetype='particle',
                dependencies=(
                    'Header_attr_HubbleParam',
                    'Header_attr_Time'
                ),
                hpath='/PartType{:d}/AExpMaximumTemperature'.format(Ti),
                attr=None,
                convert=lambda vals, raw, path, fname, hpath:
                raw * h_a_powers(vals, path, fname, hpath),
                units=U.K,
                unit_convert=None
            )

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
                'PartType{:d}_{:s}Metallicity'.format(Ti, smooth)
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

        extractors['PartType{:d}_IronMassFracFromSNIa'.format(Ti)] = extractor(
            keytype='particle{:d}'.format(Ti),
            filetype='particle',
            dependencies=(
                'Header_attr_HubbleParam',
                'Header_attr_Time'
            ),
            hpath='/PartType{:d}/IronMassFracFromSNIa'.format(Ti),
            attr=None,
            convert=lambda vals, raw, path, fname, hpath:
            raw * h_a_powers(vals, path, fname, hpath),
            units=U.dimensionless_unscaled,
            unit_convert=None
        )

        extractors['PartType{:d}_MaximumTemperature'.format(Ti)] = extractor(
            keytype='particle{:d}'.format(Ti),
            filetype='particle',
            dependencies=(
                'Header_attr_HubbleParam',
                'Header_attr_Time'
            ),
            hpath='/PartType{:d}/MaximumTemperature'.format(Ti),
            attr=None,
            convert=lambda vals, raw, path, fname, hpath:
            raw * h_a_powers(vals, path, fname, hpath),
            units=U.K,
            unit_convert=None
        )

        for src in ('AGB', 'SNII', 'SNIa'):
            extractors['PartType{:d}_MetalMassFracFrom{:s}'.format(
                Ti, src)] = extractor(
                    keytype='particle{:d}'.format(Ti),
                    filetype='particle',
                    dependencies=(
                        'Header_attr_HubbleParam',
                        'Header_attr_Time'
                    ),
                    hpath='/PartType{:d}/MetalMassFracFrom{:s}'.format(
                        Ti, src),
                    attr=None,
                    convert=lambda vals, raw, path, fname, hpath:
                    raw * h_a_powers(vals, path, fname, hpath),
                    units=U.dimensionless_unscaled,
                    unit_convert=None
                )

            extractors['PartType{:d}_TotalMassFrom{:s}'.format(Ti, src)] = \
                extractor(
                    keytype='particle{:d}'.format(Ti),
                    filetype='particle',
                    dependencies=(
                        'Header_attr_HubbleParam',
                        'Header_attr_Time'
                    ),
                    hpath='/PartType{:d}/TotalMassFrom{:s}'.format(Ti, src),
                    attr=None,
                    convert=lambda vals, raw, path, fname, hpath:
                    raw * h_a_powers(vals, path, fname, hpath),
                    units=U.dimensionless_unscaled,
                    unit_convert=None
                )

        extractors['PartType{:d}_MetalMassWeightedRedshift'.format(Ti)] = \
            extractor(
                keytype='particle{:d}'.format(Ti),
                filetype='particle',
                dependencies=(
                    'Header_attr_HubbleParam',
                    'Header_attr_Time'
                ),
                hpath='/PartType{:d}/MetalMassWeightedRedshift'.format(Ti),
                attr=None,
                convert=lambda vals, raw, path, fname, hpath:
                raw * h_a_powers(vals, path, fname, hpath),
                units=U.dimensionless_unscaled,
                unit_convert=None
            )

        extractors['PartType{:d}_SmoothedIronMassFracFromSNIa'.format(Ti)] = \
            extractor(
                keytype='particle{:d}'.format(Ti),
                filetype='particle',
                dependencies=(
                    'Header_attr_HubbleParam',
                    'Header_attr_Time'
                ),
                hpath='/PartType{:d}/SmoothedIronMassFracFromSNIa'.format(Ti),
                attr=None,
                convert=lambda vals, raw, path, fname, hpath:
                raw * h_a_powers(vals, path, fname, hpath),
                units=U.dimensionless_unscaled,
                unit_convert=None
            )

    for Ti in set((0, 4, 5)).intersection(T):
        extractors['PartType{:d}_HostHalo_TVir_Mass'.format(Ti)] = extractor(
            keytype='particle{:d}'.format(Ti),
            filetype='particle',
            dependencies=(
                'Header_attr_HubbleParam',
                'Header_attr_Time',
            ),
            hpath='/PartType{:d}/HostHalo_TVir_Mass'.format(Ti),
            attr=None,
            convert=lambda vals, raw, path, fname, hpath:
            raw * h_a_powers(vals, path, fname, hpath),
            units=U.K,
            unit_convert=None
        )

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

    extractors['PartType0_Entropy'] = extractor(
        keytype='particle0',
        filetype='particle',
        dependencies=(
            'Header_attr_HubbleParam',
            'Header_attr_Time',
            'Constants_attr_GAMMA',
            'Units_attr_UnitPressure_in_cgs',
            'Units_attr_UnitDensity_in_cgs'
        ),
        hpath='/PartType0/Entropy',
        attr=None,
        convert=lambda vals, raw, path, fname, hpath:
        raw
        * np.power(
            vals.Header_attr_HubbleParam,
            2 - 2 * vals.Constants_attr_GAMMA
        )
        * vals.Units_attr_UnitPressure_in_cgs
        * np.power(
            vals.Units_attr_UnitDensity_in_cgs,
            -vals.Constants_attr_GAMMA
        ),
        units=U.erg * U.K ** -1,
        unit_convert=None
    )

    extractors['PartType0_InternalEnergy'] = extractor(
        keytype='particle0',
        filetype='particle',
        dependencies=(
            'Header_attr_HubbleParam',
            'Header_attr_Time'
        ),
        hpath='PartType0/InternalEnergy',
        attr=None,
        convert=lambda vals, raw, path, fname, hpath:
        raw * h_a_powers(vals, path, fname, hpath)
        * np.power(vals.Units_attr_UnitVelocity_in_cm_per_s, 2),
        units=U.cm ** 2 * U.s ** -2,
        unit_convert=U.km ** 2 * U.s ** -2
    )

    extractors['PartType0_{:s}'.format(EOSstring)] = extractor(
        keytype='particle0',
        filetype='particle',
        dependencies=tuple(),
        hpath='/PartType0/{:s}'.format(EOSstring),
        attr=None,
        convert=lambda vals, raw, path, fname, hpath:
        raw,
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

    extractors['PartType4_BirthDensity'] = extractor(
        keytype='particle4',
        filetype='particle',
        dependencies=(
            'Header_attr_HubbleParam',
            'Header_attr_Time',
            'Units_attr_UnitDensity_in_cgs'
        ),
        hpath='/PartType4/BirthDensity',
        attr=None,
        convert=lambda vals, raw, path, fname, hpath:
        raw * h_a_powers(vals, path, fname, hpath)
        * vals.Units_attr_UnitDensity_in_cgs,
        units=U.g * U.cm ** -3,
        unit_convert=None
    )

    extractors['PartType4_Feedback_EnergyFraction'] = extractor(
        keytype='particle4',
        filetype='particle',
        dependencies=(
            'Header_attr_HubbleParam',
            'Header_attr_Time'
        ),
        hpath='/PartType4/Feedback_EnergyFraction',
        attr=None,
        convert=lambda vals, raw, path, fname, hpath:
        raw * h_a_powers(vals, path, fname, hpath),
        units=U.dimensionless_unscaled,
        unit_convert=None
    )

    extractors['PartType4_HostHalo_TVir'] = extractor(
        keytype='particle4',
        filetype='particle',
        dependencies=(
            'Header_attr_HubbleParam',
            'Header_attr_Time'
        ),
        hpath='/PartType4/HostHalo_TVir',
        attr=None,
        convert=lambda vals, raw, path, fname, hpath:
        raw * h_a_powers(vals, path, fname, hpath),
        units=U.K,
        unit_convert=None
    )

    extractors['PartType4_InitialMass'] = extractor(
        keytype='particle4',
        filetype='particle',
        dependencies=(
            'Header_attr_HubbleParam',
            'Header_attr_Time',
            'Units_attr_UnitMass_in_g'
        ),
        hpath='/PartType4/InitialMass',
        attr=None,
        convert=lambda vals, raw, path, fname, hpath:
        raw * h_a_powers(vals, path, fname, hpath)
        * vals.Units_attr_UnitMass_in_g,
        units=U.Msun,
        unit_convert=None
    )

    extractors['PartType4_PreviousStellarEnrichment'] = extractor(
        keytype='particle4',
        filetype='particle',
        dependencies=(
            'Header_attr_HubbleParam',
            'Header_attr_Time'
        ),
        hpath='/PartType4/PreviousStellarEnrichment',
        attr=None,
        convert=lambda vals, raw, path, fname, hpath:
        raw * h_a_powers(vals, path, fname, hpath),
        units=U.dimensionless_unscaled,
        unit_convert=None
    )

    extractors['PartType4_StellarEnrichmentCounter'] = extractor(
        keytype='particle4',
        filetype='particle',
        dependencies=(
            'Header_attr_HubbleParam',
            'Header_attr_Time'
        ),
        hpath='/PartType4/StellarEnrichmentCounter',
        attr=None,
        convert=lambda vals, raw, path, fname, hpath:
        raw * h_a_powers(vals, path, fname, hpath),
        units=U.dimensionless_unscaled,
        unit_convert=None
    )

    extractors['PartType4_StellarFormationTime'] = extractor(
        keytype='particle4',
        filetype='particle',
        dependencies=(
            'Header_attr_HubbleParam',
            'Header_attr_Time'
        ),
        hpath='/PartType4/StellarFormationTime',
        attr=None,
        convert=lambda vals, raw, path, fname, hpath:
        raw * h_a_powers(vals, path, fname, hpath),
        units=U.dimensionless_unscaled,
        unit_convert=None
    )

    extractors['PartType5_BH_AccretionLength'] = extractor(
        keytype='particle5',
        filetype='particle',
        dependencies=(
            'Header_attr_HubbleParam',
            'Header_attr_Time',
            'Units_attr_UnitLength_in_cm'
        ),
        hpath='/PartType5/BH_AccretionLength',
        attr=None,
        convert=lambda vals, raw, path, fname, hpath:
        raw * h_a_powers(vals, path, fname, hpath)
        * vals.Units_attr_UnitLength_in_cm,
        units=U.cm,
        unit_convert=U.kpc
    )

    extractors['PartType5_BH_CumlAccrMass'] = extractor(
        keytype='particle5',
        filetype='particle',
        dependencies=(
            'Header_attr_HubbleParam',
            'Header_attr_Time',
            'Units_attr_UnitMass_in_g'
        ),
        hpath='/PartType5/BH_CumlAccrMass',
        attr=None,
        convert=lambda vals, raw, path, fname, hpath:
        raw * h_a_powers(vals, path, fname, hpath)
        * vals.Units_attr_UnitMass_in_g,
        units=U.g,
        unit_convert=U.Msun
    )

    extractors['PartType5_BH_CumlNumSeeds'] = extractor(
        keytype='particle5',
        filetype='particle',
        dependencies=(
            'Header_attr_HubbleParam',
            'Header_attr_Time'
        ),
        hpath='/PartType5/BH_CumlNumSeeds',
        attr=None,
        convert=lambda vals, raw, path, fname, hpath:
        raw * h_a_powers(vals, path, fname, hpath),
        units=U.dimensionless_unscaled,
        unit_convert=None
    )

    extractors['PartType5_BH_Density'] = extractor(
        keytype='particle5',
        filetype='particle',
        dependencies=(
            'Header_attr_HubbleParam',
            'Header_attr_Time'
        ),
        hpath='/PartType5/BH_Density',
        attr=None,
        convert=lambda vals, raw, path, fname, hpath:
        raw * h_a_powers(vals, path, fname, hpath)
        * vals.Units_attr_UnitMass_in_g
        * vals.Units_attr_UnitLength_in_cm ** -3,
        units=U.g * U.cm ** -3,
        unit_convert=None
    )

    extractors['PartType5_BH_EnergyReservoir'] = extractor(
        keytype='particle5',
        filetype='particle',
        dependencies=(
            'Header_attr_HubbleParam',
            'Header_attr_Time',
            'Units_attr_UnitMass_in_g',
            'Units_attr_UnitLength_in_cm',
            'Units_attr_UnitTime_in_s'
        ),
        hpath='/PartType5/BH_EnergyReservoir',
        attr=None,
        convert=lambda vals, raw, path, fname, hpath:
        raw * h_a_powers(vals, path, fname, hpath)
        * vals.Units_attr_UnitMass_in_g
        * vals.Units_attr_UnitLength_in_cm ** 2
        * vals.Units_attr_UnitTime_in_s ** -2,
        units=U.erg,
        unit_convert=None
    )

    extractors['PartType5_BH_FormationTime'] = extractor(
        keytype='particle5',
        filetype='particle',
        dependencies=(
            'Header_attr_HubbleParam',
            'Header_attr_Time'
        ),
        hpath='/PartType5/BH_FormationTime',
        attr=None,
        convert=lambda vals, raw, path, fname, hpath:
        raw * h_a_powers(vals, path, fname, hpath),
        units=U.dimensionless_unscaled,
        unit_convert=None
    )

    extractors['PartType5_BH_Mass'] = extractor(
        keytype='particle5',
        filetype='particle',
        dependencies=(
            'Header_attr_HubbleParam',
            'Header_attr_Time'
        ),
        hpath='/PartType5/BH_Mass',
        attr=None,
        convert=lambda vals, raw, path, fname, hpath:
        raw * h_a_powers(vals, path, fname, hpath)
        * vals.Units_attr_UnitMass_in_g,
        units=U.g,
        unit_convert=U.Msun
    )

    extractors['PartType5_BH_Mdot'] = extractor(
        keytype='particle5',
        filetype='particle',
        dependencies=(
            'Header_attr_HubbleParam',
            'Header_attr_Time'
        ),
        hpath='/PartType5/BH_Mdot',
        attr=None,
        convert=lambda vals, raw, path, fname, hpath:
        raw * h_a_powers(vals, path, fname, hpath)
        * vals.Units_attr_UnitMass_in_g
        / vals.Units_attr_UnitTime_in_s,
        units=U.g / U.s,
        unit_convert=U.Msun / U.yr
    )

    extractors['PartType5_BH_MostMassiveProgenitorID'] = extractor(
        keytype='particle5',
        filetype='particle',
        dependencies=tuple(),
        hpath='/PartType5/BH_MostMassiveProgenitorID',
        attr=None,
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.dimensionless_unscaled,
        unit_convert=None
    )

    extractors['PartType5_BH_Pressure'] = extractor(
        keytype='particle5',
        filetype='particle',
        dependencies=(
            'Header_attr_HubbleParam',
            'Header_attr_Time',
            'Units_attr_UnitMass_in_g',
            'Units_attr_UnitVelocity_in_cm_per_s',
            'Units_attr_UnitLength_in_cm'
        ),
        hpath='/PartType5/BH_Pressure',
        attr=None,
        convert=lambda vals, raw, path, fname, hpath:
        raw * np.power(vals.Header_attr_HubbleParam, 2)
        * np.power(vals.Header_attr_Time, -3 * vals.Constants_attr_GAMMA)
        * vals.Units_attr_UnitMass_in_g
        * vals.Units_attr_UnitVelocity_in_cm_per_s ** 2
        * vals.Units_attr_UnitLength_in_cm ** -3,
        units=U.g * U.cm ** -1 * U.s ** -2,
        unit_convert=None
    )

    extractors['PartType5_BH_SoundSpeed'] = extractor(
        keytype='particle5',
        filetype='particle',
        dependencies=(
            'Header_attr_HubbleParam',
            'Header_attr_Time',
            'Units_attr_UnitVelocity_in_cm_per_s'
        ),
        hpath='/PartType5/BH_SoundSpeed',
        attr=None,
        convert=lambda vals, raw, path, fname, hpath:
        raw * h_a_powers(vals, path, fname, hpath)
        * vals.Units_attr_UnitVelocity_in_cm_per_s,
        units=U.cm * U.s ** -1,
        unit_convert=U.km * U.s ** -1
    )

    extractors['PartType5_BH_SurroundingGasVel'] = extractor(
        keytype='particle5',
        filetype='particle',
        dependencies=(
            'Header_attr_HubbleParam',
            'Header_attr_Time',
            'Units_attr_UnitVelocity_in_cm_per_s'
        ),
        hpath='/PartType5/BH_SurroundingGasVel',
        attr=None,
        convert=lambda vals, raw, path, fname, hpath:
        raw * h_a_powers(vals, path, fname, hpath)
        * vals.Units_attr_UnitVelocity_in_cm_per_s,
        units=U.cm * U.s ** -1,
        unit_convert=U.km * U.s ** -1
    )

    extractors['PartType5_BH_TimeLastMerger'] = extractor(
        keytype='particle5',
        filetype='particle',
        dependencies=(
            'Header_attr_HubbleParam',
            'Header_attr_Time'
        ),
        hpath='/PartType5/BH_TimeLastMerger',
        attr=None,
        convert=lambda vals, raw, path, fname, hpath:
        raw * h_a_powers(vals, path, fname, hpath),
        units=U.dimensionless_unscaled,
        unit_convert=None
    )

    extractors['PartType5_BH_WeightedDensity'] = extractor(
        keytype='particle5',
        filetype='particle',
        dependencies=(
            'Header_attr_HubbleParam',
            'Header_attr_Time',
            'Units_attr_UnitMass_in_g',
            'Units_attr_UnitLength_in_cm'
        ),
        hpath='/PartType5/BH_WeightedDensity',
        attr=None,
        convert=lambda vals, raw, path, fname, hpath:
        raw * h_a_powers(vals, path, fname, hpath)
        * vals.Units_attr_UnitMass_in_g
        * vals.Units_attr_UnitLength_in_cm ** -3,
        units=U.g * U.cm ** -3,
        unit_convert=None
    )

    extractors['RuntimePars_attr_AGB_EnergyTransferOn'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='AGB_EnergyTransferOn',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=None,
        unit_convert=None
    )

    extractors['RuntimePars_attr_AGB_MassTransferOn'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='AGB_MassTransferOn',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=None,
        unit_convert=None
    )

    extractors['RuntimePars_attr_ArtBulkViscConst'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='ArtBulkViscConst',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.dimensionless_unscaled,
        unit_convert=None
    )

    extractors['RuntimePars_attr_ArtBulkViscConstMin'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='ArtBulkViscConstMin',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.dimensionless_unscaled,
        unit_convert=None
    )

    extractors['RuntimePars_attr_ArtDiffConst'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='ArtDiffConst',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.dimensionless_unscaled,
        unit_convert=None
    )

    extractors['RuntimePars_attr_ArtDiffConstMin'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='ArtDiffConstMin',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.dimensionless_unscaled,
        unit_convert=None
    )

    extractors['RuntimePars_attr_BH_ConstantHeatTemp'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='BH_ConstantHeatTemp',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.K,
        unit_convert=None
    )

    extractors['RuntimePars_attr_BH_MaxHeatLimit'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='BH_MaxHeatLimit',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.dimensionless_unscaled,
        unit_convert=None
    )

    extractors['RuntimePars_attr_BH_MaxHeatTemp'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='BH_MaxHeatTemp',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.K,
        unit_convert=None
    )

    extractors['RuntimePars_attr_BH_MaxMergingDistanceFactor'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='BH_MaxMergingDistanceFactor',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.dimensionless_unscaled,
        unit_convert=None
    )

    extractors['RuntimePars_attr_BH_MaxRepositionDistanceFactor'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='BH_MaxRepositionDistanceFactor',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.dimensionless_unscaled,
        unit_convert=None
    )

    extractors['RuntimePars_attr_BH_MinHeatLimit'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='BH_MinHeatLimit',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.dimensionless_unscaled,
        unit_convert=None
    )

    extractors['RuntimePars_attr_BH_MinHeatTemp'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='BH_MinHeatTemp',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.K,
        unit_convert=None
    )

    extractors['RuntimePars_attr_BH_feedback_mode'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='BH_feedback_mode',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=None,
        unit_convert=None
    )

    extractors['RuntimePars_attr_BH_maxHeatingProbability'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='BH_maxHeatingProbability',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.dimensionless_unscaled,
        unit_convert=None
    )

    extractors['RuntimePars_attr_BlackHoleAccretionFactor'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='BlackHoleAccretionFactor',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.dimensionless_unscaled,
        unit_convert=None
    )

    extractors['RuntimePars_attr_BlackHoleAccretionSlope'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='BlackHoleAccretionSlope',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.dimensionless_unscaled,
        unit_convert=None
    )

    extractors['RuntimePars_attr_BlackHoleEddingtonFactor'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='BlackHoleEddingtonFactor',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.dimensionless_unscaled,
        unit_convert=None
    )

    extractors['RuntimePars_attr_BlackHoleFeedbackFactor'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='BlackHoleFeedbackFactor',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.dimensionless_unscaled,
        unit_convert=None
    )

    extractors['RuntimePars_attr_BlackHoleMaxAccretionRadius'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='BlackHoleMaxAccretionRadius',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.dimensionless_unscaled,
        unit_convert=None
    )

    extractors['RuntimePars_attr_BlackHoleNgbFactor'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='BlackHoleNgbFactor',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.dimensionless_unscaled,
        unit_convert=None
    )

    extractors['RuntimePars_attr_BlackHoleNumberOfNeighboursToHeat'] = \
        extractor(
            keytype='header',
            filetype='snapshot',
            dependencies=tuple(),
            hpath='/RuntimePars',
            attr='BlackHoleNumberOfNeighboursToHeat',
            convert=lambda vals, raw, path, fname, hpath:
            raw,
            units=U.dimensionless_unscaled,
            unit_convert=None
        )

    extractors['RuntimePars_attr_BlackHoleRadiativeEfficiency'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='BlackHoleRadiativeEfficiency',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.dimensionless_unscaled,
        unit_convert=None
    )

    extractors['RuntimePars_attr_BlackHoleViscousAlpha'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='BlackHoleViscousAlpha',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.dimensionless_unscaled,
        unit_convert=None
    )

    extractors['RuntimePars_attr_BoxSize'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=(
            'Units_attr_UnitLength_in_cm',
            'Header_attr_HubbleParam'
        ),
        hpath='/RuntimePars',
        attr='BoxSize',
        convert=lambda vals, raw, path, fname, hpath:
        raw / vals.Header_attr_HubbleParam * vals.Units_attr_UnitLength_in_cm,
        units=U.cm,
        unit_convert=U.kpc
    )

    extractors['RuntimePars_attr_CalciumOverSilicon'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='CalciumOverSilicon',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.dimensionless_unscaled,
        unit_convert=None
    )

    extractors['RuntimePars_attr_CoolingOn'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='CoolingOn',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=None,
        unit_convert=None
    )

    extractors['RuntimePars_attr_CourantFac'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='CourantFac',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.dimensionless_unscaled,
        unit_convert=None
    )

    extractors['RuntimePars_attr_DesLinkNgb'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='DesLinkNgb',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.dimensionless_unscaled,
        unit_convert=None
    )

    extractors['RuntimePars_attr_DesNumNgb'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='DesNumNgb',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.dimensionless_unscaled,
        unit_convert=None
    )

    extractors['RuntimePars_attr_DesNumNgbStar'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='DesNumNgbStar',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.dimensionless_unscaled,
        unit_convert=None
    )

    extractors['RuntimePars_attr_DesNumNgbYoungStar'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='DesNumNgbYoungStar',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.dimensionless_unscaled,
        unit_convert=None
    )

    extractors['RuntimePars_attr_EOS_Cool_GammaEffective'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='EOS_Cool_GammaEffective',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.dimensionless_unscaled,
        unit_convert=None
    )

    extractors['RuntimePars_attr_EOS_Cool_MinOverDens'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='EOS_Cool_MinOverDens',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.dimensionless_unscaled,
        unit_convert=None
    )

    extractors['RuntimePars_attr_EOS_Cool_MinPhysDens_HpCM3'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=(
            'Constants_attr_PROTONMASS',
        ),
        hpath='/RuntimePars',
        attr='EOS_Cool_MinPhysDens_HpCM3',
        convert=lambda vals, raw, path, fname, hpath:
        vals.Constants_attr_PROTONMASS.to(U.g).value * raw,
        units=U.g * U.cm ** -3,
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

    extractors['RuntimePars_attr_EOS_Jeans_MinOverDens'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='EOS_Jeans_MinOverDens',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.dimensionless_unscaled,
        unit_convert=None
    )

    extractors['RuntimePars_attr_EOS_Jeans_MinPhysDens_HpCM3'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=(
            'Constants_attr_PROTONMASS',
        ),
        hpath='/RuntimePars',
        attr='EOS_Jeans_MinPhysDens_HpCM3',
        convert=lambda vals, raw, path, fname, hpath:
        vals.Constants_attr_PROTONMASS.to(U.g).value * raw,
        units=U.g * U.cm ** -3,
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

    extractors['RuntimePars_attr_EOS_NormPhysDens_HpCM3'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=(
            'Constants_attr_PROTONMASS',
        ),
        hpath='/RuntimePars',
        attr='EOS_NormPhysDens_HpCM3',
        convert=lambda vals, raw, path, fname, hpath:
        vals.Constants_attr_PROTONMASS.to(U.g).value * raw,
        units=U.g * U.cm ** -3,
        unit_convert=None
    )

    extractors['RuntimePars_attr_ErrTolForceAcc'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='ErrTolForceAcc',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.dimensionless_unscaled,
        unit_convert=None
    )

    extractors['RuntimePars_attr_ErrTolIntAccuracy'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='ErrTolIntAccuracy',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.dimensionless_unscaled,
        unit_convert=None
    )

    extractors['RuntimePars_attr_ErrTolTheta'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='ErrTolTheta',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.dimensionless_unscaled,
        unit_convert=None
    )

    extractors['RuntimePars_attr_ErrTolThetaSubfind'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='ErrTolThetaSubfind',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.dimensionless_unscaled,
        unit_convert=None
    )

    extractors['RuntimePars_attr_Generations'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='Generations',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.dimensionless_unscaled,
        unit_convert=None
    )

    extractors['RuntimePars_attr_HubbleParam'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='HubbleParam',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.dimensionless_unscaled,
        unit_convert=None
    )

    extractors['RuntimePars_attr_IMF_Exponent'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='IMF_Exponent',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.dimensionless_unscaled,
        unit_convert=None
    )

    extractors['RuntimePars_attr_IMF_LifetimeModel'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='IMF_LifetimeModel',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=None,
        unit_convert=None
    )

    extractors['RuntimePars_attr_IMF_MaxMass_MSUN'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='IMF_MaxMass_MSUN',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.Msun,
        unit_convert=None
    )

    extractors['RuntimePars_attr_IMF_MinMass_MSUN'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='IMF_MinMass_MSUN',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.Msun,
        unit_convert=None
    )

    extractors['RuntimePars_attr_IMF_Model'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='IMF_Model',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=None,
        unit_convert=None
    )

    for element in elements:
        extractors['RuntimePars_attr_InitAbundance_{:s}'.format(element)] = \
            extractor(
                keytype='header',
                filetype='snapshot',
                dependencies=tuple(),
                hpath='/RuntimePars',
                attr='InitAbundance_{:s}'.format(element),
                convert=lambda vals, raw, path, fname, hpath:
                raw,
                units=U.dimensionless_unscaled,
                unit_convert=None
            )

    extractors['RuntimePars_attr_InitGasTemp'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='InitGasTemp',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.K,
        unit_convert=None
    )

    extractors['RuntimePars_attr_InitMetallicity'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='InitMetallicity',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.dimensionless_unscaled,
        unit_convert=None
    )

    extractors['RuntimePars_attr_MaxNumNgbDeviation'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='MaxNumNgbDeviation',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.dimensionless_unscaled,
        unit_convert=None
    )

    extractors['RuntimePars_attr_MaxRMSDisplacementFac'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='MaxRMSDisplacementFac',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.dimensionless_unscaled,
        unit_convert=None
    )

    extractors['RuntimePars_attr_MaxSizeTimestep'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='MaxSizeTimestep',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.dimensionless_unscaled,
        unit_convert=None
    )

    extractors['RuntimePars_attr_MaxSmoothingLengthChange'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='MaxSmoothingLengthChange',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.dimensionless_unscaled,
        unit_convert=None
    )

    extractors['RuntimePars_attr_MetDepCoolingOn'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='MetDepCoolingOn',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=None,
        unit_convert=None
    )

    extractors['RuntimePars_attr_MinFoFMassForNewSeed_Msun'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='MinFoFMassForNewSeed_Msun',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.Msun,
        unit_convert=None
    )

    extractors['RuntimePars_attr_MinGasHsmlFractional'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='MinGasHsmlFractional',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.dimensionless_unscaled,
        unit_convert=None
    )

    extractors['RuntimePars_attr_MinGasTemp'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='MinGasTemp',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.K,
        unit_convert=None
    )

    extractors['RuntimePars_attr_MinSizeTimestep'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='MinSizeTimestep',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.dimensionless_unscaled,
        unit_convert=None
    )

    extractors['RuntimePars_attr_Omega0'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='Omega0',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.dimensionless_unscaled,
        unit_convert=None
    )

    extractors['RuntimePars_attr_OmegaBaryon'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='OmegaBaryon',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.dimensionless_unscaled,
        unit_convert=None
    )

    extractors['RuntimePars_attr_OmegaLambda'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='OmegaLambda',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.dimensionless_unscaled,
        unit_convert=None
    )

    extractors['RuntimePars_attr_REION_H_Heating_EVpH'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=(
            'Constants_attr_PROTONMASS',
        ),
        hpath='/RuntimePars',
        attr='REION_H_Heating_EVpH',
        convert=lambda vals, raw, path, fname, hpath:
        raw / vals.Constants_attr_PROTONMASS.to(U.g).value,
        units=U.eV * U.g ** -1,
        unit_convert=U.erg * U.g ** -1
    )

    extractors['RuntimePars_attr_REION_H_ZCenter'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='REION_H_ZCenter',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.dimensionless_unscaled,
        unit_convert=None
    )

    extractors['RuntimePars_attr_REION_He_Heating_EVpH'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=(
            'Constants_attr_PROTONMASS',
        ),
        hpath='/RuntimePars',
        attr='REION_He_Heating_EVpH',
        convert=lambda vals, raw, path, fname, hpath:
        raw / vals.Constants_attr_PROTONMASS.to(U.g).value,
        units=U.eV * U.g ** -1,
        unit_convert=U.erg * U.g ** -1
    )

    extractors['RuntimePars_attr_REION_He_ZCenter'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='REION_He_ZCenter',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.dimensionless_unscaled,
        unit_convert=None
    )

    extractors['RuntimePars_attr_REION_He_ZSigma'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='REION_He_ZSigma',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.dimensionless_unscaled,
        unit_convert=None
    )

    extractors['RuntimePars_attr_SF_SchmidtLawCoeff_MSUNpYRpKPC2'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='SF_SchmidtLawCoeff_MSUNpYRpKPC2',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.Msun * U.yr ** -1 * U.kpc ** -2,
        unit_convert=None
    )

    extractors['RuntimePars_attr_SF_SchmidtLawExponent'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='SF_SchmidtLawExponent',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.dimensionless_unscaled,
        unit_convert=None
    )

    extractors['RuntimePars_attr_SF_SchmidtLawHighDensExponent'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='SF_SchmidtLawHighDensExponent',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.dimensionless_unscaled,
        unit_convert=None
    )

    extractors['RuntimePars_attr_SF_SchmidtLawHighDensThresh_HpCM3'] = \
        extractor(
            keytype='header',
            filetype='snapshot',
            dependencies=tuple(),
            hpath='/RuntimePars',
            attr='SF_SchmidtLawHighDensThresh_HpCM3',
            convert=lambda vals, raw, path, fname, hpath:
            raw,
            units=U.cm ** -3,
            unit_convert=None
        )

    extractors['RuntimePars_attr_SF_THRESH_MaxPhysDensOn'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='SF_THRESH_MaxPhysDensOn',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=None,
        unit_convert=None
    )

    extractors['RuntimePars_attr_SF_THRESH_MaxPhysDens_HpCM3'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='SF_THRESH_MaxPhysDens_HpCM3',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.cm ** -3,
        unit_convert=None
    )

    extractors['RuntimePars_attr_SF_THRESH_MetDepSFThreshMaxThresh_HpCM3'] = \
        extractor(
            keytype='header',
            filetype='snapshot',
            dependencies=tuple(),
            hpath='/RuntimePars',
            attr='SF_THRESH_MetDepSFThreshMaxThresh_HpCM3',
            convert=lambda vals, raw, path, fname, hpath:
            raw,
            units=U.cm ** -3,
            unit_convert=None
        )

    extractors['RuntimePars_attr_SF_THRESH_MetDepSFThreshNorm_HpCM3'] = \
        extractor(
            keytype='header',
            filetype='snapshot',
            dependencies=tuple(),
            hpath='/RuntimePars',
            attr='SF_THRESH_MetDepSFThreshNorm_HpCM3',
            convert=lambda vals, raw, path, fname, hpath:
            raw,
            units=U.cm ** -3,
            unit_convert=None
        )

    extractors['RuntimePars_attr_SF_THRESH_MetDepSFThreshSlope'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='SF_THRESH_MetDepSFThreshSlope',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.dimensionless_unscaled,
        unit_convert=None
    )

    extractors['RuntimePars_attr_SF_THRESH_MinOverDens'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='SF_THRESH_MinOverDens',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.dimensionless_unscaled,
        unit_convert=None
    )

    extractors['RuntimePars_attr_SF_THRESH_MinPhysDens_HpCM3'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='SF_THRESH_MinPhysDens_HpCM3',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.cm ** -3,
        unit_convert=None
    )

    extractors['RuntimePars_attr_SF_THRESH_TempMargin_DEX'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='SF_THRESH_TempMargin_DEX',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.dimensionless_unscaled,
        unit_convert=None
    )

    extractors['RuntimePars_attr_SNII_Delta_T_Divided_By_T_Vir'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='SNII_Delta_T_Divided_By_T_Vir',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.dimensionless_unscaled,
        unit_convert=None
    )

    extractors['RuntimePars_attr_SNII_Delta_T_K'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='SNII_Delta_T_K',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.K,
        unit_convert=None
    )

    extractors['RuntimePars_attr_SNII_EnergyTransferOn'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='SNII_EnergyTransferOn',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=None,
        unit_convert=None
    )

    extractors['RuntimePars_attr_SNII_Energy_ERG'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='SNII_Energy_ERG',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.erg,
        unit_convert=None
    )

    for element in elements:
        extractors['RuntimePars_attr_SNII_Factor_{:s}'.format(element)] = \
            extractor(
                keytype='header',
                filetype='snapshot',
                dependencies=tuple(),
                hpath='/RuntimePars',
                attr='SNII_Factor_{:s}'.format(element),
                convert=lambda vals, raw, path, fname, hpath:
                raw,
                units=None,
                unit_convert=None
            )

    extractors['RuntimePars_attr_SNII_MassTransferOn'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='SNII_MassTransferOn',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=None,
        unit_convert=None
    )

    extractors['RuntimePars_attr_SNII_MaxEnergyFraction'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='SNII_MaxEnergyFraction',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.dimensionless_unscaled,
        unit_convert=None
    )

    extractors['RuntimePars_attr_SNII_MaxMass_MSUN'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='SNII_MaxMass_MSUN',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.Msun,
        unit_convert=None
    )

    extractors['RuntimePars_attr_SNII_Max_Delta_T_K'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='SNII_Max_Delta_T_K',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.K,
        unit_convert=None
    )

    extractors['RuntimePars_attr_SNII_MinEnergyFraction'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='SNII_MinEnergyFraction',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.dimensionless_unscaled,
        unit_convert=None
    )

    extractors['RuntimePars_attr_SNII_MinMass_MSUN'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='SNII_MinMass_MSUN',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.Msun,
        unit_convert=None
    )

    extractors['RuntimePars_attr_SNII_Min_Delta_T_K'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='SNII_Min_Delta_T_K',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.K,
        unit_convert=None
    )

    extractors['RuntimePars_attr_SNII_Tvir0_K'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='SNII_Tvir0_K',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.K,
        unit_convert=None
    )

    extractors['RuntimePars_attr_SNII_Width_logTvir_dex'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='SNII_Width_logTvir_dex',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.dimensionless_unscaled,
        unit_convert=None
    )

    extractors['RuntimePars_attr_SNII_WindDelay_YR'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='SNII_WindDelay_YR',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.yr,
        unit_convert=None
    )

    extractors['RuntimePars_attr_SNII_WindIsotropicOn'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='SNII_WindIsotropicOn',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=None,
        unit_convert=None
    )

    extractors['RuntimePars_attr_SNII_exponent_Delta_T'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='SNII_exponent_Delta_T',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.dimensionless_unscaled,
        unit_convert=None
    )

    extractors['RuntimePars_attr_SNII_normalisation_Delta_T_K'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='SNII_normalisation_Delta_T_K',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.K,
        unit_convert=None
    )

    extractors['RuntimePars_attr_SNII_rhogas_physdensnormfac'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='SNII_rhogas_physdensnormfac',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.dimensionless_unscaled,
        unit_convert=None
    )

    extractors['RuntimePars_attr_SNII_rhogas_power'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='SNII_rhogas_power',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.dimensionless_unscaled,
        unit_convert=None
    )

    extractors['RuntimePars_attr_SNII_zdep_power'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='SNII_zdep_power',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.dimensionless_unscaled,
        unit_convert=None
    )

    extractors['RuntimePars_attr_SNIa_Efficiency'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='SNIa_Efficiency',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.dimensionless_unscaled,
        unit_convert=None
    )

    extractors['RuntimePars_attr_SNIa_EjectaVelocity_KMpS'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='SNIa_EjectaVelocity_KMpS',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.km * U.s ** -1,
        unit_convert=None
    )

    extractors['RuntimePars_attr_SNIa_EnergyFraction'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='SNIa_EnergyFraction',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.dimensionless_unscaled,
        unit_convert=None
    )

    extractors['RuntimePars_attr_SNIa_EnergyTransferOn'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='SNIa_EnergyTransferOn',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=None,
        unit_convert=None
    )

    extractors['RuntimePars_attr_SNIa_EnergyTransferStochastic'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='SNIa_EnergyTransferStochastic',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=None,
        unit_convert=None
    )

    extractors['RuntimePars_attr_SNIa_Energy_ERG'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='SNIa_Energy_ERG',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.erg,
        unit_convert=None
    )

    extractors['RuntimePars_attr_SNIa_MassTransferOn'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='SNIa_MassTransferOn',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=None,
        unit_convert=None
    )

    extractors['RuntimePars_attr_SNIa_Model'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='SNIa_Model',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=None,
        unit_convert=None
    )

    extractors['RuntimePars_attr_SNIa_TimeScale'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='SNIa_TimeScale',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.dimensionless_unscaled,
        unit_convert=None
    )

    extractors['RuntimePars_attr_SeedBlackHoleMass_Msun'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='SeedBlackHoleMass_Msun',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.Msun,
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

    extractors['RuntimePars_attr_StarformationOn'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='StarformationOn',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=None,
        unit_convert=None
    )

    extractors['RuntimePars_attr_StellarEnergyFeedbackOn'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='StellarEnergyFeedbackOn',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=None,
        unit_convert=None
    )

    extractors['RuntimePars_attr_StellarEvol_FeedbackOn'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='StellarEvol_FeedbackOn',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=None,
        unit_convert=None
    )

    extractors['RuntimePars_attr_StellarEvolutionCut_Gyr'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='StellarEvolutionCut_Gyr',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.Gyr,
        unit_convert=None
    )

    extractors['RuntimePars_attr_StellarEvolutionTimestepInterval'] = \
        extractor(
            keytype='header',
            filetype='snapshot',
            dependencies=tuple(),
            hpath='/RuntimePars',
            attr='StellarEvolutionTimestepInterval',
            convert=lambda vals, raw, path, fname, hpath:
            raw,
            units=U.dimensionless_unscaled,
            unit_convert=None
        )

    extractors['RuntimePars_attr_StellarMetalFeedbackOn'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='StellarMetalFeedbackOn',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=None,
        unit_convert=None
    )

    extractors['RuntimePars_attr_SulphurOverSilicon'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='SulphurOverSilicon',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.dimensionless_unscaled,
        unit_convert=None
    )

    extractors['RuntimePars_attr_TreeDomainUpdateFrequency'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='TreeDomainUpdateFrequency',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.dimensionless_unscaled,
        unit_convert=None
    )

    extractors['RuntimePars_attr_TypeOfOpeningCriterion'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='TypeOfOpeningCriterion',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=None,
        unit_convert=None
    )

    extractors['RuntimePars_attr_TypeOfTimestepCriterion'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='TypeOfTimestepCriterion',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=None,
        unit_convert=None
    )

    extractors['RuntimePars_attr_stellar_feedback_DeltaT'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='stellar_feedback_DeltaT',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=None,
        unit_convert=None
    )

    extractors['RuntimePars_attr_stellar_feedback_mode'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='stellar_feedback_mode',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=None,
        unit_convert=None
    )

    extractors['RuntimePars_attr_stellar_feedback_tvir'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/RuntimePars',
        attr='stellar_feedback_tvir',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=None,
        unit_convert=None
    )

    extractors['Units_attr_UnitDensity_in_cgs'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/Units',
        attr='UnitDensity_in_cgs',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.dimensionless_unscaled,
        unit_convert=None
    )

    extractors['Units_attr_UnitEnergy_in_cgs'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/Units',
        attr='UnitEnergy_in_cgs',
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

    extractors['Units_attr_UnitPressure_in_cgs'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/Units',
        attr='UnitPressure_in_cgs',
        convert=lambda vals, raw, path, fname, hpath:
        raw,
        units=U.dimensionless_unscaled,
        unit_convert=None
    )

    extractors['Units_attr_UnitTime_in_s'] = extractor(
        keytype='header',
        filetype='snapshot',
        dependencies=tuple(),
        hpath='/Units',
        attr='UnitTime_in_s',
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

    return extractors
