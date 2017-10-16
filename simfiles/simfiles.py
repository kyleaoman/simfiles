from hdf5_io import get as hdf5_get
import numpy as np
from vic_paths import victoria_LG_paths as lgpaths
import z_a_t_conversions as zat
from rahmati2013_neutral_frac import rahmati2013_neutral_frac as HI_frac, molecular_frac
import warnings
from collections import namedtuple
suffix = np.genfromtxt('/astro/koman/utilities/eagle_suffix.txt',dtype='string')

T = {
    'g': '0',
    'dm': '1',
    'b2': '2',
    'b3': '3',
    's': '4',
    'bh': '5'
}

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

softstrings = {
    'g': 'Gas', 
    'dm': 'Halo', 
    'b2': 'Disk', 
    'b3': 'Bulge', 
    's': 'Stars', 
    'bh': 'Bndry'
}

class ApostleFileset(dict):
    def __init__(self, res=None, vol=None, phys=None, snap=None):
        
        if (res not in ('L', 'M', 'H')):
            raise ValueError("ApostleFileset: res must be in 'L', 'M', 'H'")
        if (vol not in [str(i) for i in range(1,13)]):
            raise ValueError("ApostleFileset: vol must be in '1', '2', ..., '12'")
        if (phys not in ('fix', 'DMO')):
            raise ValueError("ApostleFileset: phys must be 'fix' or 'DMO'")
        if (type(snap) != int):
            raise ValueError('ApostleFileset: provide snapshot number as integer')

        self.res, self.vol, self.phys, self.snap = res, vol, phys, snap

        self._init_extractors()
        
        return

    def __setattr__(self, key, value):
        return self.__setitem__(key, value)

    def __getattr__(self, key):
        try:
            return self.__getitem__(key)
        except KeyError:
            raise AttributeError("'ApostleFileset' object has no attribute '"+str(key)+"'")
    
    def load(self, ftype, keys=tuple()):

        loaded_keys = set()

        if (ftype not in ('group', 'particle', 'snapshot')):
            raise ValueError("ApostleFileset: ftype must be in 'group', 'particle', 'snapshot'")
        if (type(keys) != tuple):
            raise ValueError('ApostleFileset: keys must be tuple')

        if ftype == 'group':
            path = lgpaths[self.res][self.phys][self.vol]+'/groups_'+suffix[self.snap]
            fname = 'eagle_subfind_tab_'+suffix[self.snap]
        if ftype == 'particle':
            path = lgpaths[self.res][self.phys][self.vol]+'/particledata_'+suffix[self.snap]
            fname = 'eagle_subfind_particles_'+suffix[self.snap]
        if ftype == 'snapshot':
            path = lgpaths[self.res][self.phys][self.vol]+'/snapshot_'+suffix[self.snap]
            fname = 'snap_'+suffix[self.snap]
        
        for key in keys:
            loaded_keys.update(self._load_key(path, fname, key))

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
        if E.attr is None:
            self[key] = E.convert(hdf5_get(path, fname, E.hpath), path, fname, E.hpath)
        else:
            self[key] = E.convert(hdf5_get(path, fname, E.hpath, attr=E.attr), path, fname, E.hpath)
            
        loaded_keys.update((key, ))
        
        return loaded_keys

    def _init_extractors(self):

        Extractor = namedtuple('Extractor', ['keytype', 'dependencies', 'hpath', 'attr', 'convert'])

        self._Extractors = {}

        h_a_powers = lambda path, fname, hpath: np.power(self.h, hdf5_get(path, fname, hpath, attr='h-scale-exponent')) * np.power(self.a, hdf5_get(path, fname, hpath, attr='aexp-scale-exponent'))

        #G
        self._Extractors['G'] = Extractor(
            keytype = 'header',
            dependencies = ('cm_per_kpc', 'solar_mass'),
            hpath = '/Constants',
            attr = 'GRAVITY',
            convert = lambda raw, path, fname, hpath: raw * self.solar_mass * np.power(1.E5, -2) / self.cm_per_kpc 
        )

        #kB
        self._Extractors['kB'] = Extractor(
            keytype = 'header',
            dependencies = tuple(),
            hpath = '/Constants',
            attr = 'BOLTZMANN',
            convert = lambda raw, path, fname, hpath: raw
        )
        
        #a
        self._Extractors['a'] = Extractor(
            keytype = 'header',
            dependencies = tuple(),
            hpath = '/Header',
            attr = 'Time',
            convert = lambda raw, path, fname, hpath: raw
        )

        #h
        self._Extractors['h'] = Extractor(
            keytype = 'header',
            dependencies = tuple(),
            hpath = '/Header',
            attr = 'HubbleParam',
            convert = lambda raw, path, fname, hpath: raw
        )

        #code_to_g
        self._Extractors['code_to_g'] = Extractor(
            keytype = 'header',
            dependencies = tuple(),
            hpath = '/Units',
            attr = 'UnitMass_in_g',
            convert = lambda raw, path, fname, hpath: raw
        )

        #code_to_cm
        self._Extractors['code_to_cm'] = Extractor(
            keytype = 'header',
            dependencies = tuple(),
            hpath = '/Units',
            attr = 'UnitLength_in_cm',
            convert = lambda raw, path, fname, hpath: raw
        )

        #code_to_Msun
        self._Extractors['code_to_Msun'] = Extractor(
            keytype = 'header',
            dependencies = ('solar_mass',),
            hpath = '/Units',
            attr = 'UnitMass_in_g',
            convert = lambda raw, path, fname, hpath: raw / self.solar_mass
        )

        #solar_mass
        self._Extractors['solar_mass'] = Extractor(
            keytype = 'header',
            dependencies = tuple(),
            hpath = '/Constants',
            attr = 'SOLAR_MASS',
            convert = lambda raw, path, fname, hpath: raw
        )

        #code_to_kpc
        self._Extractors['code_to_kpc'] = Extractor(
            keytype = 'header',
            dependencies = ('cm_per_kpc',),
            hpath = '/Units',
            attr = 'UnitLength_in_cm',
            convert = lambda raw, path, fname, hpath: raw / self.cm_per_kpc
        )

        #cm_per_kpc
        self._Extractors['cm_per_kpc'] = Extractor(
            keytype = 'header',
            dependencies = tuple(),
            hpath = '/Constants',
            attr = 'CM_PER_MPC',
            convert = lambda raw, path, fname, hpath: raw / 1.E3
        )

        #code_to_km_s
        self._Extractors['code_to_km_s'] = Extractor(
            keytype = 'header',
            dependencies = tuple(),
            hpath = '/Units',
            attr = 'UnitVelocity_in_cm_per_s',
            convert = lambda raw, path, fname, hpath: raw / 1.E5
        )

        #Lbox
        self._Extractors['Lbox'] = Extractor(
            keytype = 'header',
            dependencies = ('code_to_kpc', 'h'),
            hpath = '/Header',
            attr = 'BoxSize',
            convert = lambda raw, path, fname, hpath: raw * self.code_to_kpc / self.h
        )

        #proton_mass
        self._Extractors['proton_mass'] = Extractor(
            keytype = 'header',
            dependencies = tuple(),
            hpath = '/Constants',
            attr = 'PROTONMASS',
            convert = lambda raw, path, fname, hpath: raw
        )

        #redshift
        self._Extractors['redshift'] = Extractor(
            keytype = 'header',
            dependencies = tuple(),
            hpath = '/Header',
            attr = 'Redshift',
            convert = lambda raw, path, fname, hpath: raw
        )

        #fH
        self._Extractors['fH'] = Extractor(
            keytype = 'header',
            dependencies = tuple(),
            hpath = '/RuntimePars',
            attr = 'InitAbundance_Hydrogen',
            convert = lambda raw, path, fname, hpath: raw
        )

        #fHe
        self._Extractors['fHe'] = Extractor(
            keytype = 'header',
            dependencies = tuple(),
            hpath = '/RuntimePars',
            attr = 'InitAbundance_Helium',
            convert = lambda raw, path, fname, hpath: raw
        )

        #mu
        self._Extractors['mu'] = Extractor(
            keytype = 'header',
            dependencies = ('fHe',),
            hpath = '/RuntimePars',
            attr = 'InitAbundance_Hydrogen',
            convert = lambda raw, path, fname, hpath: 1. / (raw + .25 * self.fHe)
        )

        #gamma
        self._Extractors['gamma'] = Extractor(
            keytype = 'header',
            dependencies = tuple(),
            hpath = '/RuntimePars',
            attr = 'EOS_Jeans_GammaEffective',
            convert = lambda raw, path, fname, hpath: raw
        )

        #T0
        self._Extractors['T0'] = Extractor(
            keytype = 'header',
            dependencies = tuple(),
            hpath = '/RuntimePars',
            attr = 'EOS_Jeans_TempNorm_K',
            convert = lambda raw, path, fname, hpath: raw
        )

        #p_mass
        self._Extractors['p_mass'] = Extractor(
            keytype = 'header',
            dependencies = ('code_to_Msun', 'h'),
            hpath = '/Header',
            attr = 'MassTable',
            convert = lambda raw, path, fname, hpath: self.code_to_Msun / self.h * raw[1]
        )

        #eps_g
        self._Extractors['eps_g'] = Extractor(
            keytype = 'header',
            dependencies = ('code_to_kpc', 'h', 'a', 'eps_maxphys_g'),
            hpath = '/RuntimePars',
            attr = 'Softening' + softstrings['g'],
            convert = lambda raw, path, fname, hpath: min(raw * self.code_to_kpc * self.a / self.h, self.eps_maxphys_g)
        )

        #eps_dm
        self._Extractors['eps_dm'] = Extractor(
            keytype = 'header',
            dependencies = ('code_to_kpc', 'h', 'a', 'eps_maxphys_dm'),
            hpath = '/RuntimePars',
            attr = 'Softening' + softstrings['dm'],
            convert = lambda raw, path, fname, hpath: min(raw * self.code_to_kpc * self.a / self.h, self.eps_maxphys_dm)
        )

        #eps_b2
        self._Extractors['eps_b2'] = Extractor(
            keytype = 'header',
            dependencies = ('code_to_kpc', 'h', 'a', 'eps_maxphys_b2'),
            hpath = '/RuntimePars',
            attr = 'Softening' + softstrings['b2'],
            convert = lambda raw, path, fname, hpath: min(raw * self.code_to_kpc * self.a / self.h, self.eps_maxphys_b2)
        )

        #eps_b3
        self._Extractors['eps_b3'] = Extractor(
            keytype = 'header',
            dependencies = ('code_to_kpc', 'h', 'a', 'eps_maxphys_b3'),
            hpath = '/RuntimePars',
            attr = 'Softening' + softstrings['b3'],
            convert = lambda raw, path, fname, hpath: min(raw * self.code_to_kpc * self.a / self.h, self.eps_maxphys_b3)
        )

        #eps_s
        self._Extractors['eps_s'] = Extractor(
            keytype = 'header',
            dependencies = ('code_to_kpc', 'h', 'a', 'eps_maxphys_s'),
            hpath = '/RuntimePars',
            attr = 'Softening' + softstrings['s'],
            convert = lambda raw, path, fname, hpath: min(raw * self.code_to_kpc * self.a / self.h, self.eps_maxphys_s)
        )

        #eps_bh
        self._Extractors['eps_bh'] = Extractor(
            keytype = 'header',
            dependencies = ('code_to_kpc', 'h', 'a', 'eps_maxphys_bh'),
            hpath = '/RuntimePars',
            attr = 'Softening' + softstrings['bh'],
            convert = lambda raw, path, fname, hpath: min(raw * self.code_to_kpc * self.a / self.h, self.eps_maxphys_bh)
        )

        #eps_maxphys_*
        for ptype in T.keys():
            self._Extractors['eps_maxphys_' + ptype] = Extractor(
                keytype = 'header',
                dependencies = ('code_to_kpc', 'h', 'a'),
                hpath = '/RuntimePars',
                attr = 'Softening' + softstrings[ptype] + 'MaxPhys',
                convert = lambda raw, path, fname, hpath: raw * self.code_to_kpc / self.h
            )

        #contamination
        self._Extractors['contamination'] = Extractor(
            keytype = 'group',
            dependencies = tuple(),
            hpath = '/FOF/ContaminationCount',
            attr = None,
            convert = lambda raw, path, fname, hpath: raw
        )

        #nsubhalos
        self._Extractors['nsubhalos'] = Extractor(
            keytype = 'group',
            dependencies = tuple(),
            hpath = '/FOF/NumOfSubhalos',
            attr = None,
            convert = lambda raw, path, fname, hpath: raw
        )
            
        #gns
        self._Extractors['gns'] = Extractor(
            keytype = 'group',
            dependencies = tuple(),
            hpath = '/Subhalo/GroupNumber',
            attr = None,
            convert = lambda raw, path, fname, hpath: raw
        )

        #sgns
        self._Extractors['sgns'] = Extractor(
            keytype = 'group',
            dependencies = tuple(),
            hpath = '/Subhalo/SubGroupNumber',
            attr = None,
            convert = lambda raw, path, fname, hpath: raw
        )

        #cops
        self._Extractors['cops'] = Extractor(
            keytype = 'group',
            dependencies = ('code_to_kpc', 'h', 'a'),
            hpath = '/Subhalo/CentreOfPotential',
            attr = None,
            convert = lambda raw, path, fname, hpath: raw * self.code_to_kpc * h_a_powers(path, fname, hpath)
        )

        #vcents
        self._Extractors['vcents'] = Extractor(
            keytype = 'group',
            dependencies = ('code_to_km_s', 'h', 'a'),
            hpath = '/Subhalo/Velocity',
            attr = None,
            convert = lambda raw, path, fname, hpath: raw * self.code_to_km_s * h_a_powers(path, fname, hpath)
        )
        
        #nID
        self._Extractors['nID'] = Extractor(
            keytype = 'group',
            dependencies = tuple(),
            hpath = '/Subhalo/SubLength',
            attr = None,
            convert = lambda raw, path, fname, hpath: raw
        )

        #offID
        self._Extractors['offID'] = Extractor(
            keytype = 'group',
            dependencies = tuple(),
            hpath = '/Subhalo/SubOffset',
            attr = None,
            convert = lambda raw, path, fname, hpath: raw
        )

        #msubfind_*
        for ptype in T.keys():
            self._Extractors['msubfind_' + ptype] = Extractor(
                keytype = 'group',
                dependencies = ('code_to_Msun', 'h', 'a'),
                hpath = '/Subhalo/MassType',
                attr = None,
                convert = lambda raw, path, fname, hpath: raw[:, int(T[ptype])] * self.code_to_Msun * h_a_powers(path, fname, hpath)
            )

        #nfof
        self._Extractors['nfof'] = Extractor(
            keytype = 'header',
            dependencies = tuple(),
            hpath = '/FOF',
            attr = 'TotNgroups',
            convert = lambda raw, path, fname, hpath: raw
        )
            
        #M200
        self._Extractors['M200'] = Extractor(
            keytype = 'fofgroup',
            dependencies = ('code_to_Msun', 'h', 'a'),
            hpath = '/FOF/Group_M_Crit200',
            attr = None,
            convert = lambda raw, path, fname, hpath: raw * self.code_to_Msun * h_a_powers(path, fname, hpath)
        )

        #R200
        self._Extractors['R200'] = Extractor(
            keytype = 'fofgroup',
            dependencies = ('code_to_kpc', 'h', 'a'),
            hpath = '/FOF/Group_R_Crit200',
            attr = None,
            convert = lambda raw, path, fname, hpath: raw * self.code_to_kpc * h_a_powers(path, fname, hpath)
        )

        #ids
        self._Extractors['ids'] = Extractor(
            keytype = 'idgroup',
            dependencies = tuple(),
            hpath = '/IDs/ParticleID',
            attr = None,
            convert = lambda raw, path, fname, hpath: raw
        )

        #ids_*
        for ptype in T.keys():
            self._Extractors['ids_' + ptype] = Extractor(
                keytype = 'particle',
                dependencies = tuple(),
                hpath = '/PartType' + T[ptype] + '/ParticleIDs',
                attr = None,
                convert = lambda raw, path, fname, hpath: raw
            )

        #xyz_*
        for ptype in T.keys():
            self._Extractors['xyz_' + ptype] = Extractor(
                keytype = 'particle',
                dependencies = ('code_to_kpc', 'h', 'a'),
                hpath = '/PartType' + T[ptype] + '/Coordinates',
                attr = None,
                convert = lambda raw, path, fname, hpath: raw * self.code_to_kpc * h_a_powers(path, fname, hpath)
            )

        #vxyz_*
        for ptype in T.keys():
            self._Extractors['vxyz_' + ptype] = Extractor(
                keytype = 'particle',
                dependencies = ('code_to_km_s', 'h', 'a'),
                hpath = '/PartType' + T[ptype] + '/Velocities',
                attr = None,
                convert = lambda raw, path, fname, hpath: raw * self.code_to_km_s * h_a_powers(path, fname, hpath)
            )

        #ng_*
        for ptype in T.keys():
            self._Extractors['ng_' + ptype] = Extractor(
                keytype = 'particle',
                dependencies = tuple(),
                hpath = '/PartType' + T[ptype] + '/GroupNumber',
                attr = None,
                convert = lambda raw, path, fname, hpath: raw
            )

        #nsg_*
        for ptype in T.keys():
            self._Extractors['nsg_'+ptype] = Extractor(
                keytype = 'particle',
                dependencies = tuple(),
                hpath = '/PartType' + T[ptype] + '/SubGroupNumber',
                attr = None,
                convert = lambda raw, path, fname, hpath: raw
            )

        #m_g, m_b2, m_b3, m_s, m_bh
        for ptype in ['g', 'b2', 'b3', 's', 'bh']:
            self._Extractors['m_' + ptype] = Extractor(
                keytype = 'particle',
                dependencies = ('code_to_Msun', 'h', 'a'),
                hpath = '/PartType' + T[ptype] + '/Masses',
                attr = None,
                convert = lambda raw, path, fname, hpath: raw * self.code_to_Msun * h_a_powers(path, fname, hpath)
            )

        #m_dm
        self._Extractors['m_dm'] = Extractor(
            keytype = 'particle',
            dependencies = ('p_mass',),
            hpath = '/PartType1/ParticleIDs',
            attr = None,
            convert = lambda raw, path, fname, hpath: np.ones(raw.shape, dtype=np.float) * self.p_mass
        )

        #T_g
        self._Extractors['T_g'] = Extractor(
            keytype = 'particle',
            dependencies = ('h', 'a'),
            hpath = '/PartType0/Temperature',
            attr = None,
            convert = lambda raw, path, fname, hpath: raw * h_a_powers(path, fname, hpath)
        )

        #rho_g
        self._Extractors['rho_g'] = Extractor(
            keytype = 'particle',
            dependencies = ('code_to_g', 'code_to_cm', 'h', 'a'),
            hpath = '/PartType0/Density',
            attr = None,
            convert = lambda raw, path, fname, hpath: raw * self.code_to_g * np.power(self.code_to_cm, -3) * h_a_powers(path, fname, hpath) 
        )

        #*abundance_g, *abundance_s, sm*abundance_g, sm*abundance_s
        for typesuffix in ['g', 's']:
            for element in elements.keys():
                for prefix, smooth in {'sm': 'Smoothed', '': ''}.items():
                    self._Extractors[prefix+element+'abundance_'+typesuffix] = Extractor(
                        keytype = 'particle',
                        dependencies = ('h', 'a'),
                        hpath = '/PartType'+T[typesuffix]+'/'+smooth+elements[element],
                        attr = None,
                        convert = lambda raw, path, fname, hpath: raw * h_a_powers(path, fname, hpath)
                    )
        #SFR_g
        self._Extractors['SFR_g'] = Extractor(
            keytype = 'particle',
            dependencies = ('h', 'a'),
            hpath = '/PartType0/StarFormationRate',
            attr = None,
            convert = lambda raw, path, fname, hpath: raw * h_a_powers(path, fname, hpath)
        )

        #hsm_g, hsm_s
        for typesuffix in ['g', 's']:
            self._Extractors['hsm_' + typesuffix] = Extractor(
                keytype = 'particle',
                dependencies = ('code_to_kpc', 'h', 'a'),
                hpath = '/PartType'+T[typesuffix]+'/SmoothingLength',
                attr = None,
                convert = lambda raw, path, fname, hpath: raw * self.code_to_kpc * h_a_powers(path, fname, hpath)
            )

        #age_s
        self._Extractors['age_s'] = Extractor(
            keytype = 'particle',
            dependencies = tuple(),
            hpath = '/PartType4/StellarFormationTime',
            attr = None,
            convert = lambda raw, path, fname, hpath: zat.a_to_lb(raw)
        )

        #mHI_g
        self._Extractors['mHI_g'] = Extractor(
            keytype = 'particle',
            dependencies = ('redshift', 'rho_g', 'Habundance_g', 'proton_mass', 'SFR_g', 'fH', 'T_g', 'code_to_Msun', 'mu', 'T0', 'gamma'),
            hpath = '/PartType0/Masses',
            attr = None,
            convert = lambda raw, path, fname, hpath: raw * self.code_to_Msun * self.Habundance_g * h_a_powers(path, fname, hpath) * (1. - molecular_frac(self.SFR_g, self.T_g, self.rho_g, self.Habundance_g, mu=self.mu, proton_mass=self.proton_mass, gamma=self.gamma, fH=self.fH, T0=self.T0)) * HI_frac(self.redshift, self.rho_g * self.Habundance_g / (self.mu * self.proton_mass), self.T_g, onlyA1=True, APOSTLE_corrections=True, SFR=self.SFR_g, mu=self.mu, proton_mass=self.proton_mass, gamma=self.gamma, fH=self.fH, Habundance=self.Habundance_g, T0=self.T0, rho=self.rho_g)
        )
