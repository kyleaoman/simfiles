from simfiles._setup_cfg import snapshots, extractors
from simfiles.configs._EAGLE_extractors import generate_eagle_extractors, \
    generate_extra_extractors
from collections import namedtuple
from os.path import dirname, join

with open(join(dirname(__file__), 'machine')) as mfile:
    machine = mfile.read()

suffix29 = [
    '000_z020p000', '001_z015p132', '002_z009p993', '003_z008p988',
    '004_z008p075', '005_z007p050', '006_z005p971', '007_z005p487',
    '008_z005p037', '009_z004p485', '010_z003p984', '011_z003p528',
    '012_z003p017', '013_z002p478', '014_z002p237', '015_z002p012',
    '016_z001p737', '017_z001p487', '018_z001p259', '019_z001p004',
    '020_z000p865', '021_z000p736', '022_z000p615', '023_z000p503',
    '024_z000p366', '025_z000p271', '026_z000p183', '027_z000p101',
    '028_z000p000'
]

suffix128 = [
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

snap_id = namedtuple('snap_id', ['model', 'snap'])

path_bases = {
    'cosma': {
        'ht_lt': '/cosma6/data/dp004/Eagle/Eagle_Cores',
        'sidm': '/cosma7/data/dp004/ndcf31/EAGLE_SIDM/HighSFthreshold',
        'wdm': '/cosma6/data/dp004/dc-foro1/run_eagle/wdm'
    }
}

aliasfile = 'EAGLE.alias'

boxes = {
    'ht_lt': [
        'DMONLY',
        'RECAL',
        'RECAL_EOS_10_NOAGN',
        'RECAL_EOS_50',
        'RECAL_EOS_50_NOAGN'
    ],
    'sidm': [
        'CDM_DMO',
        'CDM_RECAL',
        'CDM_SFT10_NOAGN',
        'CDM_SFT50',
        'SIDM10_DMO',
        'SIDM10_RECAL',
        'SIDM10_SFT10',
        'SIDM1_DMO',
        'SIDM1_RECAL',
        'vdSIDM_3.04_w560_DMO',
        'vdSIDM_3.04_w560_RECAL'
    ],
    'wdm': [
        'dmo',
        'recal',
        'recal_10'
    ]
}

for modelset, models in boxes.items():
    for model in models:
        if (modelset == 'ht_lt' and model == 'RECAL_EOS_10_NOAGN') or \
           (modelset == 'sidm' and model in (
               'CDM_SFT_NOAGN', 'SIDM10_DMO', 'SIDM10_RECAL',
               'SIDM10_SFT10', 'SIDM1_DMO', 'SIDM1_RECAL',
               'vdSIDM_3.04_w560_DMO', 'vdSIDM_3.04_w560_RECAL'
           )) or (modelset == 'wdm'):
            snaps = 128
            suffix = suffix128
        else:
            snaps = 29
            suffix = suffix29
        path_prefix = '{:s}/{:s}'.format(
            path_bases[machine][modelset], model)
        if modelset in ('sidm', 'wdm'):
            path_prefix = '{:s}/data'.format(path_prefix)
        for snap in range(snaps):
            group_path = '{:s}/groups_{:s}'.format(
                path_prefix, suffix[snap])
            group_file = 'eagle_subfind_tab_{:s}'.format(
                suffix[snap])
            fof_file = 'group_tab_{:s}'.format(
                suffix[snap])
            particle_path = '{:s}/particledata_{:s}'.format(
                path_prefix, suffix[snap])
            particle_file = 'eagle_subfind_particles_{:s}'.format(
                suffix[snap])
            snapshot_path = '{:s}/snapshot_{:s}'.format(
                path_prefix, suffix[snap])
            snapshot_file = 'snap_{:s}'.format(
                suffix[snap])

            snapshots[snap_id(model=model, snap=snap)] = {
                'group': (group_path, group_file),
                'fof': (group_path, fof_file),
                'snapshot': (snapshot_path, snapshot_file),
                'particle': (particle_path, particle_file)
            }

extractors.update(generate_eagle_extractors(default_pfiletype='particle'))
extractors.update(generate_extra_extractors(default_pfiletype='particle'))
