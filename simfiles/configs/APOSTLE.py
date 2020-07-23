from simfiles._setup_cfg import snapshots, extractors
from simfiles.configs._EAGLE_extractors import generate_eagle_extractors
from collections import namedtuple
from itertools import product
from os.path import dirname, join

with open(join(dirname(__file__), 'machine')) as mfile:
    machine = mfile.read()

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

snap_id = namedtuple('snap_id', ['res', 'phys', 'vol', 'snap'])

path_bases = {
    'cosma': '/cosma6/data/dp004/lg/snapshots_all/',
    'cavi': '/sraid14/azadehf/LG/data_fix/'
}
res_str = {1: 'HR', 2: 'MR', 3: 'LR'}
vol_str = {1: 'V1', 2: 'V2', 3: 'V3', 4: 'V4', 5: 'V5', 6: 'V6',
           7: 'S1', 8: 'S2', 9: 'S3', 10: 'S4', 11: 'S5', 12: 'S6'}
phys_str = {'hydro': 'fix', 'DMO': 'DMO'}

for res, vol, phys, snap in product(
        range(1, 4), range(1, 13), ['hydro', 'DMO'], range(128)):

    path_prefix = '{:s}{:s}_{:s}_{:s}'.format(
        path_bases[machine], vol_str[vol], res_str[res], phys_str[phys])

    group_path = '{:s}/groups_{:s}'.format(path_prefix, suffix[snap])
    group_file = 'eagle_subfind_tab_{:s}'.format(suffix[snap])
    particle_path = '{:s}/particledata_{:s}'.format(path_prefix, suffix[snap])
    particle_file = 'eagle_subfind_particles_{:s}'.format(suffix[snap])
    snapshot_path = '{:s}/snapshot_{:s}'.format(path_prefix, suffix[snap])
    snapshot_file = 'snap_{:s}'.format(suffix[snap])

    if (res == 1) \
       and (
           ((phys == 'hydro') and (vol not in [1, 4, 6, 10, 11]))
           or
           ((phys == 'DMO') and (vol not in [1, 4, 11]))
       ):
        continue

    snapshots[snap_id(res=res, phys=phys, vol=vol, snap=snap)] = {
        'group': (group_path, group_file),
        'particle': (particle_path, particle_file),
        'snapshot': (snapshot_path, snapshot_file),
    }

extractors.update(generate_eagle_extractors(
    T=range(6),
    Mstring='Masses',
    Vstring='Velocities'
))
