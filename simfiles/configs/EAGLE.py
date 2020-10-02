from simfiles._setup_cfg import snapshots, extractors
from simfiles.configs._EAGLE_extractors import generate_eagle_extractors, \
    generate_extra_extractors
from collections import namedtuple
from itertools import product
from os.path import dirname, join

with open(join(dirname(__file__), 'machine')) as mfile:
    machine = mfile.read()

suffix = [
    '000_z020p000', '001_z015p132', '002_z009p993', '003_z008p988',
    '004_z008p075', '005_z007p050', '006_z005p971', '007_z005p487',
    '008_z005p037', '009_z004p485', '010_z003p984', '011_z003p528',
    '012_z003p017', '013_z002p478', '014_z002p237', '015_z002p012',
    '016_z001p737', '017_z001p487', '018_z001p259', '019_z001p004',
    '020_z000p865', '021_z000p736', '022_z000p615', '023_z000p503',
    '024_z000p366', '025_z000p271', '026_z000p183', '027_z000p101',
    '028_z000p000'
]

snap_id = namedtuple('snap_id', ['box', 'res', 'model', 'snap'])

path_bases = {
    'cosma': '/cosma7/data/Eagle/ScienceRuns/Planck1'
}

aliasfile = 'EAGLE.alias'

boxes = {
    'L0012': {
        'N0188': ['DMONLY', 'REFERENCE'],
        'N0376': ['DMONLY', 'REFERENCE', 'RECALIBRATED']
    },
    'L0025': {
        'N0188': ['DMONLY', 'REFERENCE'],
        'N0376': ['DMONLY', 'REFERENCE'],
        'N0752': ['DMONLY', 'REFERENCE', 'RECALIBRATED']
    },
    'L0050': {
        'N0376': ['DMONLY', 'REFERENCE'],
        'N0752': ['DMONLY', 'REFERENCE']
    },
    'L0100': {
        'N0752': ['DMONLY', 'REFERENCE'],
        'N1504': ['DMONLY', 'REFERENCE']
    }
}

box_list = [(box, res, model) for box, v in boxes.items()
            for res, vv in v.items() for model in vv]

for (box, res, model), snap in product(box_list, range(29)):
    path_prefix = '{:s}/{:s}{:s}/PE/{:s}/data'.format(
        path_bases[machine], box, res, model)

    group_path = '{:s}/groups_{:s}'.format(path_prefix, suffix[snap])
    group_file = 'eagle_subfind_tab_{:s}'.format(suffix[snap])
    fof_file = 'group_tab_{:s}'.format(suffix[snap])
    particle_path = '{:s}/particledata_{:s}'.format(path_prefix, suffix[snap])
    particle_file = 'eagle_subfind_particles_{:s}'.format(suffix[snap])
    snapshot_path = '{:s}/snapshot_{:s}'.format(path_prefix, suffix[snap])
    snapshot_file = 'snap_{:s}'.format(suffix[snap])

    snapshots[snap_id(box=box, res=res, model=model, snap=snap)] = {
        'group': (group_path, group_file),
        'fof': (group_path, fof_file),
        'particle': (particle_path, particle_file),
        'snapshot': (snapshot_path, snapshot_file),
    }

extractors.update(generate_eagle_extractors())
extractors.update(generate_extra_extractors())
