from simfiles._setup_cfg import snapshots, extractors
from simfiles.configs._EAGLE_extractors import generate_eagle_extractors, \
    generate_extra_extractors
from collections import namedtuple
from os.path import dirname, join
import glob
import re

with open(join(dirname(__file__), 'machine')) as mfile:
    machine = mfile.read().strip()

snap_id = namedtuple('snap_id', ['box', 'res', 'model', 'snap'])

path_bases = {
    'cosma': '/cosma7/data/Eagle/ScienceRuns/Planck1',
}

aliasfile = 'EAGLE.alias'

boxes = {
    'L0034': {
        'N1034': ['DMONLY', 'RECALIBRATED'],
    },
}

box_list = [(box, res, model) for box, v in boxes.items()
            for res, vv in v.items() for model in vv]

for box, res, model in box_list:
    if model == 'RECALIBRATED':
        path_prefix = '{:s}/{:s}{:s}_WDM_L11p2/{:s}/data'.format(
            path_bases[machine], box, res, model)
    else:
        path_prefix = '{:s}/{:s}{:s}_WDM_L11p2/PE/{:s}/data'.format(
            path_bases[machine], box, res, model)
    suffixes = [
        '_'.join(s.split('/')[-1].split('_')[1:])
        for s in glob.glob('{:s}/snapshot_*'.format(path_prefix))
        if bool(re.match("[0-9]{3,4}_z[0-9]{3}p[0-9]{3}",
                         '_'.join(s.split('/')[-1].split('_')[1:])))
    ]
    for suffix in suffixes:
        snap = int(suffix.split('_')[0])
        group_path = '{:s}/groups_{:s}'.format(path_prefix, suffix)
        group_file = 'eagle_subfind_tab_{:s}'.format(suffix)
        fof_file = 'group_tab_{:s}'.format(suffix)
        particle_path = '{:s}/particledata_{:s}'.format(path_prefix, suffix)
        particle_file = 'eagle_subfind_particles_{:s}'.format(suffix)
        snapshot_path = '{:s}/snapshot_{:s}'.format(path_prefix, suffix)
        snapshot_file = 'snap_{:s}'.format(suffix)

        snapshots[snap_id(box=box, res=res, model=model, snap=snap)] = {
            'group': (group_path, group_file),
            'fof': (group_path, fof_file),
            'particle': (particle_path, particle_file),
            'snapshot': (snapshot_path, snapshot_file),
        }

extractors.update(generate_eagle_extractors())
extractors.update(generate_extra_extractors())
