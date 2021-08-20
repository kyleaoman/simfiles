from simfiles._setup_cfg import snapshots, extractors
from simfiles.configs._EAGLE_extractors import generate_eagle_extractors, \
    generate_extra_extractors
from collections import namedtuple
from itertools import product
from os.path import dirname, join
import glob
import re

with open(join(dirname(__file__), 'machine')) as mfile:
    machine = mfile.read().strip()

# snap_id and snip_id will be identical as far as the snapshot dict indices
# are concerned, so we need an extra flag
snap_id = namedtuple('snap_id', ['res', 'phys', 'vol', 'snap'])
snip_id = namedtuple('snip_id', ['res', 'phys', 'vol', 'snip', 'is_snip'])

path_bases = {
    'cosma': '/cosma6/data/dp004/lg/snapshots_all/',
    'cavi': '/sraid14/azadehf/LG/data_fix/'
}

aliasfile = 'APOSTLE.alias'

res_str = {1: 'HR', 2: 'MR', 3: 'LR'}
vol_str = {1: 'V1', 2: 'V2', 3: 'V3', 4: 'V4', 5: 'V5', 6: 'V6',
           7: 'S1', 8: 'S2', 9: 'S3', 10: 'S4', 11: 'S5', 12: 'S6'}
phys_str = {'hydro': 'fix', 'DMO': 'DMO'}

for res, vol, phys in product(range(1, 4), range(1, 13), ['hydro', 'DMO']):

    if (res == 1) \
       and (
           ((phys == 'hydro') and (vol not in [1, 4, 6, 10, 11]))
           or
           ((phys == 'DMO') and (vol not in [1, 4, 11]))
       ):
        continue

    path_prefix = '{:s}{:s}_{:s}_{:s}'.format(
        path_bases[machine], vol_str[vol], res_str[res], phys_str[phys])
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

        snapshots[snap_id(res=res, phys=phys, vol=vol, snap=snap)] = {
            'group': (group_path, group_file),
            'fof': (group_path, fof_file),
            'particle': (particle_path, particle_file),
            'snapshot': (snapshot_path, snapshot_file),
        }

    snip_suffixes = [
        '_'.join(s.split('/')[-1].split('_')[1:])
        for s in glob.glob('{:s}/snipshot_*'.format(path_prefix))
        if bool(re.match("[0-9]{3,4}_z[0-9]{3}p[0-9]{3}",
                         '_'.join(s.split('/')[-1].split('_')[1:])))
    ]
    for suffix in snip_suffixes:
        snip = int(suffix.split('_')[0])
        group_path = '{:s}/groups_snip_{:s}'.format(path_prefix, suffix)
        group_file = None
        fof_file = 'group_tab_{:s}'.format(suffix)
        particle_path = None
        particle_file = None
        snapshot_path = '{:s}/snipshot_{:s}'.format(path_prefix, suffix)
        snapshot_file = 'snip_{:s}'.format(suffix)
        snapshots[
            snip_id(res=res, phys=phys, vol=vol, snip=snip, is_snip=True)
        ] = {
            'group': (group_path, group_file),
            'fof': (group_path, fof_file),
            'particle': (particle_path, particle_file),
            'snapshot': (snapshot_path, snapshot_file)
        }

extractors.update(generate_eagle_extractors(
    T=range(6),
    Mstring='Masses',
    Vstring='Velocities',
    EOSstring='SfFlag',
    omit=(
        'Header_attr_ExpansionFactor',
        'RuntimePars_attr_BH_MaxRepositionDistanceFactor'
    )
))
extractors.update(generate_extra_extractors(
    T=range(6),
    Mstring='Masses',
    Vstring='Velocities',
    EOSstring='SfFlag',
))
