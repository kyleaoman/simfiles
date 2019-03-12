from collections import namedtuple

snapshots = {}  # intialize snapshot dict

# define extractor fields
extractor = namedtuple(
    'extractor',
    ['keytype', 'filetype', 'dependencies', 'hpath', 'attr', 'convert',
     'units', 'unit_convert'])

extractors = {}  # initialize extractor dict
