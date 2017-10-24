from collections import namedtuple

snapshots = {} #intialize snapshot dict

extractor = namedtuple('extractor', ['keytype', 'filetype', 'dependencies', 'hpath', 'attr', 'convert', 'units', 'unit_convert']) #define extractor fields

extractors = {} #initialize extractor dict
