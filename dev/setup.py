'''
TODO: This will ultimately move into __init__.py in app's base directory.
TODO: throw verbose errors if config file not found; if vocabulary file or dict not found, look for Excel spreadsheet and pass to vocab_from_xls if found; if neither found, throw error.
Name:           cfizer.py
Author:         Cameron Wilson, CEMAC, University of Leeds
Date:           August 2023
CFizer Version: See __version__

'''

import yaml
import os
from datetime import datetime
import numpy as np
# from utils import vocab_from_xls

# WILDCARDS = {'*', '?'}
verbose = False
quiet = False

COMPRESSION = (True, 'zlib', 1) # True/False, compression type, complevel (1:9)
    # See https://unidata.github.io/netcdf4-python/efficient-compression-of-netcdf-variables
'''
Only lossless compression is considered here; all precision is maintained.
WARNING: only tested with zlib compression.
'''
if COMPRESSION[0]:
    if COMPRESSION[2] < 1 or COMPRESSION[2] > 9:
        raise ValueError("Only compression levels from 1 to 9 are valid.")
    if COMPRESSION[1] not in {
        'zlib', 'szip', 'zstd', 'bzip2', 'blosc_lz', 'blosc_lz4', 
        'blosc_lz4hc', 'blosc_zlib', 'blosc_zstd'
    }:
        raise ValueError("Invalid compression algorithm specified.")

CF_VERSION = '1.10'

DIM_GROUPS = {
    0: '0+1d', 
    1: '0+1d', 
    2: '2d', 
    3: '3d'
}

DIM_ACTIONS = {
    '0+1d': 'merge',
    '2d': 'merge', 
    '3d': 'split'
    }

default_time_unit = 's'

CF_ATTRIBUTES = {'title', 
                 'institution', 
                 'source', 
                 'history', 
                 'references', 
                 'comment',
                 'conventions'}

MONC_ID_ATTR = 'MONC timestep'  # It is essential that this be preserved across all modified datasets, either as attribute or variable.

OPTIONS_DATABASE = {
    'variable': 'options_database',
    'dimensions': {'number_options', 'kvp', 'string'}
}

FROM_INPUT_FILE = {
    'reftime': 'time'
}

CHUNKING_DIMS = {'z': 2, 'zn':2}  # Ideally chunk sizes should be power of 2.

DEFAULT_CALENDAR = 'proleptic_gregorian'  # As per ISO 8601:2004.

DEPHY_OPTIONS = ('dephy_file',)  # 'dephy_forcings_enabled'
DROP_FOR_DEPHY = ("longitude", "latitude", "z0")

# Global attributes and their Numpy data types.
do_not_propagate = {'MONC timestep': np.int32}
split_attrs = {
    'MONC time': (
        lambda i, ds_list: ds_list[i]['time'].data.tolist(), 
        np.float64
    ),
    'Previous diagnostic write at': (
        lambda i, ds_list: ds_list[0].attrs['Previous diagnostic write at'] 
        if i == 0 else ds_list[i - 1]['time'].data.tolist(), 
        np.float64
    )
}
group_attrs = {'created': datetime}

'''
Define required fields in vocabulary file/dictionary, as key: value pairs.
If a set given as value, this contains all valid values.
If an empty dictionary is given as a value, valid values for that field are
themselves key: value pairs, where the key is the existing field value and 
the value is the new field value, e.g. 'dimension_changes': {'x': 'xu'} means 
that, for the specified variable, dimension x will be replaced by dimension 
xu.
If an empty list is given as a value, the expected value is a list.
If perturbation_to_absolute == True, reference_value must be the name of the 
variable that is added to make the conversion.
'''
VOCAB_FIELDS = {
    'dims': [],
    'updated_name': None, 
    'units': None,
    'axis': None,
    'standard_name': None,
    'long_name': None,
    'dimension_changes': {},
    'perturbation_to_absolute': {True, False},
    'reference_variable': None
}

# Define dimensions by which files & their contained variables are to be 
# divided. Dimensions exclude time, i.e. 0d variables vary only with time.
VOCAB_DIMS = {'0d', '1d', '2d', '3d'}

app_dir = os.getcwd()
if 'dev' not in app_dir:
    if 'test_data' not in app_dir:
        app_dir = os.path.join(app_dir, 'dev')
    else:
        app_dir = os.path.join(os.path.dirname(app_dir), 'dev')

config_file = open(os.path.join(app_dir, "config.yml"))
CONFIG = yaml.safe_load(config_file)
config_file.close()
for k, v in CONFIG.items():
    if isinstance(v, str) and v.lower() == 'none':
        CONFIG[k] = None

vocab_file = open(os.path.join(app_dir, "vocabulary.yml"))
vocabulary = yaml.safe_load(vocab_file)
vocab_file.close()
reference_vars = {}

# If dimensions (outermost key) aren't simply digits, assume the dimension is
# the first character:
if all([isinstance(k, str) for k in vocabulary.keys()]):
    vocabulary = {int(k[0]): v for k, v in vocabulary.items()}

# Validate entries in vocabulary:
for dim, group in vocabulary.items():
    for variable, attributes in group.items():
        to_drop = []
        for k, v in attributes.items():
            if k not in VOCAB_FIELDS:
                print(f'WARNING: {os.path.join(app_dir, "vocabulary.yml")} contains invalid field: {k}')
                to_drop.append(k)
                continue
                # raise KeyError(f'Vocabulary file contains invalid field: {k}')
            if VOCAB_FIELDS[k] is not None:
                if isinstance(VOCAB_FIELDS[k], set):
                    if v not in VOCAB_FIELDS[k]:
                        raise ValueError(
                            f'setup: Variable {variable}: {v} is not a valid '
                            f'value for {k}. Valid values: {VOCAB_FIELDS[k]}'
                        )
                if (isinstance(VOCAB_FIELDS[k], dict) and 
                    not isinstance(v, dict)):
                    raise TypeError(
                        f'setup: Variable {variable}: {k} requires a key:value '
                        'pair.'
                    )
                if isinstance(VOCAB_FIELDS[k], (list, tuple)) and not isinstance(v, (list, tuple)):
                    raise TypeError(
                        f'setup: Variable {variable}: {k} requires either a '
                        'list or tuple.'
                    )
            
        # Remove any invalid fields from variable.
        [attributes.pop(k) for k in to_drop]        
        
        # If any variables are perturbations that the user wants to convert to
        # an absolute quantity, record name of reference variable.
        if 'perturbation_to_absolute' in attributes and attributes['perturbation_to_absolute']:
            if 'reference_variable' not in attributes or not attributes['reference_variable']:
                raise KeyError(f'setup: {variable}: If '
                               'perturbation_to_absolute is True, '
                               'reference_variable must contain the name of '
                               'the variable containing reference value(s).')
            reference_vars[attributes['reference_variable']] = {
                'for': (dim, variable)
            }

# If any variables are perturbations that the user wants to convert to
# an absolute quantity, by adding a reference variable, track where to
# find each reference variable by dimension.
for var in reference_vars.keys():
    # Locate reference variable by dimension group
    for dim, dim_vars in vocabulary.items():
        if var in dim_vars.keys():
            reference_vars[var]['dim'] = dim
            vocabulary[reference_vars[var]['for'][0]][reference_vars[var]['for'][1]].update({'ref_dim': dim, 'ref_array': None})
    # Throw error if ref variable not found
    if not reference_vars[var]:
        raise ValueError(
            f'setup: Reference variable {var} not found in vocabulary.'
            )

# TODO: look for standard_name table, and if not found, download.
# Refer to cf-checker.
