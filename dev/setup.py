'''
TODO: This will ultimately move into __init__.py in app's base directory.
TODO: throw verbose errors if config file not found; if vocabulary file or dict not found, look for Excel spreadsheet and pass to vocab_from_xls if found; if neither found, throw error.
Name:           cfizer.py
Author:         Cameron Wilson, CEMAC, University of Leeds
Date:           August 2023
File Revision:  ...
CFizer Version: See __version__

'''

import yaml
import os
import sys
from cfizer import VOCAB_FIELDS


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

INPUT_FILE = {
    'reftime': 'time'
}

DEFAULT_CALENDAR = 'proleptic_gregorian'  # As per ISO 8601:2004.

DEPHY_OPTIONS = ('dephy_file',)  # 'dephy_forcings_enabled'
DROP_FOR_DEPHY = ("longitude", "latitude", "z0")

# # Categorise files, based on number of dimensions (in addition to time, 
# # i.e. time only is 0d), & how they are to be processed.
# # DIM_GROUPS = ['0+1d', '0+1d', '2d', '3d']
# DIM_GROUPS = {0: {'action': 'merge',
#                   'group': '0+1d'}, 
#               1: {'action': 'merge',
#                   'group': '0+1d'}, 
#               2: {'action': 'merge',
#                   'group': '2d'}, 
#               3: {'action': 'split',
#                   'group': '3d'}}

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

try:
    vocab_file = open(os.path.join(app_dir, "vocabulary.yml"))
except:
    sys.exit(
        f"Vocabulary file not found in application directory: "
        f"{os.path.join(app_dir, 'vocabulary.yml')}. If Excel spreadsheet "
        f"version if available, run the following command from the application "
        f"directory (app_dir):\n ./vocab_from_xlsx.py <path_to_Excel_file>"
    )
else:
    VOCABULARY = yaml.safe_load(vocab_file)
    vocab_file.close()
    for variable, attributes in VOCABULARY.items():
        # if 'dims' not in attributes:
        #     raise KeyError('All variables must be specified along with their dimensions.')
        for k, v in attributes.items():
            if k not in VOCAB_FIELDS:
                raise KeyError(f'Vocabulary file contains invalid field: {k}')
            if VOCAB_FIELDS[k] is not None:
                if isinstance(VOCAB_FIELDS[k], tuple):
                    if v not in VOCAB_FIELDS[k]:
                        raise ValueError(f'Variable {variable}: {v} is not a valid value for {k}. Valid values: {VOCAB_FIELDS[k]}')
                if isinstance(VOCAB_FIELDS[k], dict) and not isinstance(v, dict):
                    raise TypeError(f'Variable {variable}: {k} requires a key:value pair.')
                if isinstance(VOCAB_FIELDS[k], (list, tuple)) and not isinstance(v, (list, tuple)):
                    raise TypeError(f'Variable {variable}: {k} requires either a list or tuple.')
                    
        if 'perturbation_to_absolute' in attributes and attributes['perturbation_to_absolute']:
            if 'reference_variable' not in attributes or not attributes['reference_variable']:
                raise KeyError(f'{variable}: If perturbation_to_absolute is True, reference_variable must contain the name of the variable containing reference value(s).')

# TODO: look for standard_name table, and if not found, download.
# Refer to cf-checker.
