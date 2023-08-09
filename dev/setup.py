'''
TODO: This will ultimately move into __init__.py in app's base directory.
TODO: throw verbose errors if config file not found; if vocabulary file or dict not found, look for Excel spreadsheet and pass to vocab_from_xls if found; if neither found, throw error.
'''

import yaml
import os
from utils import vocab_from_xls


MONC_ID_ATTR = 'MONC timestep'

OPTIONS_DATABASE = {
    'variable': 'options_database',
    'dimensions': {'number_options', 'kvp'},
    'dx': 'dxx',
    'dy': 'dyy'
}

INPUT_FILE = {
    'reftime': 'time'
}

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

# Define required fields in vocabulary file/dictionary, as key: value pairs.
# If tuple given as value, this contains all valid values.
# If an empty dictionary is given as a value, valid values for that field are
# themselves key: value pairs, where the key is the existing field value and 
# the value is the new field value, e.g. 'dimension_changes': {'x': 'xu'} means 
# that, for the specified variable, dimension x will be replaced by dimension 
# xu.
VOCAB_FIELDS = {
    'updated_name': None, 
    'units': None,
    'standard_name': None,
    'long_name': None,
    'dimension_changes': {},
    'reftime': {'options_database', 'input_file'}
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

vocab_file = open(os.path.join(app_dir, "vocabulary.yml"))
VOCABULARY = yaml.safe_load(vocab_file)
vocab_file.close()
