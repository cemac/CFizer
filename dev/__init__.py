import os
import yaml
import sys
from utils import str_to_class, first_rest

VERSION = '0.1.1'
CF_VERSION = '1.10'

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

# WILDCARDS = {'*', '?'}
verbose = False  # Default
quiet = False  # Default
CF_ATTRIBUTES = {'title', 
                 'institution', 
                 'source', 
                 'history', 
                 'references', 
                 'comment',
                 'Conventions'}
DEFAULT_CALENDAR = 'proleptic_gregorian'  # As per ISO 8601:2004.

app_dir = os.getcwd()
"""TODO: UPDATE FOR PRODUCTION VERSION'S NEW DIRECTORY STRUCTURE"""
if 'dev' not in app_dir:
    if 'test_data' not in app_dir:
        app_dir = os.path.join(app_dir, 'dev')
    else:
        app_dir = os.path.join(os.path.dirname(app_dir), 'dev')

try:
    config_file = open(os.path.join(app_dir, "config.yml"))
except:
    raise OSError(
        f"config.yml not found in application directory, {app_dir}"
    )
CONFIG = yaml.safe_load(config_file)
config_file.close()
for k, v in CONFIG.items():
    if isinstance(v, str) and v.lower() == 'none':
        CONFIG[k] = None

'''
Check compression details.
Only lossless compression is considered here; all precision is maintained.
WARNING: only tested with zlib compression.
'''
if CONFIG['compression']['on']:
    if CONFIG['compression']['level'] < 1 or CONFIG['compression']['level'] > 9:
        raise ValueError("Only compression levels from 1 to 9 are valid.")
    if CONFIG['compression']['type'] not in {
        'zlib', 'szip', 'zstd', 'bzip2', 'blosc_lz', 'blosc_lz4', 
        'blosc_lz4hc', 'blosc_zlib', 'blosc_zstd'
    }:
        raise ValueError("Invalid compression algorithm specified.")
    
DIM_GROUPS = CONFIG['dimension_groups']
DIM_ACTIONS = CONFIG['group_actions']
MONC_ID_ATTR = CONFIG['monc_id_attribute']  # It is essential that this be preserved across all modified datasets, either as attribute or variable.
OPTIONS_DATABASE = CONFIG['options_database']
FROM_INPUT_FILE = {
    'reftime': 'time'
} if CONFIG['reference_time'] == 'input_file' else None
default_time_unit = CONFIG['time_units'] if 'time_units' in CONFIG and CONFIG['time_units'] else 's'
CHUNKING_DIMS = CONFIG['chunking']
DEPHY_OPTIONS = CONFIG['dephy_true_if']
DROP_FOR_DEPHY = CONFIG['drop_for_dephy']
do_not_propagate = {k: str_to_class(*v.split('.')) for k, v in CONFIG['do_not_propagate'].items()}  #'MONC timestep': np.int32
split_attrs = {}
for attr, details in CONFIG['attributes_to_split'].items():
    split_attrs[attr] = (
        lambda i, ds_list: first_rest(first_var=details[0][0], 
                                      rest_var=details[0][1], 
                                      module_type=details[1], 
                                      i=i, 
                                      ds_list=ds_list),
        str_to_class(*details[1].split('.'))
    )
group_attrs = {k: str_to_class('builtins', v) 
               if '.' not in v else str_to_class(*v.split('.'))
               for k, v in CONFIG['group_attributes'].items()}

try:
    vocab_file = open(os.path.join(app_dir, "vocabulary.yml"))
except:
    exit(
        f"Vocabulary file ({os.path.join(app_dir, 'vocabulary.yml')}) not "
        f"found. To derive from Excel spreadsheet, run the following command "
        f"from the application directory:"
        f"{os.path.join(app_dir, 'vocab_from_xlsx.py')}"
    )
vocabulary = yaml.safe_load(vocab_file)
vocab_file.close()
print(f"Vocabulary loaded from {os.path.join(app_dir, 'vocabulary.yml')}.")
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
            if k == 'units':
                vocabulary[dim][variable][k] = str(v)  # This ensures cfunits.Units.isvalid works correctly on fractions (units = 1).
            if VOCAB_FIELDS[k] is not None:
                if isinstance(VOCAB_FIELDS[k], set):
                    if type(list(VOCAB_FIELDS[k])[0])(v) not in VOCAB_FIELDS[k]:
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

print("Vocabulary validated.")


class ConfigError(Exception):
    "Raised when something critical is missing from, or incorrectly defined in, config.yml. This should cause the program to exit immediately."
    def __init__(self, *args: object) -> None:
        super().__init__(*args)


class ConfigWarning(Exception):
    "Raised when something non-critical is missing from, or incorrectly defined in, config.yml. When this arises, it is intended that the program will be allowed to continue, but ending with a non-zero exit code."
    def __init__(self, *args: object) -> None:
        super().__init__(*args)


class VocabError(Exception):
    "Raised when something critical is missing from, or incorrectly defined in, vocabulary.yml. This should cause the program to exit immediately."
    def __init__(self, *args: object) -> None:
        super().__init__(*args)
    

class VocabWarning(Exception):
    "Raised when something non-critical is missing from, or incorrectly defined in, vocabulary.yml. When this arises, it is intended that the program will be allowed to continue, but ending with a non-zero exit code."
    def __init__(self, *args: object) -> None:
        super().__init__(*args)


