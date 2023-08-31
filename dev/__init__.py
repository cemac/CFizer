__version__ = '1.0.0'

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
