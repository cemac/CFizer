# import re
from time import perf_counter
# from datetime import timedelta, datetime
# from cfunits import Units
import os
# from difflib import SequenceMatcher
import numpy as np


def performance_time(func):
	def wrapper(*args, **kwargs):
		start_time = perf_counter()
		response = func(*args, **kwargs)  # run the wrapped function
		end_time = perf_counter()
		duration = end_time - start_time
		print(f"{func} on process {os.getpid()} took {duration} seconds.")  # ({args}, {kwargs})
		return response
	return wrapper


def vocab_from_xls(filepath: str) -> dict:
    '''
    This will for Excel file to convert to vocabulary dictionary, if it 
    contains the right columns (fields).
    It will attempt to fill any missing data based on CF conventions, including:
    - inferring cell_method from presence of terms such as "mean" and "max" in variable name.
    - looking up standard_name, if supplied.
    It will check any units given using `cfunits`, and check against CF standard names definitions.

    The intention is to compile a full lookup table, such that the main 
    application does not need to infer anything or look anywhere else for the 
    required data.
    '''
    import openpyxl as xl

    VOCABULARY = {}
    # 
    return VOCABULARY


# def options_db_to_global_attrs():
#     pass


# def is_sequence(object):
#     return isinstance(object, (
#         list,
#         set,
#         tuple,
#         range
#     ))


# def is_sequence_of(object, type: type):
#     if is_sequence(object):
#         return all([isinstance(element, type) for element in object])
#     else:
#         return False


def type_from_str(string: str):
    '''
    Uses numpy dtypes, to ensure xarray.Dataset.to_netcdf gives predictable 
    results.
    Using double precision for floats is wasteful, but using single precision
    results in numbers being off due to precision, and their displaying with a 
    trailing 'f' in NC files.
    '''
    if not isinstance(string, str):
        # Can't work with it, so return as is.
        return np.str_(string)
    # Because xarray doesn't accept bool or np.bool_ values for variables/
    # attributes, leave any of these as strings.
    if string.lower() == 'true':
        return np.str_(string.lower())  # return True
    elif string.lower() == 'false':
        return np.str_(string.lower())  # return False
    else:
        try:
            float(string)
            # Note: string.isnumeric() returns false for numbers containing 
            # decimal points or scientific notation.
        except ValueError:
            # not a number
            return np.str_(string)  # leave as a string
        else:
            if '.' not in string:
                try:
                    int(string)
                except ValueError:
                    # not an integer;
                    # should never reach this option.
                    return np.double(string)
                else:
                    return np.int32(string)
            else:
                return np.double(string)  # floating point, whether integer or not


# def decode_time_units(units: str) -> tuple:
#     try:
#         (unit, ref_datetime) = units.lower().split(' since ')
#     except AttributeError or ValueError:
#         raise ValueError('''Argument must be a string in the format:
#                         <time-unit> since <time-origin>,
#                          with <time-origin> in ISO format.''')
    
#     try:
#         dt = datetime.fromisoformat(ref_datetime)
#     except:
#         raise ValueError('Origin date-time must be in ISO format.')
    
#     return (unit, dt)


# def timedelta_from_units(number: float|int, unit: str) -> timedelta:
    
#     try:
#         float(number)
#     except ValueError:
#         raise TypeError('First argument must be numeric')

#     time_units = dict(
#         days=0,
#         seconds=0, 
#         microseconds=0, 
#         milliseconds=0, 
#         minutes=0, 
#         hours=0, 
#         weeks=0
#     )

#     if unit not in time_units:
#         raise KeyError('Invalid unit supplied')
    
#     time_units[unit] = number

#     return timedelta(
#         days=time_units['days'],
#         seconds=time_units['seconds'],
#         microseconds=time_units['microseconds'],
#         milliseconds=time_units['milliseconds'],
#         minutes=time_units['minutes'],
#         hours=time_units['hours'],
#         weeks=time_units['weeks']
#     )

# def datetime_from_relative_units(number: float|int, units: str) -> datetime:
    
#     try:
#         (unit, ref_datetime) = decode_time_units(units)
#     except TypeError:
#         raise TypeError('''Second argument must be a string in the format:
#                         [time unit] since [date or date-time in ISO format]''')

#     delta = timedelta_from_units(number, unit)

#     return ref_datetime + delta


def generate_coords(number: int, 
                    spacing: float|int, 
                    midpoint: bool = False,
                    data_type: np.dtype = np.float64) -> np.ndarray:
    first = 0.5 if midpoint else 0
    return np.arange(first, number, 1, dtype=data_type) * spacing


# def stem_str(*args):
#         if not args:
#             raise ValueError('At least one string must be provided.')
#         if None in args:
#             return stem_str(*[arg for arg in args if arg is not None])
#         filenames = [os.path.splitext(os.path.basename(path))[0] for path in args]  # If each path contains only a filename (no extension) this will still work.
#         stem = filenames[0]
#         if len(args) > 1:
#             for f in filenames[1:]:
#                 match = SequenceMatcher(a=stem, b=f).find_longest_match()
#                 stem = stem[match.a: match.a + match.size]
#         return stem
    

def tf(string) -> bool:
    '''
    Adapted from https://stackoverflow.com/a/43357954
    '''
    if isinstance(string, bool):
        return string
    return True if string.lower() in {
        'true', 'yes', 't', 'y', '1'
    } else False


# if __name__ == '__main__':
    # value = 315960.
    # units = 'seconds since 2020-01-25'
    # decoded_units = decode_time_units(units)
    # print('decode_time_units:', decoded_units)
    # print('timedelta_from_units:', timedelta_from_units(value, decoded_units[0]))
    # print('datetime_from_relative_units:', datetime_from_relative_units(value, units))
    # print(f'stem_str("spam", None): {stem_str("spam", None)}')
    # print('kg/m^2 =', format_units('kg/m^2'))
    # print('kg/m^2/s =', format_units('kg/m^2/s'))
    # print('kg m/s^2 =', format_units('kg m/s^2'))
    # print('kg/kg =', format_units('kg/kg'))
