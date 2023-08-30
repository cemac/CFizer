from time import perf_counter
import os
import numpy as np


def performance_time(func):
	def wrapper(*args, **kwargs):
		start_time = perf_counter()
		response = func(*args, **kwargs)  # run the wrapped function
		end_time = perf_counter()
		duration = end_time - start_time
		print(f"Process {os.getpid()}: {func} took {duration} seconds.")
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


def generate_coords(number: int, 
                    spacing: float|int, 
                    midpoint: bool = False,
                    data_type: np.dtype = np.float64) -> np.ndarray:
    first = 0.5 if midpoint else 0
    return np.arange(first, number, 1, dtype=data_type) * spacing


def tf(string) -> bool:
    '''
    Adapted from https://stackoverflow.com/a/43357954
    '''
    if isinstance(string, bool):
        return string
    return True if string.lower() in {
        'true', 'yes', 't', 'y', '1'
    } else False
