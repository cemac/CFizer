# from time import perf_counter
# import os
import numpy as np
import sys
import xarray
import datetime, time


def str_to_class(module: str, class_name: str):
    '''
    Adapted from https://stackoverflow.com/a/1176180
    '''
    return getattr(sys.modules[module], class_name)


def first_rest(first_var: str, rest_var: str, module_type: str, i: int, ds_list:list[xarray.Dataset]):
    # print("first_rest received:", first_var, rest_var, np_type, i, ds_list)
    var = first_var if i == 0 else rest_var
    if var in ds_list[0].variables:
        return str_to_class(*module_type.split('.'))(ds_list[0][var].data.tolist())
    elif var in ds_list[0].attrs:
        return str_to_class(*module_type.split('.'))(ds_list[0].attrs[var])
    else:
        raise KeyError(
            "Variable(s) not found in dataset's variables or attributes."
        )
# def performance_time(func):
#     def wrapper(*args, **kwargs):
#         start_time = perf_counter()
#         response = func(*args, **kwargs)  # run the wrapped function
#         end_time = perf_counter()
#         duration = end_time - start_time
#         if ('shared' in kwargs and 
#             'verbose' in kwargs['shared'] and 
#             kwargs['shared']['verbose']):
#             print(f"Process {os.getpid()}: {func} took {duration} seconds.")
#         return response
#     return wrapper


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


# def tf(string) -> bool:
#     '''
#     Adapted from https://stackoverflow.com/a/43357954
#     '''
#     if isinstance(string, bool):
#         return string
#     return True if string.lower() in {
#         'true', 'yes', 't', 'y', '1'
#     } else False
