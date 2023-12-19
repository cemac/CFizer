import numpy as np
import sys
import xarray
from typing import Union


def str_to_class(module: str, class_name: str):
    """
    Returns a class with the name provided as a string, which can be used for type-casting.
    
    module: name of the module defining the class, e.g. numpy
    class_name: the name of the wanted class, e.g. float64

    Adapted from https://stackoverflow.com/a/1176180
    """
    return getattr(sys.modules[module], class_name)


def first_rest(
    first_var: str,
    rest_var: str,
    module_type: str,
    i: int,
    ds_list: list[xarray.Dataset],
):
    """
    Returns classes for the variables to be used for the first and remaining
    values for a global variable, in a dataset to be split.

    first_var:  name of variable providing first dataset in series' value.
    rest_var: name of variable providing values for remaining datasets in series.
    module_type: the data type to be applied, including the module name, e.g. 'numpy.float64'
    i: dataset's position in series - 0 is first.
    ds_list: list containing datasets resulting from split.
    """
    var = first_var if i == 0 else rest_var
    if var in ds_list[0].variables:
        return str_to_class(*module_type.split("."))(ds_list[0][var].data.tolist())
    elif var in ds_list[0].attrs:
        return str_to_class(*module_type.split("."))(ds_list[0].attrs[var])
    else:
        raise KeyError("Variable(s) not found in dataset's variables or attributes.")


def type_from_str(string: str):
    """
    Uses numpy dtypes, to ensure xarray.Dataset.to_netcdf gives predictable
    results.
    Using double precision for floats is wasteful, but using single precision
    results in numbers being off due to precision, and their displaying with a
    trailing 'f' in NC files.
    """
    if not isinstance(string, str):
        # Can't work with it, so return as is.
        return np.str_(string)
    # Because xarray doesn't accept bool or np.bool_ values for variables/
    # attributes, leave any of these as strings.
    if string.lower() == "true":
        return np.str_(string.lower())  # return True
    elif string.lower() == "false":
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
            if "." not in string:
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


def generate_coords(
    number: int,
    spacing: Union[float, int],
    midpoint: bool = False,
    data_type: np.dtype = np.float64,
) -> np.ndarray:
    """
    Returns a numpy array containing points for a new coordinate variable.
    """
    first = 0.5 if midpoint else 0
    return np.arange(first, number, 1, dtype=data_type) * spacing


def dict_strings(d: dict) -> dict:
    """
    Converts a potentially complex (nested) dictionary to a collection 
    (dictionary or list) of strings.
    Will also handle lists, tuples and sets. Datasets and Data Arrays will
    return None (to skip any such stored e.g. in vocabulary), while it attempts
    to return the string representation of any other types using a __str__ 
    method, if defined.
    This is used in generating the log file.
    """
    all_strings = {}
    if not isinstance(d, dict):
        if isinstance(d, (list, tuple, set)):
            return [str(element) for element in d]
        elif isinstance(d, (xarray.Dataset, xarray.DataArray)):
            return None
        else:
            return str(d)
    for k, v in d.items():
        if v is None:
            continue
        elif isinstance(v, dict):
            all_strings[str(k)] = dict_strings(v)
        elif isinstance(v, (list, tuple, set)):
            all_strings[str(k)] = []
            for element in v:
                if isinstance(element, (list, tuple, set)):
                    all_strings[str(k)].append([str(x) for x in element])
                elif isinstance(element, dict):
                    all_strings[str(k)].append(dict_strings(element))
                else:
                    all_strings[str(k)].append(str(element))
        elif isinstance(v, (xarray.Dataset, xarray.DataArray)):
            continue
        else:
            all_strings[str(k)] = str(v)
    return all_strings
