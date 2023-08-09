import xarray as xr
from by_class import TimeUnits
from setup import *


def is_monc(dataset: xr.Dataset) -> bool:
    pass


def get_n_dims(dataset: xr.Dataset) -> int:
    return max([len(dataset[v].dims) for v in dataset.variables if v != OPTIONS_DATABASE['variable']])


def add_cf_attrs(dataset: xr.Dataset) -> xr.Dataset:
    pass


def attr_to_var(dataset: xr.Dataset, attr: str = None) -> xr.Dataset:
    if attr is None:
        # convert all attributes to variables, using this function recursively.
        for attr in dataset.attrs.keys():
            attr_to_var(dataset=dataset, attr=attr)
    # TODO


def missing_coords(dataset: xr.Dataset) -> xr.Dataset:
    pass


def cfize_variables(dataset: xr.Dataset) -> xr.Dataset:
    # run cfize_variable on each variable, including coordinate variables.
    pass
        

def time_units_from_input(dataset: xr.Dataset) -> TimeUnits:
    pass


def save_cf_dataset(dataset: xr.Dataset, filepath: str, compression=None):
    # This calls xr.Dataset.to_netcdf() with any required compression options.
    pass


def merge_ds_groups(*args) -> xr.Dataset:
    '''
    This is designed to merge datasets whose time coordinates match, but contain
    mutually exclusive sets of data variables.
    '''
    # TODO: check each argument is of DsGroup type & processed attribute of each contains one dataset.

    # Check each DsGroup contains a single dataset in its processed attribute
    if len(args) == 1:
        return args[0].processed[0]

    pass