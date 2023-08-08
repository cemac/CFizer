import xarray as xr
import re
from datetime import datetime as dt
import cfunits
import os
from utils import is_sequence, is_sequence_of


def merge_ds(ds1: xr.Dataset, ds2: xr.Dataset) -> xr.Dataset:
    pass
    

def coord_from_dim(
        ds: xr.Dataset, 
        dim: str, 
        data: list|tuple|set, 
        attrs: dict = None
) -> xr.Dataset:
    '''
    '''
    # TODO: ensure all arguments are valid
    if not is_sequence(data):
        raise TypeError

    if dim in ds.coords:
        # Add any missing information, or update as required
        if attrs is not None:
            for a, v in attrs.items():
                ds[dim].attrs[a] = v
        if not all(ds[dim].data == data):
            attrs = ds[dim].attrs
            try:
                ds.update({dim: (dim, data)})
            except ValueError as e:
                raise e
            for a, v in attrs.items():
                ds[dim].attrs[a] = v
            
    else:
        new_var = xr.Variable(dims=dim, attrs=attrs, data=data)
        ds = ds.assign(variables={dim: new_var})
    return ds
    

def coord_points(n: int, delta: float|int, first: float|int = 0) -> list:
    # TODO: ensure all arguments are valid

    return [(x * delta + first) for x in range(n)]


def change_dims(
        ds: xr.Dataset, data_var: str, dim_changes: dict
) -> xr.Dataset:
    '''
    Change the dimensions of a data variable.
    ds: xarray.Dataset
    dim_changes: dictionary mapping existing dimension names to new dimension names (e.g. {'x': 'xu', 'y': 'yu'})
    Only specified dimensions are changed; others are left as is.
    '''
    # TODO: Check dim_changes is valid type, and check each key
    # is an existing dim.

    ds[data_var] = ds[data_var].swap_dims(dim_changes)  # (new_dims, ds[data_var].data, ds[data_var].attrs)
    return ds


def replace_dim(
        ds: xr.Dataset, data_var: str, current_dim: str, new_dim: str
) -> xr.Dataset:
    '''
    Change one specific dimension of a data variable.
    '''
    # Check parameters are valid types, and check current_dim &
    # new_dim are existing dims, and that current_dim is a dim of data_var.
    if not isinstance(ds, xr.Dataset):
        raise TypeError('ds parameter must be an xarray Dataset.')
    for parameter in [data_var, current_dim, new_dim]:
        if not isinstance(parameter, str):
            raise TypeError(
                f'data_var, current_dim and new_dim must all be strings.'
            )
    if data_var not in ds.data_vars:
        raise KeyError(f'{data_var} is not a data variable of dataset.')
    for d in [current_dim, new_dim]:
        if not d in ds.dims:
            raise KeyError(f'{d} is not a dimension in dataset.')
    
    # Re-assign data variable with new dimensions & all existing data.
    ds[data_var] = ds[data_var].swap_dims({current_dim: new_dim})  # (data_var_dims, ds[data_var].data, ds[data_var].attrs)
    return ds


class Options_Db:
    def __init__(self, ds: xr.Dataset):
        '''
        Currently, this hard-codes the expectation that the options
        database comprises a list (key, value) pairs, stored in the
        options_database data variable.
        This layout could be instead passed in, with parameters taken
        from the config file.
        '''
        # Get key: value pairs from options database variable, converting
        # binary data to strings.
        self.as_dict = {k.decode('utf-8'): v.decode('utf-8') 
                    for [k, v] in ds['options_database'].data}
        
        # Validation: check that number_options dimension size agrees
        # with len(self.as_dict).
        if ds.dims['number_options'] != len(self.as_dict):
            raise ValueError(
                'Mismatch between number of key:value pairs and number_options dimension size.'
            )

        # Convert number-containing strings to numbers and Boolean-containing
        # strings to Booleans.
        for k, v in self.as_dict.items():
            if v.lower() == 'true':
                self.as_dict[k] = True
            elif v.lower() == 'false':
                self.as_dict[k] = False
            else:
                try:
                    float(v)
                    # Note: string.isnumeric() returns false for numbers containing decimal points or scientific notation.
                except ValueError:
                    # not a number
                    pass  # leave as a string
                else:
                    if '.' not in v:
                        try:
                            int(v)
                        except ValueError:
                            # not an integer;
                            # should never reach this option.
                            self.as_dict[k] = float(v)
                        else:
                            self.as_dict[k] = int(v)
                    else:
                        self.as_dict[k] = float(v)  # floating point, whether integer or not

        # For convenience, create attributes from each key.
        # Alternatively, do this only for keys requested in
        # config file.
        for k, v in self.as_dict.items():
            self.__setattr__(k, v)


def test():
    test_dir = '/home/earcwi/OneDrive/EUREC4A/Code/CFizer/test_data'
    nc3d = xr.open_dataset(
        os.path.join(test_dir, 'd20200128_diagnostic_3d_313200.nc')
    )
    odb = Options_Db(nc3d)
    print(odb.dxx)


if __name__ == '__main__': test()