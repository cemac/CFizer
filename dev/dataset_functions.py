import xarray as xr
from by_class import TimeUnits
from setup import *
from by_class import DsGroup
from utils import stem_str, type_from_str


def is_monc(dataset: xr.Dataset) -> bool:
    return MONC_ID_ATTR in set(dataset.attrs).union(set(dataset.variables))


def get_n_dims(dataset: xr.Dataset) -> int:
    return max([len(dataset[v].dims) for v in dataset.variables if v != OPTIONS_DATABASE['variable']])


def add_cf_attrs(dataset: xr.Dataset) -> xr.Dataset:
    pass


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


def merge_ds_groups(*args) -> DsGroup:
    # check each argument is of DsGroup type
    if not all(isinstance(arg, DsGroup) for arg in args):
        raise TypeError('Expecting a sequence of DsGroup objects.')
    if not all([isinstance(group.processed, xr.Dataset) for group in args]):
        raise TypeError("Each DsGroup.processed attribute must contain a single dataset.")
    
    # to_merge = [arg.processed[0] for arg in args]  # grab datasets
    # Find matching time-points

    # Check each DsGroup contains a single dataset in its processed attribute
    if len(args) == 1:
        return args[0].processed

    if len(args) == 2:
        # Find common stem
        stem = stem_str(args[0].stem, args[1].stem)

        # Find common name, if names don't already match
        name = args[0].name if args[0].name == args[1].name else stem_str(
            args[0].name, args[1].name
        )

        # Merge lists of filepaths from each group
        filepaths = args[0].filepaths + args[1].filepaths

        # To preserve original datasets, create a copy of one.
        new_ds = args[0].processed.copy(deep=True)

        # Merge 2nd dataset into the copy of the first.
        new_ds.merge(
            other=args[1].processed, 
            join='exact', 
            combine_attrs='drop_conflicts'
        )

        new_group = DsGroup(name=name, action='merge_dims', filepaths=filepaths)
        new_group.process()

        return args[0].processed.merge(
            other=args[1].processed, 
            join='exact', 
            combine_attrs='drop_conflicts'
        )
    
    else:
        # call recursively until only 2 args to process
        return merge_ds_groups(args[0], merge_ds_groups(*args[1:]))
    


def merge_dimensions(*args) -> xr.Dataset:
    '''
    This is designed to merge datasets whose time coordinates match, but contain
    mutually exclusive sets of data variables.
    '''
    
    if not all(isinstance(arg, xr.Dataset) for arg in args):
        raise TypeError('All parameters must be of type xarray.Dataset')
    
    if len(args) == 1:
        return args[0].copy(deep=True)
    
    if len(args) == 2:
        if all(['created' in ds.variables for ds in args]):
            # if all(args[0].created.data > args[1].created.data):
            #     new_ds = args[0].copy(deep=True)
            #     to_merge = args[1]
            # elif all(args[1].created.data > args[0].created.data):
            #     new_ds = args[1].copy(deep=True)
            #     to_merge = args[0]
            if max(args[0].created.data) > max(args[1].created.data):
                new_ds = args[0].copy(deep=True)
                to_merge = args[1]
            else:
                new_ds = args[1].copy(deep=True)
                to_merge = args[0]
            # last_created = max([max(ds.created.data) for ds in args])
            new_ds.merge(
                other=to_merge, 
                join='exact', 
                combine_attrs='drop_conflicts',
                overwrite_vars='created'
            )  # overwrite_vars keeps base ds' variable values if conflict.
            # combine_attrs only seems to apply to attributes of variables, 
            # etc, not global attributes.
        else:
            new_ds = args[0].copy(deep=True)
            # to_merge = args[1]
            try:
                new_ds.merge(
                    other=args[1], 
                    join='exact', 
                    combine_attrs='drop_conflicts'
                )
            except xr.MergeError as e:
                raise e
        return new_ds
    
    # call recursively until only 2 args to process
    return merge_dimensions(args[0], merge_dimensions(*args[1:]))
    