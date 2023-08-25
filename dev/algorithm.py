#!/usr/bin/env python

# Initialise, incl. from VOCAB & CONFIG files.
from setup import *
import xarray as xr
import os
import os.path as op
from glob import iglob
from cfunits import Units
import numpy as np
from datetime import datetime
import typing
from difflib import SequenceMatcher
from utils import type_from_str
from dask.delayed import Delayed
from dask.diagnostics import ProgressBar
from time import perf_counter, sleep
import argparse
import re
from multiprocessing import Process, Pool  # multiprocess not available in Jaspy environment.


def performance_time(func):
	def wrapper(*args, **kwargs):
		start_time = perf_counter()
		response = func(*args, **kwargs)  # run the wrapped function
		end_time = perf_counter()
		duration = end_time - start_time
		print(f"{func} on process {os.getpid()} took {duration} seconds.")  # ({args}, {kwargs})
		return response
	return wrapper


def ds_feed(file_list: list|set|tuple|str):
    '''
    This generator replaces the DsGroup.datasets generator, because generators 
    can't be passed (as a "pickle") to a new Process as an attribute, method or 
    function within an object.
    '''
    # ENHANCEMENT: test that all paths supplied are valid & are NC files.
    if isinstance(file_list, str):
        file_list = [file_list]
    return (xr.open_dataset(path, decode_times=False) for path in file_list)
    

def ds_list(file_list: list|set|tuple|str):
    # ENHANCEMENT: test that all paths supplied are valid & are NC files.
    if isinstance(file_list, str):
        file_list = [file_list]
    return [xr.open_dataset(path, decode_times=False) for path in file_list]


def variable_glob(ds: xr.Dataset, var_glob: str) -> list:
    '''
    Takes a glob-style string and finds matching variable(s) in a dataset.
    '''
    match_str = re.compile(var_glob.replace('*', '.*').replace('?','.'))
    return [match_str.fullmatch(v).string for v in set(ds.variables) if match_str.fullmatch(v)]


def stem_str(*args: str):
    '''
    Finds common segment of a list of filepaths/names. Also works on any 
    strings.
    '''
    if not args:
        raise ValueError('stem_str: At least one string must be provided.')
    if None in args:
        return stem_str(*[arg for arg in args if arg is not None])
    
    filenames = [op.splitext(op.basename(path))[0] for path in args]
    # If each path contains only a filename (no extension), or a string other 
    # than a filepath, this will still work.
    
    stem = filenames[0]
    if len(args) > 1:
        for f in filenames[1:]:
            match = SequenceMatcher(a=stem, b=f).find_longest_match()
            stem = stem[match.a: match.a + match.size]
    
    return stem


def get_time_var(*args: int|xr.Dataset) -> str:
    '''
    TODO: update to use VOCABULARY[dim]
    '''
    if isinstance(args[0], int):
        if args[0] not in VOCABULARY:
            raise KeyError(
                'get_time_var: Vocabulary file has no specifications for '
                f'{args[0]} spatial dimensions.')
        time_vars = [v for v in VOCABULARY[args[0]].keys() 
                     if 'time' in v.lower()]
        if not time_vars:
            raise KeyError('get_time_var: No time variable found in '
                           f'vocabulary for {args[0]}d datasets.')
    elif isinstance(args[0], xr.Dataset):
        time_vars = [d for d in args[0].dims if 'time' in d.lower()]
        if not time_vars:
            raise KeyError(f'No time variable found in dataset.')
    else:
        raise AttributeError(
            'get_time_var: Either a number of spatial dimensions or an xarray.'
            'Dataset must be passed in as a parameter.')
    if len(time_vars) > 1:
        as_string = ''
        for i, s in enumerate(time_vars):
            as_string += f'\n[{i}] {s}'
        while True:
            sel = input(f'Please select time variable for {args[0]}d datasets. Enter number: {as_string}')
            if isinstance(sel, int) and sel >= 0 and sel < len(time_vars):
                break
            else:
                print('Invalid selection.')
        return time_vars[sel]
    return time_vars[0]


def is_monc(dataset: xr.Dataset) -> bool:
    return MONC_ID_ATTR in set(dataset.attrs).union(set(dataset.variables))


def get_n_dims(dataset: xr.Dataset) -> int:
    return max([
        len(dataset[v].dims) for v in dataset.variables 
        if v != OPTIONS_DATABASE['variable']
        ])


def split_ds(dataset: xr.Dataset, var: str = 'time') -> list[xr.Dataset]:
    # Because xarray.Dataset.sel only creates datasets that point to the
    # original, any changes to global attributes etc in one dataset will
    # propagate to all. To create independent datasets, each resulting ds
    # must be deep-copied.

    # TODO: check groupby creates independent datasets. If not, need to copy.

    if var not in dataset.dims:
        raise AttributeError(f'split_ds: {var} not in dataset dimensions.')
    print(f'Splitting dataset {dataset.attrs["title"]} by {var}.')
    grouped = {v: ds for (v, ds) in dataset.groupby(var)}
    for k, v in grouped.items():
        v.attrs['title'] = v.attrs['title'].strip('_ +,.&') + '_' + str(int(k))
        print(f'Created new dataset with title, {v.attrs["title"]}')
    split = list(grouped.values())
    return split


@performance_time
def ds_to_nc(ds: xr.Dataset, filepath: str, encodings: dict = None) -> None:
    ds.to_netcdf(
        path=filepath, 
        encoding=encodings
        )
    

@performance_time
def ds_to_nc_dask(ds: xr.Dataset, filepath: str, encodings: dict = None) -> Delayed:
    return ds.to_netcdf(
        path=filepath, 
        encoding=encodings, 
        compute=False
        )
    # writer.compute()


@performance_time
def perform_write(writer: Delayed) -> None:
    writer.compute()


class DsGroup:
    '''
    Each file in a DsGroup has the same variables, coordinates & dimensions.
    '''
    def __init__(self, 
                 name: str = '', 
                 n_dims: int = None,
                 time_variable: str = '',
                 action: str = None, 
                 filepaths: list|set|tuple = None,
                 groups: list|set|tuple = None) -> None:
        '''
        n_dims: Number of *spatial* dimensions (X, Y, Z), ignoring variations 
                such as z & zn.
        '''
        # If groups parameter contains existing DsGroup objects:
        if groups and all([isinstance(g, DsGroup) for g in groups]):
            # check each DsGroup has only one dataset in its processed 
            # attribute.
            if not all([isinstance(g.processed, xr.Dataset) for g in groups]):
                raise TypeError('DsGroup: To merge DsGroup objects, the '
                                'processed attribute of each must contain a '
                                'single xarray.Dataset object.')
            self.name = stem_str(*[group.name for group in groups])
            self.n_dims = sum([g.n_dims for g in groups])  # This is a little dubious, given there's no guarantee dimension groups don't overlap, but works for MONC outputs where 0d & 1d are grouped.
            self.time_var = time_variable or groups[0].time_var
            if not all([g.time_var == self.time_var for g in groups]):
                raise ValueError('DsGroup: When attempting to combine groups, '
                                 'found mismatch in time variables.')
            self.action = action or 'merge_groups'
            if 'merge' in self.action:
                self.process = self.merge_groups
            self.stem = stem_str(*[g.stem for g in groups])
            self.filepaths = None
            self.to_process = [g.processed for g in groups]  # Dataset objects
            self.processed = None  # Dataset object, when completed
        else:
            self.n_dims = n_dims    # TODO: Could set to parse files to infer 
                                    # if not given.
            self.name = name if name else f'{str(n_dims)}d'
            # If no action is specified, will only be processed for CF 
            # compliance, but files won't be merged/split.
            self.action = action.lower() if action else None
            # if self.action == 'split':
            #     self.process = self.split_times
            # elif self.action == 'merge':
            #     self.process = self.merge_times
            # else:
            #     self.process = self.keep_time_series

            self.filepaths = filepaths if filepaths else []
            # Cannot set default valueas an empty list in arguments, or it 
            # assigns SAME empty list to every instance.

            self.processed = None  # Filepath(s) of output
            self.stem = stem_str(*self.filepaths) if self.filepaths else None

            # Although it would be best to verify that all member datasets have 
            # the same time variable, for now, just use the first one.
            try:
                self.time_var = time_variable or get_time_var(self.n_dims)
            except KeyError:
                if self.filepaths:
                    try:
                        self.time_var = get_time_var(
                            xr.open_dataset(
                            self.filepaths[0], decode_times=False
                            ))
                    except:
                        self.time_var = None    # Search in datasets as files 
                                                # are added.       
                else: 
                    self.time_var = None    # Search in datasets as files are 
                                            # added.
            # self.datasets is a generator, so any changes made to a member
            # object must be saved elsewhere, e.g. self.processed or NC file.
            # self.datasets = ds_feed(self.filepaths)
                # (xr.open_dataset(
                # path, decode_times=False
                # ) for path in self.filepaths)

    def datasets(self):
        return ds_feed(self.filepaths)
    
    def add(self, filepath: str = '') -> None:
        
        if not op.exists(filepath):
            raise OSError(f'DsGroup.add: Filepath {filepath} not found.')
        
        # Add filepath to collection
        self.filepaths.append(filepath)
        
        # Update name stem common to all files in group.
        self.stem = stem_str(self.stem, filepath)

        # print(f"Added {op.basename(filepath)} to {self.name}; stem: {self.stem}")

        # If don't have a time variable yet, attempt to find in new file.
        if not self.time_var:
            try:
                self.time_var = get_time_var(
                    xr.open_dataset(filepath), 
                    decode_times=False
                    )
            except KeyError:
                print(f'DsGroup.add: Time variable not found in {filepath}.')
    
    def merge_times(self) -> None:
        pass

    def merge_groups(self) -> None:
        pass


class TimeUnits(Units):
    '''
    Uses proleptic_gregorian as default calendar, as per ISO 8601:2004,
    even though CF standard calendar is a mixed Gregorian/Julian.
    '''
    def __init__(self, units:str=None, calendar:str='proleptic_gregorian', formatted=False, names=False, definition=False, _ut_unit=None):
        super().__init__(units, calendar, formatted, names, definition, _ut_unit)
        if not self.isvalid:
            raise ValueError('Units not valid according to cfunits module.')
    
    def time_unit(self) -> str:
        return self.formatted().split(' since ')[0]
    
    def since(self) -> str:
        return self.reftime.isoformat()
    
    def cf(self) -> str:
        return self.formatted()
    
    def base_units_match(self, other: Units) -> bool:
        return self.formatted().split(' since ')[0] == other.formatted().split(' since ')[0]
    
    def base_units_equivalent(self, other) -> bool:
        return Units(units=self.formatted().split(' since ')[0]).equivalent(
            Units(units=other.formatted().split(' since ')[0])
            )
    
    def ref_dates_match(self, other: Units) -> bool:
        return self.reftime == other.reftime
    
    def has_calendar(self) -> bool:
        try:
            return self.calendar is not None
        except AttributeError:
            return False

    def calendars_match(self, other: Units) -> bool:
        try:
            return self.calendar == other.calendar
        except AttributeError:
            return False
        except Exception as e:
            print('calendars_match function:', e)
            return False


# class DirectoryParser:
    # def __init__(self,
    #              directory: str = '',
    #              target: str = '') -> None:
    #     '''
    #     directory:          path of directory containing NC files to be made CF 
    #                         compliant.
    #     target (optional):  path to which processed files should be written.
    #     '''
        # if directory and not op.exists(directory):
        #     raise OSError(f'Directory {directory} not found.')
        # self.directory = directory or os.getcwd()
        # self.target_dir = target or op.join(op.dirname(
        #       self.directory), f'{op.basename(self.directory)}+processed'
        #     )
        # # If target directory doesn't exist, create it.
        # if not op.exists(self.target_dir):
        #     print(f'{self.target_dir} does not yet exist: making...')
        #     try:
        #         os.makedirs(self.target_dir)
        #         print(f'{self.target_dir} now exists? {op.exists(self.target_dir)}')
        #     except Exception as e:
        #         print(e)
        #         raise e
        # self.input_files = []
        # # Set up dict to categorize MONC outputs by number of spatial dimensions
        # self.by_dim = {n_dim: DsGroup(
        #     name=group, n_dims=n_dim, action=DIM_ACTIONS[group]
        # ) for n_dim, group in DIM_GROUPS.items()}
        # # ID time variable of each group from VOCAB.
        
        # # Parse directory & categorise NC files by type and dimensions
        # self.sort_nc()

        # # Attempt to derive time unit info.
        # # Assume time units will be common to all.
        # self.time_units = self.time_units_from_input() if self.input_files else None

    # @performance_time
    # def sort_nc(self) -> None:
    #     '''
    #     Sorts into MONC output files and other NC files, the latter listed as
    #     possible input files.

    #     It also categorizes the MONC outputs according to their number of 
    #     dimensions.

    #     It would be best if all files of a given run were in one directory, 
    #     rather than split between 0-2d and 3d.
    #     '''
    #     print('Categorising fields by dimension')
    #     # For each NC file in source directory (recursive parsing or not?):
    #     for filepath in iglob(f'{op.join(self.directory, "*.nc")}', 
    #                       root_dir=os.pathsep, 
    #                       recursive=False):
    #         with xr.open_dataset(filepath, 
    #                              decode_times=False) as ds:  # concat_characters=False
    #             # Using concat_characters=False preserves the string dimension 
    #             # of the options_database variable. However, it makes dealing
    #             # with the binary -> string conversion painful and messy.

    #             # Categorise MONC output / other
    #             if is_monc(ds):
    #                 # Find number of (non-options-database) dimensions.
    #                 n_dims = get_n_dims(ds) - 1  # Subtract 1 for time.
    #                 # Add filepath to relevant dimension group.
    #                 self.by_dim[n_dims].add(filepath=filepath)
    #                 # print(f"{op.basename(filepath)}: {n_dims} spatial dimensions.")
    #             else:
    #                 # Categorise as potential input file
    #                 self.input_files.append(filepath)
    #                 # print(f"{op.basename}: possible input file.")
    
    

    # @performance_time
    # def process_by_dim(self, dim: int = None) -> str:

    #     processes = []

    #     to_process = {dim: self.by_dim[dim]} if (
    #         dim and dim in self.by_dim
    #         ) else self.by_dim
        
    #     for dim, group in to_process.items():

    #         print(f"Processing {dim}:{group.name}.")
            
    #         # If group is to be merged:
    #         if group.action == 'merge':
    #             processing = list(group.datasets())  # Open all datasets in group.
    #             for ds in processing:
    #                 # Assign any required attributes to variables
    #                 '''
    #                 Assume these are only correct for last time-point in series.
    #                 '''

    #                 # Assign any 'one per group' global attributes to temporary 
    #                 # variables, if larger in current file than current values.

    #             # Merge all datasets in group along time variable, & assign 
    #             # resulting dataset to new group attribute. Use chunking if 
    #             # large.

    #             # Close individual datasets, keeping only resultant.

    #             # Apply 'once per group' attributes back to resulting dataset.

    #             # Derive title/filename from group's stem & dimension
                
    #             # Call CF compliance function on dataset. 
                
    #             # Set filepath and save as self.processed

    #             # Set encoding

    #             # Save dataset to NetCDF

    #             # Close dataset

    #             pass

    #         # If not merging:
    #         else:
    #             # Derive base title/filename from group's stem & dimension
    #             title = group.stem + group.name if group.name not in group.stem else group.stem
                
    #             writers = []
                
    #             if group.action == 'split':
    #                 # For each filepath in group, process in a separate Process
    #                 processes += [Process(
    #                     target=process_large, 
    #                     kwargs={'filepath':path, 
    #                             'dim':dim,
    #                             'title':title, 
    #                             'time_var':group.time_var, 
    #                             'target_dir':self.target_dir}
    #                     ) for path in group.filepaths[1:]]
    #                 [p.start() for p in processes[-2:]]  # TODO: this needs to know how many time points are in the file.

    #                 # Continue on local process
    #                 process_large(
    #                     filepath=group.filepaths[0],
    #                     dim=dim,
    #                     title=title,
    #                     time_var=group.time_var,
    #                     target_dir=self.target_dir
    #                 )
    #             else:
    #                 # For each filepath in group:
    #                 for i, path in enumerate(group.filepaths[:4]):
    #                     # Open as dataset
    #                     with xr.open_dataset(path, 
    #                                         decode_times=False, 
    #                                         chunks={'z': 11, 'zn':11}) as ds:

    #                         # Update dataset's title
    #                         ds.attrs['title'] = title
                        
    #                         # apply chunking if large (e.g. ds.nbytes >= 1e9).

    #                         #Call CF compliance function
                            
    #                         if group.action == 'split':
    #                             # Split dataset by time-point, yielding multiple 
    #                             # new datasets.
    #                             group.processed = split_ds(dataset=ds,
    #                                                     var=group.time_var)
    #                             # group.processed only ever needs to hold latest
    #                             # collection of datasets.

    #                             ds.close()
    #                         else:
    #                             # Copy dataset as-is into self.processed
    #                             group.processed = [ds]  # This should be wrapped in a DsGroup function

    #                         # For each new dataset:
    #                         for new_ds in group.processed:
    #                             # Derive new title

    #                             # Update required global attributes (MONC time, title, 
    #                             # MONC timestep, (previous diagnostic write at?)).

    #                             # set filepath
    #                             filepath = op.join(
    #                                 self.target_dir, f"{new_ds.attrs['title']}.nc"
    #                                 )
    #                             # set encodings
    #                             encodings = {
    #                                 k:{
    #                                     'dtype': v.dtype,
    #                                     '_FillValue': None
    #                                 } for k, v in new_ds.variables.items()
    #                             }

    #                             # Export to NetCDF (set as single-command function, 
    #                             # so can wrap in performance_time).
    #                             # TODO: Use compute=False option, then call 
    #                             # compute() on resulting dask.delayed.Delayed 
    #                             # object?
    #                             if i<2:
    #                                 # Test regular save function
    #                                 print(f"Saving new dataset as {filepath}.")
    #                                 ds_to_nc(ds=new_ds,
    #                                         filepath=filepath,
    #                                         encodings=encodings)
    #                             else:
    #                                 # Test dask version
    #                                 print(f"Preparing delayed writers for {filepath}.")
    #                                 writers.append(ds_to_nc_dask(ds=new_ds,
    #                                         filepath=filepath,
    #                                         encodings=encodings))
    #                             # new_ds.to_netcdf(
    #                             #     path=filepath, 
    #                             #     encoding=encodings)

    #                             # Close new dataset: NEED TO CHECK 
    #                             # dask.delayed.Delayed object doesn't need dataset 
    #                             # to remain open until it completes.
    #                             new_ds.close()

    #                             # Run CF checker on output file

    #                         if i > 1:
    #                             print(f"Computing last {len(group.processed)} writers.")
    #                             [perform_write(writer) 
    #                             for writer in writers[-len(group.processed):]
    #                             ]
                                
    #                         # If worked on single dataset, close it.
    #                         if len(group.processed) == 1:
    #                             group.processed[0].close()

    #     [p.join() for p in processes]
    #     return # Success/fail message with resulting filepath.


def globals_to_vars(ds: xr.Dataset, 
                    time_var: str, 
                    last_time_point: int|float) -> xr.Dataset:
    
    vars = {}
    for global_attr, var_attrs in CONFIG['global_to_variable'].items():

        # Any advantage of name = global_attr.replace(' ', '_')?

        if global_attr not in ds.attrs:
            raise AttributeError
        
        if global_attr in do_not_propagate:
            # assign `np.nan` to each data point except last, and then 
            # set coords as per other case below.
            data = type_from_str(ds.attrs[global_attr])
            coords = [last_time_point]

        else:
            data = [type_from_str(ds.attrs[global_attr])]*len(ds[group.time_var].data)
            coords = ds[group.time_var]

        try:
            vars[global_attr] = xr.DataArray(
                name=global_attr,
                data=data,
                coords={time_var: coords},
                dims=(time_var,),
                attrs=var_attrs
            )
        except Exception as e:
            raise e
        
    return vars


def cf_merge(
        group: DsGroup, 
        time_units: TimeUnits, 
        target_dir: str
    ) -> str:
    '''
    Should print & return name of saved file(s)
    '''

    from cfize_ds import cfize_dataset  # To avoid circular imports
    
    print(f"Processing {group.n_dims}:{group.name}.")

    processing = list(group.datasets())  # Open all datasets in group.
    for i, ds in enumerate(processing):
        # Find dataset's time coordinate variable
        time_var = variable_glob(ds=ds, var_glob=group.time_var)
        if len(time_var) != 1:
            # TODO: disambiguation
            raise AttributeError(
                'Multiple variables found in dataset that match time '
                f'variable pattern, {group.time_var}.')
        else:
            time_var = time_var[0]

        # Assign any required global attributes to variables, 
        # with associated attributes. Remove global attribute.
        '''
        Assume these are only correct for last time-point in series.
        '''
        
        last_time_point = ds[time_var].data[-1]
        try:
            new_vars = globals_to_vars(
                        ds=ds, 
                        time_var=time_var, 
                        last_time_point=last_time_point
            )
        except Exception as e:
            raise e
        for name, array in new_vars.items():
            ds.attrs.pop(name)
            processing[i] = ds.assign({name: array})  # Need to assign to original dataset in list rather than placeholder ds.
            
        # Assign any 'one per group' global attributes to temporary 
        # variables, if larger in current file than current values.
        

    # Merge all datasets in group along time variable, & assign 
    # resulting dataset to new group attribute. Use chunking if 
    # large.

    # Close individual datasets, keeping only resultant.

    # Apply 'once per group' attributes back to resulting dataset.

    # Derive title/filename from group's stem & dimension
    
    # Call CF compliance function on dataset. 
    
    # Set filepath and save as processed
    processed = f"{target_dir}/filename.nc"

    # Set encoding

    # Save dataset to NetCDF

    # Close dataset

    return processed


def process_large(
        filepath: str, 
        dim: int, 
        title: str, 
        time_var: str, 
        time_units: TimeUnits, 
        target_dir: str, 
        split: bool
    ) -> str:
    '''
    Should return name of saved file(s)
    '''

    from cfize_ds import cfize_dataset
    
    print(op.basename(filepath), "- process_large - process id:", os.getpid())
    
    with xr.open_dataset(filepath, 
                         decode_times=False
                         ) as ds:
        
        # apply chunking if large (e.g. ds.nbytes >= 1e9).
        # chunks=CHUNKING_DIMS
        
        # Update dataset's title
        ds.attrs['title'] = title
    
        

        #Call CF compliance function
        ds = cfize_dataset(
            dataset=ds,
            dimension=dim,
            title=title,
            time_units=time_units
        )
            
        
        if split:
            # Split dataset by time-point, yielding multiple 
            # new datasets.
            processed = split_ds(dataset=ds, 
                                var=time_var)
            # processed only ever needs to hold latest collection of datasets.
        else:
            processed = ds  # Will this persist after context of ds ends?

    writers = []
    # For each new dataset:
    for ds in processed:
        # Derive new title

        # Update required global attributes (MONC time, title, 
        # MONC timestep, (previous diagnostic write at?)).

        # set filepath
        filepath = op.join(
            target_dir, f"{ds.attrs['title']}.nc"
            )
        # set encodings
        encodings = {
            k:{
                'dtype': v.dtype,
                '_FillValue': None
            } for k, v in ds.variables.items()
        }

        # Export to NetCDF (set as single-command function, 
        # so can wrap in performance_time).
        # TODO: Use compute=False option, then call 
        # compute() on resulting dask.delayed.Delayed 
        # object?
        # if i == 0:
        #     # Test regular save function
        #     print(f"Saving new dataset as {filepath}.")
        #     ds_to_nc(ds=ds,
        #             filepath=filepath,
        #             encodings=encodings)
        # else:
        # Test dask version
        print(f"Preparing delayed writer for {filepath}.")
        writers.append(ds_to_nc_dask(ds=ds,
                                filepath=filepath,
                                encodings=encodings
                                ))
        ds.close()
    print(f"Computing writes.")
    [perform_write(writer) for writer in writers]
    
    # Run CF checker on output files
    pass


def process_parallel(
        groups: dict, 
        target_dir: str, 
        n_proc: int, 
        time_units: TimeUnits = None
):
    '''
    groups:     dimension-specific groups to be processed.
    target_dir: directory in which to write processed NC files.
    n_proc:     number of processes in process pool (doesn't include controller)
    time_units: CF-compliant unit to use for time coordinate.
    '''

    # Verify n_proc + controller doesn't exceed available cores
    # TODO: confirm this works on Jasmin!
    if n_proc + 1 > os.cpu_count(): n_proc = os.cpu_count() - 1

    results = {dim: [] for dim in groups.keys()}  # set up empty dict for results from process pool

    # Set up process pool
    with Pool(processes=n_proc) as pool:
        
        # Allocate datasets to be processed in increasing order of dimensions
        # Work with 0-2d files first, because 3d processing depends on 0d/1d.
        '''
        ENHANCEMENT: this processing in order of increasing dimension is 
        critical for converting perturbations to absolutes. For non-MONC data, 
        it may not be safe to assume that reference variables will be in
        lower-dimension datasets, so waiting for their completion would hang 
        execution or fail. A safer alternative would be to export any reference 
        variables to a standalone NetCDF file, which can then be accessed when 
        needed. In this case, there should be a flag to indicate when the
        required variable is available.
        '''
        for dim in sorted(list(groups.keys())):
            group = groups[dim]
            print(f"Processing {dim}:{group.name}.")
            
            if group.action == 'split' or dim == 3:
                # The 2nd condition ensures that large files are allocated to
                # the process pool, whether they are to be split or not.
                
                # Derive base title/filename from group's stem & dimension
                title = group.stem + group.name if group.name not in group.stem else group.stem
                
                for_controller = []

                # Wait for 0+1d processing to finish, so reference variables
                # are available for perturbations.
                for d in set(reference_vars.values()):
                    groups[d].processed = [
                        result.get() for result in results[d]
                    ]
                
                # Round-robin allocation of file processing to controller and 
                # workers:
                for i, filepath in enumerate(group.filepaths):
                    # Allocate one file per round to controller
                    if i%(n_proc+1) == 0:
                        for_controller.append(filepath)
                    else:
                        results[dim].append(
                            pool.apply_async(
                                func=process_large, 
                                kwds={
                                    'filepath':filepath, 
                                    'dim':dim,
                                    'title':title, 
                                    'time_var':group.time_var, 
                                    'time_units':time_units,
                                    'target_dir':target_dir,
                                    'split':group.action == 'split'
                                    }
                                )
                            )
                
                # Now work on remaining files on controller process
                group.processed = [process_large(
                        filepath=filepath,
                        dim=dim,
                        title=title,
                        time_var=group.time_var,
                        time_units=time_units,
                        target_dir=target_dir,
                        split=group.action == 'split'
                    ) for filepath in for_controller]

                # Gather completed jobs
                [group.processed.append(result.get()) 
                    for result in results[dim]]

            # If group is to be merged:
            elif group.action == 'merge':
                # results[dim] = pool.apply_async(
                #     func=cf_merge,
                #     kwds={
                #         'group': group,
                #         'time_units': time_units,
                #         'target_dir': target_dir
                #     }
                # )
                results[dim].append(
                    pool.apply_async(
                    func=cf_merge,
                    kwds={
                        'group': group,
                        'time_units': time_units,
                        'target_dir': target_dir
                    })
                )
        
            # Otherwise, process for CF compliance, but leave dataset as a 
            # single, standalone file.
            else:
                # ENHANCEMENT: Could still evaluate by size whether it's 
                # processed locally or passed to the Pool.

                # Round-robin allocation of file processing to controller and 
                # workers:
                
                pass

        # Check all other processes are complete
        # for result in results:
        #     print(result.get())  # assuming here that and DirectoryParser.process_by_dim and process_large each return values.
        for dim in list(groups.keys()).sort:
            # Sort because expect smaller to finish first
            if not group[dim].processed:
                group.processed = [result.get() for result in results[dim]]


def process_serial(
        groups: dict, 
        target_dir: str, 
        time_units: TimeUnits = None
):
    '''
    groups:     dimension-specific groups to be processed.
    target_dir: directory in which to write processed NC files.
    time_units: CF-compliant unit to use for time coordinate.
    '''

    # Processe datasets in increasing order of dimensions.
    # Work with 0-2d files first, because 3d processing depends on 0d/1d.
    '''
    ENHANCEMENT: this processing in order of increasing dimension is 
    critical for converting perturbations to absolutes. For non-MONC data, 
    it may not be safe to assume that reference variables will be in
    lower-dimension datasets, so waiting for their completion would hang 
    execution or fail. A safer alternative would be to export any reference 
    variables to a standalone NetCDF file, which can then be accessed when 
    needed. In this case, there should be a flag to indicate when the
    required variable is available.
    '''
    for dim in sorted(list(groups.keys())):
        group = groups[dim]
        print(f"Processing {dim}:{group.name}.")
        
        if group.action == 'split':
            
            # Derive base title/filename from group's stem & dimension
            title = group.stem + group.name if group.name not in group.stem else group.stem
            
            group.processed = [process_large(
                    filepath=filepath,
                    dim=dim,
                    title=title,
                    time_var=group.time_var,
                    time_units=time_units,
                    target_dir=target_dir,
                    split=group.action == 'split'
                ) for filepath in group.filepaths]

        # If group is to be merged:
        elif group.action == 'merge':
            group.processed = cf_merge(
                    group=group,
                    time_units=time_units,
                    target_dir=target_dir
            )
    
        # Otherwise, process for CF compliance, but leave dataset as a 
        # single, standalone file.
        else:
            pass


@performance_time
def sort_nc(directory) -> [dict, list]:
    '''
    Sorts into MONC output files and other NC files, the latter listed as
    possible input files.

    It also categorizes the MONC outputs according to their number of 
    dimensions.

    It would be best if all files of a given run were in one directory, 
    rather than split between 0-2d and 3d.
    '''
    print('Categorising files by dimension')

    by_dim = {n_dim: DsGroup(
        name=group, n_dims=n_dim, action=DIM_ACTIONS[group]
    ) for n_dim, group in DIM_GROUPS.items()}

    input_files = []

    # For each NC file in source directory (recursive parsing or not?):
    for filepath in iglob(f'{op.join(directory, "*.nc")}', 
                          root_dir=os.pathsep, 
                          recursive=False):
        with xr.open_dataset(filepath, 
                             decode_times=False) as ds:  
                            # concat_characters=False
            # Using concat_characters=False preserves the string dimension 
            # of the options_database variable. However, it makes dealing
            # with the binary -> string conversion painful and messy.

            # Categorise MONC output / other
            if is_monc(ds):
                # Find number of (non-options-database) dimensions.
                n_dims = get_n_dims(ds) - 1  # Subtract 1 for time.
                # Add filepath to relevant dimension group.
                by_dim[n_dims].add(filepath=filepath)
                # print(f"{op.basename(filepath)}: {n_dims} spatial dimensions.")
            else:
                # Categorise as potential input file
                input_files.append(filepath)
                # print(f"{op.basename}: possible input file.")
    
    return [by_dim, input_files]


@performance_time
def time_units_from_input(filepaths: list) -> TimeUnits:
    '''
    Attempt to find time unit data in possible input file(s).
    Refers to FROM_INPUT_FILE['reftime'] from setup.py, for variable name
    under which time units are expected.
    '''

    input_data = {}
    for f in filepaths:
        with xr.open_dataset(f, decode_times=False) as ds:  
            # , concat_characters=False
            try:
                input_data[f] = ds[FROM_INPUT_FILE['reftime']].attrs
            except KeyError:
                # Variable not found in prospective input file; move on to
                # next input file, if any.
                continue
    if len(input_data) != 1:
        # Request user input for tie-breaking.
        if len(input_data) > 1:
            print('Multiple time units found:')
            for i, u in enumerate(input_data.values()):
                print(f'[{i}] {u}')
            while True:
                user_input = input('Please select number or enter units for time, including reference date in ISO format.')
                if user_input.isnumeric():
                    try:
                        input_data = list(input_data.values())[int(user_input)]
                    except:
                        print('Invalid number; try again.')
                        continue
                    else:
                        break
                else:
                    input_data = {'units': user_input}
        else:
            input_data = {
                'units': 
                input('Enter units for time, including reference date in ISO format:')
            }
    else:
        input_data = list(input_data.values())[0]
    while True:
        try:
            time_units = TimeUnits(
                units=input_data['units'], 
                calendar=input_data['calendar'] if 'calendar' in input_data else DEFAULT_CALENDAR
            )
        except ValueError or KeyError as e:
            input_data['units'] = input(e + 'Enter units for time, including reference date in ISO format:')
            continue
        else:
            return time_units
            
            
def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'source_dir',
        help='Directory containing MONC output NetCDF files to process.'
    )
    parser.add_argument(
        '--target_dir', '-t', 
        dest='target_dir', 
        help='Directory for processed NetCDF files.',
        required=False
    )
    parser.add_argument(
        '--reference_time', '-r',
        dest='ref_time',
        help='Date or date-time of origin for time units, in ISO format (yyyy-mm-dd[Thh:mm:ss][{+/-UTC offset, hh:mm}]).',
        required=False
    )
    parser.add_argument(
        '--calendar', '-c',
        dest='calendar',
        choices={'gregorian', 'standard', 'none', 'proleptic_gregorian',
                 '360_day', 'noleap', '365_day', 'all_leap', '366_day',
                 'julian'},
        help='Calendar to use for time units.',
        required=False
    )
    parser.add_argument(
        '--cpus', '-p',
        dest='n_proc',
        type=int,
        choices=range(1, 17),
        help='Total number of cores to use for parallel processing (including controller process).',
        required=False
    )
    return parser.parse_args()


def main():
    
    print("Main app process id:", os.getpid())
    
    # get command line arguments
    args = parse_arguments()
    
    # Validate supplied directory to parse.
    if not op.exists(args.source_dir):
        raise OSError(f'Source directory {op.abspath(source_dir)} not found.')
    source_dir = op.abspath(args.source_dir)
    
    time_units = None
    # Validate & set reference time if supplied.
    if args.ref_time:
        try:
            ref_time = datetime(args.ref_time)
        except Exception as e:
            raise TypeError(f"{args.ref_time} is an invalid reference time. {e}")
    else:
        ref_time = None

    # Validate calendar and if possible set time units.
    calendar = DEFAULT_CALENDAR
    if args.calendar:
        try:
            calendar = Units(calendar=args.calendar).calendar
        except Exception as e:
            print(
                f"Calendar {args.calendar} not valid."
            )  # This should not be possible
        if ref_time:
            try:
                time_units = TimeUnits(
                    units=default_time_unit + ' since ' + ref_time.isoformat(),
                    calendar=calendar)
            except ValueError as e:
                print(
                    f"Could not create valid time units from reference date-time {ref_time} and calendar {calendar}. {e}"
                )
    
    # Get/set directory for output files.
    target_dir = op.abspath(args.target_dir) if args.target_dir else op.abspath(
        op.join(op.dirname(source_dir), 
                f'{op.basename(source_dir)}+processed')
    )
    # If target directory doesn't exist, create it.
    if not op.exists(target_dir):
        try:
            os.makedirs(target_dir)
        except OSError as e:
            raise OSError(
                f"Unable to create target directory, {target_dir}. {e}"
            )
    if not os.access(target_dir, os.W_OK):
        raise OSError(
            f"Write permission denied for target directory, {target_dir}."
        )
    if not args.target_dir:
        print(f"Processed files will be saved to: {target_dir}")

    # Set number of processes to place in process pool.
    # Subtract 1 from number of processes specified, to use one as controller.
    n_proc = args.n_proc - 1 if args.n_proc else 0
    if n_proc > os.cpu_count() - 1:
        raise OSError(
            f"Not enough cores available for {args.n_proc} processes. Maximum: {os.cpu_count()}."
        )
        # TODO: check os.cpu_count() works correctly on JASMIN, when mutliple
        # cores are allocated.
    
    # Parse directory & categorise NC files by type and dimensions
    [group_by_dim, input_files] = sort_nc(source_dir)
    
    # ID time variable of each group from VOCAB.
    
    # Attempt to derive time unit info, if not already supplied at command line.
    # Assume time units will be common to all datasets in directory.
    if not time_units:
        time_units = time_units_from_input(input_files) if input_files else None

    # NOTE: by this stage, time_units should contain a valid reference date /
    # datetime and calendar. The base units will be set to the default 
    # 'seconds', but this will be overridden by any units found in time 
    # coordinate variable(s).
            
    # For each dimension group:
    # If n_proc>0, use process pool; otherwise process sequentially.
    if n_proc > 0:
        process_parallel(
            groups=group_by_dim, 
            time_units=time_units, 
            n_proc=n_proc, 
            target_dir=target_dir
        )
    else:
        process_serial(
            groups=group_by_dim, 
            time_units=time_units, 
            target_dir=target_dir
        )

    # parser.process_by_dim()
    # parser.process_by_dim(dim=3)  # Test split+save components alone.

    # Find dimension groups to be merged

    # For each set of dimension groups to be merged:

        # Create new group, comprising the groups to be merged.

        # Check all datasets to be merged have matching time coordinates.
        # If not, give error message and keep datasets separate.

        # Merge each group's resultant single dataset:

            # Open each group.processed as xr.Dataset

            # Derive combined name from stem of each group's filename stem and 
            # dimension.

            # Check & identify common time variable.

            # Get dataset from each group's processed attribute OR load from
            # saved interim file.

            # Merge dataset of each group:

                # Take maximum value of an 'one per group' attributes.

                # xarray.Dataset.merge -> new dataset.

                # Close source datasets

                # Add 'one per group' attributes to new dataset.

            # Assign new title to new dataset

        # Set filepath

        # Set encoding

        # Save dataset to NetCDF

        # Close dataset

        # Delete interim NC files, unless flagged to do otherwise.

    # Garbage collection if required.

    exit(0)


if __name__ == '__main__': main()