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

            self.processed = []
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
            self.datasets = (xr.open_dataset(
                path, decode_times=False
                ) for path in self.filepaths)

    def add(self, filepath: str = '') -> None:
        
        if not op.exists(filepath):
            raise OSError(f'DsGroup.add: Filepath {filepath} not found.')
        
        # Add filepath to collection
        self.filepaths.append(filepath)
        
        # Update name stem common to all files in group.
        self.stem = stem_str(self.stem, filepath)

        print(f"Added {op.basename(filepath)} to {self.name}; stem: {self.stem}")

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


class DirectoryParser:
    def __init__(self,
                 directory: str = '',
                 target: str = '') -> None:
        '''
        directory:          path of directory containing NC files to be made CF 
                            compliant.
        target (optional):  path to which processed files should be written.
        '''
        if directory and not op.exists(directory):
            raise OSError(f'Directory {directory} not found.')
        self.directory = directory or os.getcwd()
        self.target_dir = target or op.join(op.dirname(
              self.directory), f'{op.basename(self.directory)}+processed'
            )
        # If target directory doesn't exist, create it.
        if not op.exists(self.target_dir):
            print(f'{self.target_dir} does not yet exist: making...')
            try:
                os.makedirs(self.target_dir)
                print(f'{self.target_dir} now exists? {op.exists(self.target_dir)}')
            except Exception as e:
                print(e)
                raise e
        self.input_files = []
        # Set up dict to categorize MONC outputs by number of spatial dimensions
        self.by_dim = {n_dim: DsGroup(
            name=group, n_dims=n_dim, action=DIM_ACTIONS[group]
        ) for n_dim, group in DIM_GROUPS.items()}
        # ID time variable of each group from VOCAB.
        
        # Parse directory & categorise NC files by type and dimensions
        self.sort_nc()

        # Attempt to derive time unit info.
        # Assume time units will be common to all.
        self.time_units = self.time_units_from_input() if self.input_files else None

    def sort_nc(self) -> None:
        '''
        Sorts into MONC output files and other NC files, the latter listed as
        possible input files.

        It also categorizes the MONC outputs according to their number of 
        dimensions.

        It would be best if all files of a given run were in one directory, 
        rather than split between 0-2d and 3d.
        '''
        print('Categorising fields by dimension')
        # For each NC file in source directory (recursive parsing or not?):
        for filepath in iglob(f'{op.join(self.directory, "*.nc")}', 
                          root_dir=os.pathsep, 
                          recursive=False):
            with xr.open_dataset(filepath, 
                                 decode_times=False) as ds:  # concat_characters=False
                # Using concat_characters=False preserves the string dimension 
                # of the options_database variable. However, it makes dealing
                # with the binary -> string conversion painful and messy.

                # Categorise MONC output / other
                if is_monc(ds):
                    # Find number of (non-options-database) dimensions.
                    n_dims = get_n_dims(ds) - 1  # Subtract 1 for time.
                    # Add filepath to relevant dimension group.
                    self.by_dim[n_dims].add(filepath=filepath)
                    # print(f"{op.basename(filepath)}: {n_dims} spatial dimensions.")
                else:
                    # Categorise as potential input file
                    self.input_files.append(filepath)
                    # print(f"{op.basename}: possible input file.")
        
    def time_units_from_input(self) -> TimeUnits:
        '''
        Attempt to find time unit data in possible input file(s).
        Refers to FROM_INPUT_FILE['reftime'] from setup.py, for variable name
        under which time units are expected.
        '''

        input_data = {}
        for f in self.input_files:
            with xr.open_dataset(op.join(self.directory, f), decode_times=False) as ds:  # , concat_characters=False
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

    def process_by_dim(self, dim: int = None) -> None:

        to_process = {dim: self.by_dim[dim]} if (
            dim and dim in self.by_dim
            ) else self.by_dim
        
        for dim, group in to_process.items():

            print(f"Processing {dim}:{group.name}.")
            
            # If group is to be merged:
            if group.action == 'merge':
                processing = list(group.datasets)  # Open all datasets in group.
                for ds in processing:
                    # Assign any required attributes to variables
                    '''
                    Assume these are only correct for last time-point in series.
                    '''

                    # Assign any 'one per group' global attributes to temporary 
                    # variables, if larger in current file than current values.

                # Merge all datasets in group along time variable, & assign 
                # resulting dataset to new group attribute. Use chunking if 
                # large.

                # Close individual datasets, keeping only resultant.

                # Apply 'once per group' attributes back to resulting dataset.

                # Derive title/filename from group's stem & dimension
                
                # Call CF compliance function on dataset. 
                
                # Set filepath and save as self.processed

                # Set encoding

                # Save dataset to NetCDF

                # Close dataset

                pass

            # If not merging:
            else:
                # Derive base title/filename from group's stem & dimension
                title = group.stem + group.name if group.name not in group.stem else group.stem
                
                # For each filepath in group:
                for path in group.filepaths:

                    # Open as dataset
                    with xr.open_dataset(path, decode_times=False) as ds:

                        # Update dataset's title
                        ds.attrs['title'] = title
                    
                        # apply chunking if large (e.g. ds.nbytes >= 1e9).

                        #Call CF compliance function
                        
                        if group.action == 'split':
                            # Split dataset by time-point, yielding multiple 
                            # new datasets.
                            group.processed = split_ds(dataset=ds,
                                                       var=group.time_var)
                            # group.processed only ever needs to hold latest
                            # collection of datasets.

                            ds.close()
                        else:
                            # Copy dataset as-is into self.processed
                            group.processed = [ds]  # This should be wrapped in a DsGroup function

                        # For each new dataset:
                        for new_ds in group.processed:
                            # Derive new title

                            # Update required global attributes (MONC time, title, 
                            # MONC timestep, (previous diagnostic write at?)).

                            # set filepath
                            filepath = op.join(
                                self.target_dir, f"{new_ds.attrs['title']}.nc"
                                )
                            print(f"Saving new dataset as {filepath}.")

                            # set encodings
                            encodings = {
                                k:{
                                    'dtype': v.dtype,
                                    '_FillValue': None
                                } for k, v in new_ds.variables.items()
                            }

                            # Export to NetCDF (set as single-command function, so 
                            # can wrap in performance_time). Use compute=False 
                            # option, then call compute() on resulting 
                            # dask.delayed.Delayed object.
                            new_ds.to_netcdf(
                                path=filepath, 
                                encoding=encodings)

                            # Close new dataset: NEED TO CHECK 
                            # dask.delayed.Delayed object doesn't need dataset 
                            # to remain open until it completes.
                            new_ds.close()

                            # Run CF checker on output file

                        # If worked on single dataset, close it.
                        if len(group.processed) == 1:
                            group.processed[0].close()


def main():
    # TODO: Validate supplied arguments/directory

    # TODO: Set directory to parse & target directory.
    source = '/gws/nopw/j04/eurec4auk/monc_prelim_output/jan_28_3d'
    target = '~/cfizer/testing'
    # source = os.path.join(os.path.dirname(app_dir), 'test_data')
    # target = None
    try:
        print("Creating directory parser, to process", source)
        parser = DirectoryParser(directory=source,
                                 target=target)
    except OSError as e:
        print(e)
        exit(1)
    except Exception as e:  # Catch anything not anticipated
        print(e)
        exit(1)
            
    # For each dimension group:
    # parser.process_by_dim()
    parser.process_by_dim(dim=3)  # Test split+save components alone.

    # For each set of dimension groups to be merged:

        # Create new group, comprising the groups to be merged.

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