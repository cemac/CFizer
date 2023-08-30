from difflib import SequenceMatcher
import os.path as op
import os
import xarray as xr
from glob import iglob
from setup import *
import re
from units import TimeUnits
from cfize_ds import MoncDs
from utils import type_from_str


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
    
    '''
    # If argument is an integer, assume it's the number of spatial dimensions
    if isinstance(args[0], int):
        if args[0] not in vocabulary:
            raise KeyError(
                'get_time_var: Vocabulary file has no specifications for '
                f'{args[0]} spatial dimensions.')
        time_vars = [v for v in vocabulary[args[0]].keys() 
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


def merge_dimensions(*args) -> xr.Dataset:
    '''
    This is designed to merge datasets whose time coordinates match, but contain
    mutually exclusive sets of data variables.

    Each argument should be an xarray.Dataset object.
    '''
    
    if not all(isinstance(arg, xr.Dataset) for arg in args):
        raise TypeError('All parameters must be of type xarray.Dataset')
    
    if len(args) == 1:  # Should never happen
        return args[0].copy(deep=True)
    
    if len(args) == 2:
        try:
            new_ds = args[0].merge(
                other=args[1], 
                join='exact', 
                combine_attrs='drop_conflicts'
            )
        except xr.MergeError as e:
            raise e
        return new_ds
    
    # call recursively until only 2 args to process
    return merge_dimensions(args[0], merge_dimensions(*args[1:]))


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


class DsGroup:
    '''
    Each file in a DsGroup has the same variables, coordinates & dimensions.
    '''
    def __init__(self, 
                 shared: dict,
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

        # global vocabulary, reference_vars

        # If groups parameter contains existing DsGroup objects:
        if groups and all([isinstance(g, DsGroup) for g in groups]):
            # check each DsGroup has only one dataset in its processed 
            # attribute.
            if not all([
                (isinstance(g.processed, list) and 
                len(g.processed) == 1 and 
                g.processed[0][-3:] == '.nc')
                for g in groups
            ]):
                raise TypeError('DsGroup: To merge DsGroup objects, the '
                                'processed attribute of each must contain a '
                                'list whose single element is a string '
                                'containing a NetCDF file path.')
            
            self.filepaths = [g.processed[0] for g in groups]
            
            # Open each group.processed as xr.Dataset
            try:
                if shared['verbose']: print(
                    f"Process {os.getpid()}: DsGroup: opening datasets {[op.basename(path) for path in self.filepaths]}."
                )
                self.to_process = [
                    xr.open_dataset(path, decode_times=False) 
                    for path in self.filepaths
                ]  # Dataset objects
            except OSError:
                raise OSError(
                    f"DsGroup: Invalid filepath in processed attribute(s) of "
                    f"group(s) passed as argument(s): {self.filepaths}."
                )
            except Exception as e:
                raise e
            
            self.name = stem_str(*[group.name for group in groups])
            self.stem = stem_str(*[g.stem for g in groups])
            # Derive combined name from stem of each group's filename stem and 
            # name.
            self.name = self.stem.strip(' _') + '_' + self.name

            
            # Check & identify common time variable.
            self.n_dims = sum([g.n_dims for g in groups])  # TODO: This is a little dubious, given there's no guarantee dimension groups don't overlap, but works for MONC outputs where 0d & 1d are grouped.
            
            self.time_var = time_variable or groups[0].time_var
            # Check all datasets to be merged have matching time coordinates.
            if not all([g.time_var == self.time_var for g in groups]):
                # If not, give error message and keep datasets separate.
                print(
                    f'DsGroup: When attempting to combine groups, '
                    f'found mismatch in time variables: '
                    f'{[g.time_var for g in groups]}.'
                )
                self.action = None
                # raise ValueError('DsGroup: When attempting to combine groups, '
                #                  'found mismatch in time variables.')
            else:
                self.action = action.lower() if action else 'merge_groups'
            self.processed = None  # Filepath(s) of output
        else:
            self.n_dims = n_dims    # TODO: Could set to parse files to infer 
                                    # if not given.
            self.name = name if name else f'{str(n_dims)}d'
            # If no action is specified, will only be processed for CF 
            # compliance, but files won't be merged/split.
            self.action = action.lower() if action else None
            self.filepaths = filepaths if filepaths else []
            # Cannot set default valueas an empty list in arguments, or it 
            # assigns SAME empty list to every instance.

            self.to_process = None
            self.processed = None  # Filepath(s) of output
            self.stem = stem_str(*self.filepaths) if self.filepaths else None

            # Although it would be best to verify that all member datasets have 
            # the same time variable, for now, just use the first one.
            try:
                self.time_var = time_variable or get_time_var(self.n_dims)
            except KeyError:
                if self.filepaths:
                    # Although it would be best to verify that all member 
                    # datasets have the same time variable, for now, just use 
                    # the first one.
                    try:
                        self.time_var = get_time_var(
                            xr.open_dataset(
                            self.filepaths[0], decode_times=False
                            ))
                    except:
                        self.time_var = None    # Search in datasets as files 
                                                # are added.

                    # TODO: verify all datasets in group have same time 
                    # variable & time points.
                    # if all([tv==time_var[0] for tv in time_var.values()]):
                    #     time_var = time_var[0]
                    # else:
                    #     #TODO: disambiguate.
                    #     raise ValueError(
                    #         f"Inconsistent time variables across dataset group {group.name}."
                    #         )
                else: 
                    self.time_var = None    # Search in datasets as files are 
                                            # added.

    def yield_datasets(self):
        return ds_feed(self.filepaths)
    
    def get_datasets(self):
        return ds_list(self.filepaths)
    
    def add(self, filepath: str) -> None:

        # global vocabulary  # Allow time variable names to be updated
        
        if not op.exists(filepath):
            raise OSError(f'DsGroup.add: Filepath {filepath} not found.')
        
        if op.isdir(filepath):
            for path in iglob(op.join(filepath, '*.nc'), 
                              os.pathsep, 
                              recursive=False):
                self.add(path)

        # Add filepath to collection
        self.filepaths.append(filepath)
        
        # Update name stem common to all files in group.
        self.stem = stem_str(self.stem, filepath)

        # print(f"Added {op.basename(filepath)} to {self.name}; stem: {self.stem}")

        # If don't have a time variable yet, attempt to find in new file.
        if not self.time_var:
            try:
                self.time_var = get_time_var(
                    xr.open_dataset(
                        filepath, 
                        decode_times=False
                    )
                )
            except KeyError:
                print(f'DsGroup.add: Time variable not found in {filepath}.')
        # If time variable contains wildcards, seek matching variable in ds. If
        # found, update vocabulary accordingly.
        elif len({'*', '?'}.intersection(self.time_var)) > 0:
            exact_time_var = variable_glob(
                xr.open_dataset(
                    filepath, 
                    decode_times=False
                ), self.time_var
            )
            if len(exact_time_var) != 1:
                # TODO: disambiguation
                raise AttributeError(
                    'Multiple variables found in dataset that match time '
                    f'variable pattern, {self.time_var}.')
            exact_time_var = exact_time_var[0]
            vocabulary[self.n_dims][exact_time_var] = vocabulary[self.n_dims][self.time_var]
            vocabulary[self.n_dims].pop(self.time_var)
            self.time_var = exact_time_var
        else:
            # TODO: confirm time variable in new dataset matches self.time_var
            pass

    def merge_times(self, shared: dict) -> None:
        '''
        Should print & return name of saved file(s)
        '''

        # from cfize_ds import cfize_dataset  # To avoid circular imports
        # global vocabulary
        if shared['verbose']: print(f"Process {os.getpid()}: Merging time series in group {self.n_dims}:{self.name}, containing datasets {[op.basename(path) for path in self.filepaths]}.")

        hold_attrs = {}  # Global attributes to be set once for merged group.
        # time_var = {}

        processing = self.get_datasets()  # list(self.datasets())  # Open all datasets in group.
        for i, ds in enumerate(processing):
            # Find dataset's time coordinate variable
            # time_var[i] = variable_glob(ds=ds, var_glob=self.time_var)
            # if len(time_var[i]) != 1:
            #     # TODO: disambiguation
            #     raise AttributeError(
            #         'Multiple variables found in dataset that match time '
            #         f'variable pattern, {self.time_var}.')
            # else:
            #     time_var[i] = time_var[i][0]

            # Assign any required global attributes to variables, 
            # with associated attributes. Remove global attribute.
            '''
            Assume these are only correct for last time-point in series.
            '''
            
            last_time_point = ds[self.time_var].data[-1]  # time_var[i]
            try:
                new_vars = globals_to_vars(
                            ds=ds, 
                            time_var=self.time_var, 
                            last_time_point=last_time_point
                )  # time_var[i]
            except Exception as e:
                raise e
            for name, array in new_vars.items():
                ds.attrs.pop(name)
                processing[i] = ds.assign({name: array})  # Need to assign to original dataset in list rather than placeholder ds.
                
            # Store largest value per group of any 'one per group' global 
            # attributes, as value to assign to merged dataset.
            for attr, attr_type in group_attrs.items():
                value = ds.attrs[attr]
                try:
                    value = attr_type(value)
                except TypeError as e:
                    if attr_type == datetime:
                        try:
                            value = datetime.fromisoformat(value)
                        except:
                            d, t = value.split()
                            d = [int(n) for n in re.split('[/-]', d)]
                            t = [int(n) for n in t.split(':')]
                            try:
                                value = datetime(
                                    year=d[2],
                                    month=d[1],
                                    day=d[0],
                                    hour=t[0],
                                    minute=t[1],
                                    second=t[2]
                                )
                            except ValueError:
                                value = datetime(
                                    year=d[0],
                                    month=d[1],
                                    day=d[2],
                                    hour=t[0],
                                    minute=t[1],
                                    second=t[2]
                                )
                    else:
                        raise e
                if attr in hold_attrs:
                    if value > hold_attrs[attr]:
                        hold_attrs[attr] = value
                else:
                    hold_attrs[attr] = value

        # Merge all datasets in group along time variable. 
        # TODO: Use chunking if large.
        self.to_process = xr.combine_by_coords(
            processing,
            combine_attrs='drop_conflicts',
            join='exact',
            coords=[self.time_var],
            data_vars='minimal'
        )  # Assigned to self.to_process, because this becomes dataset to process for CF compliance.
        
        # Close individual datasets, keeping only resultant.
        del processing

        # Apply 'once per group' attributes back to resulting dataset.
        # Convert any datetime objects to strings, so xarray.Dataset.to_netcdf 
        # can work with them. This can't be done until now, because finding the 
        # maximum requires the correct data type.
        for k, v in hold_attrs.items():
            if isinstance(v, (datetime)):
                hold_attrs[k] = v.isoformat(sep=' ')
        self.to_process.attrs.update(hold_attrs)

        # Derive title/filename from group's stem & dimension
        # merged.attrs.update({'title': self.stem.strip('_ ') + self.n_dims if f'{self.n_dims}d' not in self.stem else self.stem.strip('_ ')})
        self.to_process.attrs['title'] = self.stem.strip('_ ') + self.n_dims if f'{self.n_dims}d' not in self.stem else self.stem.strip('_ ')

    def merge_groups(self, shared: dict) -> None:
        """
        This function assumes the member subgroups have already been 'CFized'.
        That means any date(time) attributes will be expected in ISO format
        (yyyy-mm-dd[[T]hh:mm:ss])
        """
        
        if shared['verbose']: print(f"Process {os.getpid()}: merging groups - {self.name} with {self.n_dims} dimensions, containing {[op.basename(path) for path in self.filepaths]}.")
        

        if self.action != 'merge_groups':
            print(
                'merge_groups is not valid if group.action != "merge_groups".'
            )
            return
        
        hold_attrs = {}  # Global attributes to be set once for merged group.

        self.to_process = self.get_datasets()

        for i, ds in enumerate(self.to_process):
        # Store largest value per group of any 'one per group' global 
        # attributes, as value to assign to merged dataset.
            for attr, attr_type in group_attrs.items():
                value = ds.attrs[attr]
                try:
                    value = attr_type(value)
                except TypeError as e:
                    if attr_type == datetime:
                        try:
                            value = datetime.fromisoformat(value)
                        except:
                            d, t = value.split()
                            d = [int(n) for n in re.split('[/-]', d)]
                            t = [int(n) for n in t.split(':')]
                            value = datetime(
                                year=d[0],
                                month=d[1],
                                day=d[2],
                                hour=t[0],
                                minute=t[1],
                                second=t[2]
                            )
                    else:
                        raise e
                if attr in hold_attrs:
                    if value > hold_attrs[attr]:
                        hold_attrs[attr] = value
                else:
                    hold_attrs[attr] = value

        # xarray.Dataset.merge -> new dataset.
        merged = merge_dimensions(*self.to_process)
        
        # Close source datasets
        [ds.close() for ds in self.to_process]
        self.to_process = []

        # Apply 'once per group' attributes back to resulting dataset.
        # Convert any datetime objects to strings, so xarray.Dataset.to_netcdf 
        # can work with them. This can't be done until now, because finding the 
        # maximum requires the correct data type.
        for k, v in hold_attrs.items():
            if isinstance(v, (datetime)):
                hold_attrs[k] = v.isoformat(sep=' ')
        merged.attrs.update(hold_attrs)

        # Assign new title to new dataset
        merged.attrs['title'] = self.name
        
        # Set filepath and save as processed
        filepath = op.join(shared['target_dir'], f"{merged.attrs['title']}.nc")

        # Set encoding
        # TODO: encoding needs to be set, with each variable's encoding 
        # specified. Otherwise, _FillValue is applied to all, including 
        # coordinates, the latter contravening CF Conventions.
        encodings = {
            k:{
                'dtype': v.dtype,
                '_FillValue': None
            } for k, v in merged.variables.items()
        }  # if k == 'options_database' or k in ds.ds.coords

        # Save dataset to NetCDF
        # xarray docs report engine='h5netcdf' may sometimes be 
        # faster. However, it doesn't natively handle the 
        # various string length types used here.
        merged.to_netcdf(
            path=filepath, 
            encoding=encodings)  # , engine='h5netcdf'

        # TODO: Run each file through cf-checker?
        
        # Close dataset
        merged.close()

        self.processed = [filepath]
        return self.processed
        
    def cf_only(self, 
                shared: dict):
        '''
        TODO: This function still to be tested.
        '''
        if shared['verbose']: print(f"Processing {group.n_dims}:{group.name}.")

        # self.to_process = list(self.datasets())  # Open all datasets in group.
        self.cfize_and_save(shared=shared)
        
    
    def cfize_and_save(self, 
                       shared: dict
                       ):
        
        if shared['verbose']: print(f"Process {os.getpid()}: cfizing & saving datasets in group {self.name} with {self.n_dims} dimensions.")

        update_globals = {}

        target_dir = shared['target_dir']
        # time_units = shared['time_units']
        # vocab = shared['vocabulary']

        # global reference_vars  # Allow reference_vars to be updated.
        self.processed = []

        if self.to_process and isinstance(self.to_process, xr.Dataset):
            self.to_process = [self.to_process]
        else:
            self.to_process = self.get_datasets()

        for dataset in self.to_process:
            # Set up MoncDs object for further processing
            ds = MoncDs(
                dataset=dataset,
                time_variable=self.time_var,
                n_dims=self.n_dims,
                title=dataset.attrs['title']
            )
            
            # Call CF compliance function on dataset.
            cf_ds = ds.cfize(shared=shared)  # xarray.Dataset

            # Check whether time variable name has changed
            self.time_var = ds.time_var # Update group's time_var to match 
                                        # any update to merged dataset's.

            # Set filepath and save as processed
            filepath = op.join(target_dir, f"{cf_ds.attrs['title']}.nc")

            # Check for any reference variables for perturbation variables, and 
            # note filepath if found:
            # [
            #     reference_vars[v].update({'filepath':filepath})
            #     for v in reference_vars.keys()
            #     if v in cf_ds.variables
            # ]
            [
                shared['vocabulary'][reference_vars[v]['for'][0]][reference_vars[v]['for'][1]].update({'ref_array': cf_ds[v]})
                for v in reference_vars.keys()
                if v in cf_ds.variables
            ]
            for v in reference_vars.keys():
                if v in cf_ds.variables:
                    [perturbation_dim, 
                     perturbation_var] = reference_vars[v]['for']
                    add_to_var = {'ref_array': cf_ds[v]}
                    if 'vocabulary' in update_globals:
                        if perturbation_dim in update_globals['vocabulary']:
                            update_globals['vocabulary'][perturbation_dim].update({perturbation_var: add_to_var})
                        else:
                            update_globals['vocabulary'].update({
                                perturbation_dim: {perturbation_var: add_to_var}
                            })
                    else:
                        update_globals['vocabulary'] = {
                            perturbation_dim: {perturbation_var: add_to_var}
                        }
            
            # Set encoding
            # Encoding needs to be set, with each variable's encoding 
            # specified. Otherwise, _FillValue is applied to all, including 
            # coordinates, the latter contravening CF Conventions.
            # encodings = {
            #     k:{
            #         'dtype': v.dtype,
            #         '_FillValue': None
            #     } for k, v in cf_ds.variables.items()
            # }  # if k == 'options_database' or k in ds.ds.coords
            encodings = {
            k:{
                'dtype': v.dtype,
                '_FillValue': None,
                COMPRESSION[1]: COMPRESSION[0],
                'complevel': COMPRESSION[2]
            } for k, v in cf_ds.variables.items()
        } if COMPRESSION[0] else {
            k:{
                'dtype': v.dtype,
                '_FillValue': None
            } for k, v in cf_ds.variables.items()
        }  # if k == 'options_database' or k in ds.ds.coords

            # Save dataset to NetCDF
            # xarray docs report engine='h5netcdf' may sometimes be 
            # faster. However, it doesn't natively handle the 
            # various string length types used here.
            if COMPRESSION[0]:
                cf_ds.to_netcdf(
                    path=filepath, 
                    encoding=encodings,
                    engine='netcdf4',
                    format='netCDF4'
                )
            else:
                cf_ds.to_netcdf(
                    path=filepath, 
                    encoding=encodings)  # , engine='h5netcdf'

            # TODO: Run each file through cf-checker?
            
            # Close dataset
            cf_ds.close()

            self.processed.append(filepath)

        # if len(self.processed) == 1:
        #     self.processed = self.processed[0]

        return update_globals
    