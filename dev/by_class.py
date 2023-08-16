import xarray as xr
from setup import *
import os
from glob import glob
from cfunits import Units
from dataset_functions import *
from variable_functions import *
from difflib import SequenceMatcher
from utils import type_from_str, generate_coords, stem_str
import numpy as np
from datetime import datetime
from time import perf_counter
# import h5netcdf


class DirectoryParser:
    DIM_GROUPS = {
        0: '0+1d', 
        1: '0+1d', 
        2: '2d', 
        3: '3d'
    }

    DIM_ACTIONS = {
        '0+1d': 'merge',
        '2d': 'merge', 
        '3d': 'split'
    }

    def __init__(self,
                 directory: str = '',
                 target: str = '') -> None:
        '''
        directory: path of directory containing NC files to be made CF compliant.
        target (optional): path to which processed files should be written.
        '''

        # validate parameters
        if not all([isinstance(param, str) for param in (directory, target)]):
            raise TypeError('Source (and optionally target) directories must be passed as strings.')
        if not os.path.exists(directory):
            raise OSError(f'Directory {directory} not found.')

        self.directory = directory or os.getcwd()
        self.monc_files = []
        self.input_files = []
        self.variables = {}
        self.by_dim = None  # Purges any DsGroup instances from previous run/debug
        self.by_dim = {n_dim: DsGroup(
            name=group, action=self.DIM_ACTIONS[group]
        ) for n_dim, group in self.DIM_GROUPS.items()}
        self.sort_nc()
        self.time_units = self.time_units_from_input() if self.input_files else None
        self.target_dir = target or os.path.join(os.path.dirname(self.directory), f'{os.path.basename(self.directory)}+processed')
        if not os.path.exists(self.target_dir):
            try:
                os.makedirs(self.target_dir)
            except Exception as e:
                raise e
        self.processed = {group: None for group in self.DIM_ACTIONS}

    def sort_nc(self):
        '''
        Sorts into MONC output files and other NC files, the latter listed as
        possible input files.
        It also categorizes the MONC outputs according to their number of 
        dimensions.
        '''
        for file in glob('*.nc', root_dir=self.directory, recursive=False):
            filepath = os.path.join(self.directory, file)
            with xr.open_dataset(filepath, 
                                 decode_times=False) as ds:  # concat_characters=False
                # Using concat_characters=False preserves the string dimension 
                # of the options_database variable. However, it makes dealing
                # with the binary -> string conversion painful and messy.
                if is_monc(ds):
                    self.monc_files.append(filepath)
                    # Find number of (non-options-database) dimensions, and add
                    # filepath to the appropriate dataset group.
                    n_dims = get_n_dims(ds) - 1
                    self.by_dim[n_dims].add(filepath=filepath)
                    if len(self.by_dim[n_dims].filepaths) == 1:
                        # If first dataset in group, add variables to 
                        # dictionary linking each variable to its group.
                        self.variables.update({
                            v: self.by_dim[n_dims]
                            for v in ds.variables if v not in self.variables
                            })
                        self.variables.update({
                            v: self.by_dim[n_dims] 
                            for v in ds.dims if v not in self.variables
                            })
                else:
                    self.input_files.append(filepath)

    def time_units_from_input(self) -> Units:
        # TODO: this will use TimeUnits class
        '''Attempt to find specified input data in possible input file(s)'''
        input_data = {}
        for f in self.input_files:
            with xr.open_dataset(os.path.join(self.directory, f), decode_times=False) as ds:  # , concat_characters=False
                try:
                    input_data[f] = ds[INPUT_FILE['reftime']].attrs
                except KeyError:
                    # Variable not found in prospective input file; move on to
                    # next input file, if any.
                    continue
        if len(input_data) != 1:
            # TODO: work out how to decide between potentially conflicting data,
            # and/or request user input.
            if len(input_data) > 1:
                print('Multiple time units found:')
                while True:
                    user_input = input('Please select number or enter units for time, including reference date in ISO format.')
                    if user_input.isnumeric():
                        try:
                            input_data = list(input_data.values())[int(user_input)]
                        except:
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
            except ValueError as e:
                input_data['units'] = input(e + 'Enter units for time, including reference date in ISO format:')
                continue
            else:
                return time_units
    
    def datasets(self):
        # Generator
        return (xr.open_dataset(path, decode_times=False) for path in self.monc_files)  # , concat_characters=False
    
    def cfize(self):
        '''
        TODO: This function should probably be separate from DirectoryParser,
        as it should be possible for it to accept a group or individual 
        dataset & still work.
        '''
        for dim, group in self.by_dim.items():
            print(perf_counter(), 'Merging datasets in group', group.name)
            # First, merge time-points within dimension group, if required.
            # Any splitting into separate time points is done after other 
            # processing.
            if group.action == 'merge':
                # extract global attributes for each file, and assign each as a `data_var` with dimension of time.
                for ds in group.datasets:
                    ds.attr_to_var(CONFIG['global_to_variable'])

                # Merge all datasets in group to create single time-series dataset
                group.process()
        
        # Merge any groups with the same name: these comprise different
        # dimensions that should be merged, i.e. merging the 0d and 1d
        # time series into 0+1d.
        # Collect the merged groups, and any groups not needing to be merged,
        # into self.processed dict.
        for group in self.processed.keys():
            groups_to_merge = [g for g in self.by_dim.values() if g.name == group]
            if len(groups_to_merge) > 1:
                # Create new DsGroup containing groups to be merged
                self.processed[group] = DsGroup(groups=groups_to_merge)
                # Perform the merge
                print(f'{perf_counter()}: Merging groups {" & ".join([g.name for g in groups_to_merge])}')
                self.processed[group].process()
                # Update self.variables to point now to the merged group
                self.variables.update({
                    v: self.processed[group]
                    for v in self.processed[group].processed.variables
                    })
                self.variables.update({
                    v: self.processed[group]
                    for v in self.processed[group].processed.dims
                    })
            else:
                # If only one group per group name, point self.processed to that group for ongoing processing.
                self.processed[group] = groups_to_merge[0]

        # Update global attributes and variables. Split time-points if required.
        for name, group in self.processed.items():
            print(perf_counter(), 'Working on group', name)

            # combine filename stem & name from relevant group, to use as title
            title = group.stem + group.name if group.name not in group.stem else group.stem.strip('_ ,+')

            if group.action == 'split':
                # Process each dataset in group. Group won't have a processed attribute yet.
                to_process = group.datasets
            else:
                # Process group.processed only
                to_process = [MoncDataset(dataset=group.processed)]
                
            for ds in to_process:
                print(f'{ds.ds.attrs["title"]}:')
                # Convert required `options_database` items to global attrs.
                print(perf_counter(), 'Adding attributes from options database.')
                ds.options_to_attrs()  # This ds is the same object as to_process[0]

                # Add CF-required globals, using data from config file.
                print(perf_counter(), 'Adding attributes required by CF.')
                ds.add_cf_attrs(title=title)

                # add missing coordinate variables
                print(perf_counter(), 'Adding any missing coordinate variables.')
                ds.missing_coords()

                # Apply CF-required updates to all variables, using vocabulary.
                # Need to pass in self, to give access to variable dictionary,
                # if any variables might need converting from perturbations to
                # absolutes.
                print(perf_counter(), 'Making variables CF-compliant.')
                ds.cf_var(time_units=self.time_units, parser=self)
                if group.time_var != ds.time_var:
                    # group.time_var = ds.time_var
                    self.processed[name].time_var = ds.time_var  # This assigns to the original group, rather than the ephemeral reference, group.
            
                if group.action == 'split':
                    # Remove any attribute contained in ds.do_not_duplicate.

                    # Separate datasets by timepoint
                    print(perf_counter(), 'Splitting dataset by time-point.')
                    n_ds = group.split_times(ds)
                    
                    # Finish processing the new datasets created by split.
                    # TODO: this may be better off in the split_times function, or in between.
                    for i in range(len(group.processed) - n_ds, len(group.processed)):  #, new_ds in enumerate(group.processed[-n_ds:]):
                        new_ds = group.processed[i]
                        # update or drop `Previous diagnostic write at` & `MONC time`.
                        if 'Previous diagnostic write at' in ds.split_attrs:
                            if i == 0:
                                # First time point
                                prev_write = 0.
                            else:
                                prev_write = group.processed[i-1].time.data.tolist()
                            new_ds.attrs['Previous diagnostic write at'] = prev_write
                        
                        if 'MONC time' in ds.split_attrs:
                            new_ds.attrs['MONC time'] = new_ds.time.data.tolist()

                        # Update title
                        title = new_ds.attrs['title']
                        print(perf_counter(), 'Creating file from dataset', title)

                        if i < len(group.processed) - 1:
                            # Remove unknown attribute values for all but last 
                            # dataset.
                            for attr in ds.do_not_duplicate:
                                new_ds.attrs[attr] = np.nan

                        # Save each ds as NetCDF file, using required lossless 
                        # compression if specified.
                        filepath = os.path.join(
                            self.target_dir, 
                            title + '.nc'
                        )

                        # TODO: encoding needs to be set, with each variable's encoding specified. Otherwise, _FillValue is applied to all, including coordinates, the latter contravening CF Conventions.
                        encodings = {
                            k:{
                                'dtype': v.dtype,
                                '_FillValue': None
                            } for k, v in new_ds.variables.items()
                        }  # if k == 'options_database' or k in ds.ds.coords

                        # xarray docs report engine='h5netcdf' may sometimes be 
                        # faster. However, it doesn't natively handle the 
                        # various string length types used here.
                        new_ds.to_netcdf(
                            path=filepath, 
                            encoding=encodings)  # , engine='h5netcdf'

                        # TODO: Run each file through cf-checker


                        print(perf_counter(), 'Closing new dataset.')
                        new_ds.close()  # Free up some memory; will this interfere if managing write using dask?
                
                else:
                    print(perf_counter(), 'Creating file from dataset', title)
                    # Save resulting single dataset as NetCDF, using required lossless compression if specified.
                    # TODO: encoding attribute to include required compression settings, some of which require engine to be specified too.
                    # TODO: compute allows build & write operations to be delayed using dask.
                    filepath = os.path.join(
                        self.target_dir, 
                        ds.ds.title + '.nc'
                    )
                    # TODO: encoding needs to be set, with each variable's encoding specified. Otherwise, _FillValue is applied to all, including coordinates, the latter contravening CF Conventions.
                    encodings = {
                        k:{
                            'dtype': v.dtype,
                            '_FillValue': None
                        } for k, v in ds.ds.variables.items()
                    }  # if k == 'options_database' or k in ds.ds.coords
                    
                    ds.ds.to_netcdf(path=filepath, encoding=encodings)

                    # TODO: Run each file through cf-checker
                    # This would be easiest if import the checker itself & can pass dataset or file to it.
                    # CFChecker.checker(filepath)

                    # Update group.processed with current version of dataset
                    group.processed = ds.ds  # Expects xr.Dataset
        
        print(perf_counter(), 'Done processing groups.')

    def get_ref_var(self, variable: str) -> xr.DataArray:
        if variable not in self.variables:
            raise KeyError(f'{variable} not found in datasets processed.')
        ref_group = self.variables[variable]
        # ref_group = self.processed[ref_group_name]
        ref_ds = ref_group.processed
        return ref_ds[variable]


class TimeUnits(Units):
    '''
    Uses proleptic_gregorian as default calendar, as per ISO 8601:2004,
    even though CF standard calendar is a mixed Gregorian/Julian.
    '''
    def __init__(self, units=None, calendar='proleptic_gregorian', formatted=False, names=False, definition=False, _ut_unit=None):
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
        return self.time_unit() == other.time_unit()
    
    def base_units_equivalent(self, other: Units) -> bool:
        return Units(units=self.time_unit()).equivalent(Units(units=other.time_unit()))
    
    def ref_dates_match(self, other: Units) -> bool:
        return self.since() == other.since()
    
    def has_calendar(self) -> bool:
        try:
            return self.calendar is not None
        except AttributeError:
            return False

    def calendars_match(self, other: Units) -> bool:
        if all([u.has_calendar for u in [self, other]]):
            return self.calendar == other.calendar
        else:
            return False
    

class MoncDataset:
    ''''
    Notes:
        This has been set up as a standalone class, rather than inheriting from xr.Dataset, to avoid some complexities of inheriting from xarray data classes.
    '''
    do_not_duplicate = {'MONC timestep'}
    split_attrs = {
        'title',
        'MONC time',
        'Previous diagnostic write at'
    }
    
    def __init__(self, filepath: str = None, dataset: xr.Dataset = None):
        if dataset:
            if not isinstance(dataset, xr.Dataset):
                raise TypeError
            self.ds = dataset if is_monc(dataset) else None
        elif filepath:
            if not isinstance(filepath, str):
                raise TypeError('filepath must be a string')
            try:
                self.ds = xr.open_dataset(filepath, decode_times=False)  # , concat_characters=False
            except Exception as e:
                raise e
            else:
                if not is_monc(self.ds):
                    self.ds = None
        else:
            raise ValueError('Either a filepath (string) or xarray.Dataset must be passed in.')
        if self.ds:
            self.options = None
            self.set_options(fields=CONFIG['options_to_attrs'])
            time_vars = [d for d in self.ds.dims if 'time' in d]
            if len(time_vars) != 1:
                # Throw error or prompt user to select which dimension to use
                raise AttributeError(f'More than one time dimension found: {time_vars}')
            self.time_var = time_vars[0]
    
    def set_options(self, fields: list|set|tuple = None):
        # Get required fields from argument; import all if none supplied.
        if fields is None:
            if self.ds.options_database.dtype == 'S1':
                # This is used if concat_characters=False options is used in
                # xarray.open_dataset. It seems to create some problems down
                # the line, so is not advised.
                self.options = {
                    ''.join([c.decode('utf-8') for c in k]): 
                    type_from_str(''.join([c.decode('utf-8') for c in v]))
                    for k, v in self.ds.options_database.data
                }  # If dataset opened with concat_characters=False
            else:
                self.options = {
                    k.decode('utf-8'): type_from_str(v.decode('utf-8'))
                    for [k, v] in self.ds.options_database.data}
        else:
            if self.ds.options_database.dtype == 'S1':
                self.options = {
                    ''.join([c.decode('utf-8') for c in k]): 
                    type_from_str(''.join([c.decode('utf-8') for c in v]))
                    for k, v in self.ds.options_database.data
                    if ''.join([c.decode('utf-8') for c in k]) in fields
                }  # If dataset opened with concat_characters=False
            else:
                self.options = {
                    k.decode('utf-8'): type_from_str(v.decode('utf-8'))
                    for [k, v] in self.ds.options_database.data
                    if k.decode('utf-8') in fields}
        
        # Drop any options that are reset by DEPHY, if used.
        if (all([opt in self.options for opt in DEPHY_OPTIONS]) and
            all([self.options[opt] for opt in DEPHY_OPTIONS])):
            for a in DROP_FOR_DEPHY:
                if a in self.options:
                    self.options.pop(a)
        
    def options_to_attrs(self):
        self.ds = self.ds.assign_attrs(**self.options)
        
    def add_cf_attrs(self, **kwargs):
        # TODO: <version> and <url> and any other required data will be assigned during packaging.
        defaults = {attr: None for attr in CF_ATTRIBUTES}
        defaults['title'] = kwargs['title'] if 'title' in kwargs and kwargs['title'] else self.title
        defaults['history'] = f'{datetime.now().isoformat(timespec="minutes")}: output files processed using CFizer version <version>, <url>.'
        defaults['conventions'] = 'CF-1.10'
        for attr in CF_ATTRIBUTES:
            self.ds.attrs[attr] = CONFIG[attr] if attr in CONFIG else defaults[attr]


    def attr_to_var(self, attr: str|list|set|tuple = None):
        
        if attr is None:
            # convert all attributes to variables
            attr = set(self.ds.attrs.keys())
        
        if isinstance(attr, str):
            attr = {attr}
        
        last_time_point = self.ds[self.time_var].data[-1]
        globals = {}
        if all([isinstance(a, str) for a in attr]):
            for g in attr:
                name = g.replace(' ', '_')  # TODO: Replace this with a full "safe_string" function
                if g in self.do_not_duplicate:  # name
                    # TODO: this should assign `np.nan` to each data point except last, and then set coords as per other case below.
                    data = type_from_str(self.ds.attrs[g])
                    coords = [last_time_point]
                else:
                    data = [type_from_str(self.ds.attrs[g])]*len(self.ds[self.time_var].data)
                    coords = self.ds[self.time_var]
                globals[g] = xr.DataArray(
                    name=name,
                    data=data,
                    coords={self.time_var: coords},
                    dims=(self.time_var,)
                )  # TODO: Test reverting to globals[name] = ...
        else:
            raise TypeError('attr must be a string, a collection of strings, or None.')
        self.ds = self.ds.assign(
            {name: array for name, array in globals.items()}
        )

    def missing_coords(self):
        '''
        Look for any dimensions that currently don't also exist as coordinate 
        variables.
        '''
        dim_type = np.dtype('float32')
        missing = []
        for dim in self.ds.dims:
            # Set type to match existing spatial coordinates (assume for now all are the same)
            if dim in self.ds.coords:
                if 'time' not in dim:
                    dim_type = self.ds.coords[dim].dtype
            elif dim not in OPTIONS_DATABASE['dimensions']:
                missing.append(dim)
        
        # Look for any coordinates listed in config file that weren't already found in missing coordinates
        # [missing.append(dim) for dim in CONFIG['new_coordinate_variables'].keys() if dim not in missing]  # Remove this: it causes coordinates to be added when they are not present in any variables' dimensions
        
        for dim in missing:
            attributes = {}
            if dim in CONFIG['new_coordinate_variables']:
                # Get required information from config file, if available
                config = CONFIG['new_coordinate_variables'][dim]
                spacing = self.options[config['spacing']]
                midpoint = 'cent' in config['position'] or 'mid' in config['position']
                for k, v in config['attributes'].items():
                    if k.lower() == 'units':
                        attributes['units'] = Units(v).formatted()
                    else:
                        attributes[k.lower().replace(' ', '_')] = v if k != 'standard_name' else v.replace(' ', '_')
            else:
                # If not in config file, assume it's x, y, xu or yu, with associated default attributes etc.
                if dim[0] == 'x':
                    spacing = self.options['dxx']
                    attributes = {
                        'axis': 'X',
                        'units': 'm'
                    }
                elif dim[0] == 'y':
                    spacing = self.options['dyy']
                    attributes = {
                        'axis': 'Y',
                        'units': 'm'
                    }
                else:
                    # TODO: prompt user for grid spacing
                    raise AttributeError(f'Unknown grid spacing for {dim}.')
                attributes['units'] = 'm'
                # 'xu' and 'yv' are at endpoints, by default. 'x' and 'y' at 
                # cell centres.
                midpoint = dim[-1] not in ('u', 'v')
                attributes['long_name'] = f'{dim[0]}-coordinate in Cartesian system (cell-{"centers" if midpoint else "edges"})'
                
            # Generate data points for coordinate variable
            # First convert spacing to required np.dtype
            # spacing = np.array(spacing).astype(dim_type).tolist()
            points = generate_coords(number=self.ds.dims[dim],
                                     spacing=spacing,
                                     midpoint=midpoint,
                                     data_type=dim_type)
            
            # create new xarray variable based on the new coordinate
            new_var = xr.Variable(dims=dim, attrs=attributes, data=points)

            # Add new coordinate variable to dataset
            self.ds = self.ds.assign(variables={dim: new_var})

    def cf_var(self, variable: str = None, time_units: TimeUnits = None, parser: DirectoryParser = None) -> None:
        # If no variable specified, attempt to process all variables present.
        if not variable:
            [self.cf_var(variable=var, time_units=time_units, parser=parser) for var in self.ds.variables.keys() if var != OPTIONS_DATABASE['variable']]
            return

        # Check variable exists in self.ds.
        if not variable in self.ds.variables.keys():
            raise KeyError(f'{variable} not a variable of dataset.')

        # Find variable in vocabulary
        try:
            updates = VOCABULARY[variable]
        except KeyError:
            # Check whether variable is already CF compliant. Skip rest of function if so.
            if cf_compliant(self.ds[variable]):
                return
            
            # TODO: attempt to find in `standard_names` table
            if 'standard_name' in self.ds[variable].attrs:
                pass
            
            # Temporarily allow program to continue
            return # raise KeyError(f'{variable} not found in vocabulary.')

        '''Apply changes to variable where necessary'''
        # TODO: split these into separate functions

        # Update dims/coordinates
        if 'dimension_changes' in updates:
            # Don't need to check new_dim exists, if missing_coords() is run
            # afterwards. But because current workflow has missing_coords before
            # updating variables, do the check anyway.
            if not set(updates['dimension_changes'].keys()).issubset(set(self.ds[variable].dims)):
                raise KeyError(
                        f'{set(updates["dimension_changes"].keys())} are not all dimensions of {variable}.'
                    )
            # Apply change(s)
            self.ds[variable] = self.ds[variable].swap_dims(updates['dimension_changes'])
            # Check new_dim exists
            if not set(self.ds[variable].dims).issubset(self.ds.coords):
                self.missing_coords()
                
        # Update/add name attributes
        if 'standard_name' in updates:
            # TODO: check standard_name is in list of standard names.

            self.ds[variable].attrs['standard_name'] = updates['standard_name']

        if 'long_name' in updates:
            self.ds[variable].attrs['long_name'] = updates['long_name']

        # For CF compliance, at least one name attribute is required
        if 'standard_name' not in self.ds[variable].attrs and 'long_name' not in self.ds[variable].attrs:
            raise KeyError('CF requires at least one of standard_name and long_name to be assigned to each variable.')
        
        # Check units cf-compliant & assign
        if 'units' in updates:
            if 'standard_name' in self.ds[variable].attrs:
                # TODO: check units supplied are consistent with those in standard_name table
                pass

            # If time, derive corrected units
            if all([
                ('time' in variable or (
                'axis' in self.ds[variable].attrs and 
                self.ds[variable].attrs['axis'] == 'T')
                ),
                TimeUnits(units=updates['units']).base_units_equivalent(TimeUnits(units='s')),
                variable in self.ds.coords]):
                    
                # Check provided unit against that obtained from input file or user.
                new_unit = TimeUnits(units=updates['units'])
                if not new_unit.equals(time_units):
                    # Can we derive workable combination?
                    if new_unit.isreftime:
                        if new_unit.iscalendartime:
                            # Incompatible; either ask for user confirmation
                            # or use units from vocabulary.
                            pass
                        else:
                            new_unit = TimeUnits(units=new_unit.units, calendar=time_units.calendar)
                    else:
                        if time_units.base_units_equivalent(other=new_unit):
                            new_unit = time_units
                        else:
                            # TODO: query user, with time_units.units and time_units.calendar as prompts
                            new_unit = TimeUnits(
                                units=f'{new_unit.units} since {time_units.since}',
                                calendar=time_units.calendar
                            )
                if new_unit.has_calendar():
                    self.ds[variable].attrs['calendar'] = new_unit.calendar
                
                # Assign axis attribute
                self.ds[variable].attrs['axis'] = updates['axis'] if 'axis' in updates else 'T'

            else:
                new_unit = Units(units=updates['units'])
                if not new_unit.isvalid:
                    # TODO: prompt user for valid value
                    raise ValueError(f'Invalid unit specified in vocabulary: {variable}: {updates["units"]}')
                
            # If unit is different from current, and it's equivalent (e.g. km -> m), apply conversion as well as new unit.
            if 'units' in self.ds[variable].attrs:
                current_unit = Units(self.ds[variable].attrs['units'])
                if new_unit.equivalent(current_unit) and not new_unit.equals(current_unit):
                    Units.conform(
                        x=self.ds[variable].data, 
                        from_units=current_unit, 
                        to_units=new_unit, 
                        inplace=True
                    )

            # Apply new units
            self.ds[variable].attrs['units'] = new_unit.formatted()

        elif 'units' not in self.ds[variable].attrs:
            # Look for unit in standard_names table

            raise AttributeError(f'units attribute missing for {variable}.')
        
        else:
            # Check existing unit is valid
            if 'units' not in self.ds[variable].attrs:
                raise AttributeError(f'units attribute must be assigned for each variable. None found for {variable}.')
            current_unit = Units(self.ds[variable].units)
            if current_unit.isvalid:
                self.ds[variable].attrs['units'] = current_unit.formatted()
            else:
                raise ValueError(f'Current unit of {variable}, {current_unit}, is not valid.')

        # Assign axis attribute if coordinate variable
        if variable in self.ds.coords:
            if 'axis' in updates:
                self.ds[variable].attrs['axis'] = updates['axis']
            elif 'axis' not in self.ds[variable].attrs:
                # Infer or query user
                if variable[0].lower() in {'x', 'y', 'z'}:
                    if input(f'Apply axis = {variable[0].upper()} to {variable}? (y/n)')[0].lower() == 'y':
                        self.ds[variable].attrs['axis'] = variable[0].upper()
                if any(['time' in x for x in {
                    variable, 
                    self.ds[variable].attrs['standard_name'],
                    self.ds[variable].attrs['long_name']
                    }]):
                    if input(f'Apply axis = T to {variable}? (y/n)')[0].lower() == 'y':
                        self.ds[variable].attrs['axis'] = 'T'
            # TODO: Ideally, this would check whether `positive` attribute is also needed. If provided, it would ideally check it is consistent with e.g. standard_name (standard_name = altitude, positive = up)
            if 'positive' in updates:
                self.ds[variable].attrs['positive'] = updates['positive']

        
        # Update name; this should be done last, so the existing name can still
        # be used as a key.
        if 'updated_name' in updates:
            if variable == self.time_var:
                self.time_var = updates['updated_name']
            self.ds = self.ds.rename({variable: updates['updated_name']})

        # perturbation to absolute
        if 'perturbation_to_absolute' in updates and updates['perturbation_to_absolute']:
            # This will require the reference variable to be brought in from a different dataset, in most if not all cases.
            # Any shared dimensions must have the same names, so this requires all processing to be done in ascending order of dimensions.

            if 'reference_variable' not in updates:
                raise KeyError(f'{variable}: If perturbation_to_absolute is True, reference_variable must contain the name of the variable containing reference value(s).')
            if updates['reference_variable'] not in self.ds.variables:
                # Look for reference variable in other dimension groups
                if parser is not None:
                    try:
                        ref_var = parser.get_ref_var(updates['reference_variable'])
                    except AttributeError as e:
                        raise AttributeError(f'{variable}: {e}')
                else: 
                    raise KeyError(
                        f'{variable}: Reference variable {updates["reference_variable"]} does not exist in dataset. If it is in another dataset, parser must be specified, to give access to the DirectoryParser.variables dictionary.'
                    )
            else:
                ref_var = self.ds[updates['reference_variable']]
            absolute = self.ds[variable] + ref_var
            self.ds[variable] = absolute


class DsGroup:
    '''
    Each file in a DsGroup has the same variables, coordinates & dimensions.
    '''
    def __init__(self, 
                 name: str = '', 
                 action: str = None, 
                 filepaths: list|set|tuple = None,
                 groups: list|set|tuple = None) -> None:
        
        # TODO: validate parameters

        if groups and all([isinstance(g, DsGroup) for g in groups]):
            # If group is to comprise multiple existing groups, check each 
            # has only one dataset in its processed attribute.
            if not all([isinstance(group.processed, xr.Dataset) for group in groups]):
                raise TypeError('To merge DsGroup objects, the processed attribute of each must contain a single xarray.Dataset object.')
            # merge existing groups
            self.name = stem_str(*[group.name for group in groups])
            self.time_var = groups[0].time_var
            if not all([g.time_var == self.time_var for g in groups[1:]]):
                raise ValueError('When attempting to combine groups, found mismatch in time variables.')
            self.action = 'merge_groups'
            self.process = self.merge_groups
            self.stem = stem_str(*[group.stem for group in groups])
            self.filepaths = []
            for group in groups:
                self.filepaths += group.filepaths
            self.to_merge = [group.processed for group in groups]
            self.processed = None
        else:
            self.name = name
            self.action = action.lower()
            if self.action == 'split':
                self.process = self.split_times
            elif self.action == 'merge':
                self.process = self.merge_times
            else:
                self.process = self.keep_time_series
            self.filepaths = filepaths if filepaths else []  # Cannot set default valueas an empty list  in arguments, or it assigns SAME empty list to every instance.
            self.processed = []
            self.stem = stem_str(*self.filepaths) if self.filepaths else None
            # Although it would be best to verify that all member datasets have the same time variable, for now, just use the first one.
            self.time_var = MoncDataset(filepath=self.filepaths[0]).time_var if self.filepaths else None
            # Creating a generator for datasets is not appropriate, because it
            # newly creates each MoncDataset every time it is used, meaning any
            # changes made to a given MoncDataset object are lost afterwards.
            self.datasets = [MoncDataset(path) for path in self.filepaths]

    # def datasets(self):
    #     '''
    #     This returns a generator that yields the datasets contained in each
    #     filepath.
    #     '''
    #     # return (xr.open_dataset(path, decode_times=False) for path in self.filepaths)
    #     return (MoncDataset(path) for path in self.filepaths)

    def add(self, filepath):
        # TODO: validate filepath
        self.filepaths.append(filepath)
        self.stem = stem_str(filepath, self.stem)
        self.datasets.append(MoncDataset(filepath))
        if not self.time_var:
            self.time_var = self.datasets[-1].time_var

    def keep_time_series(self):
        # copy all datasets/files into self.processed
        self.processed = [ds for ds in self.datasets]

    def merge_times(self):
        '''
        Using the option, data_vars='minimal', and merging only on
        coords=[self.time_var] ensures the options_database variable is not
        replicated for every time-point.
        combine_attrs='drop_conflicts' can be used because any necessary global
        attributes can be preserved by MoncDataset.attr_to_var().
        '''
        # Take latest 'created' time-point as the value for the combined ds.
        created = max(
            [monc_ds.ds.attrs['created'] for monc_ds in self.datasets]
        )
        self.processed = xr.combine_by_coords(
            [monc_ds.ds for monc_ds in self.datasets],
            combine_attrs='drop_conflicts',
            join='exact',
            coords=[self.time_var],
            data_vars='minimal'
        )
        self.processed.attrs['created'] = np.str_(created)

    def split_times(self, dataset: MoncDataset) -> int:
        # Because xarray.Dataset.sel only creates datasets that point to the
        # original, any changes to global attributes etc in one dataset will
        # propagate to all. To create independent datasets, each resulting ds
        # must be deep-copied.
        # split = None
        # times = [t for t in dataset.ds[self.time_var].data]
        # split = [dataset.ds.sel({self.time_var: t}).copy(deep=True) 
        #          for t in times]
        grouped = {t: ds for (t, ds) in dataset.ds.groupby(self.time_var)}
        for k, v in grouped.items():
            v.attrs['title'] = v.attrs['title'] + '_' + str(int(k))
        split = list(grouped.values())
        if self.processed is not None:
            self.processed += split
        else:
            self.processed = split
        dataset.ds.close()  # Original no longer needed
        return len(split)

    def merge_groups(self):
        self.processed = merge_dimensions(*self.to_merge)
        self.processed.attrs['title'] = self.name
        

def test():
    data_dir = os.path.join(os.path.dirname(app_dir), 'test_data')
    
    # print("Testing MoncDataset instantiation: passing in filepath")
    # filename = 'd20200128_diagnostic_1d_3600.nc'
    # print(f"Is {filename} MONC output? {MoncDataset(filepath=os.path.join(data_dir, filename)).is_monc()}")

    # print("Testing MoncDataset instantiation: passing in xarray.Dataset")
    # filename = 'd20200128_diagnostic_1d_3600.nc'
    # ds = xr.open_dataset(os.path.join(data_dir, filename), decode_times=False)
    # print(f"Is {filename} MONC output? {MoncDataset(dataset=ds).is_monc()}")

    print('Testing DirectoryParser initialisation.')
    parser = DirectoryParser(directory=data_dir)
    parser.cfize()
    


if __name__ == '__main__': test()