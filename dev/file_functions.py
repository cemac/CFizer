import os
from glob import glob, iglob
import re
import xarray as xr
from time import perf_counter
from difflib import SequenceMatcher
from setup import CONFIG, VOCABULARY, app_dir
import utils
from cfunits import Units
from datetime import datetime, timedelta
from cf_functions import *


DIM_GROUPS = ['0+1d', '0+1d', '2d', '3d']


def performance_time(func):
	def wrapper(*args, **kwargs):
		start_time = perf_counter()
		response = func(*args, **kwargs)  # run the wrapped function
		end_time = perf_counter()
		duration = end_time - start_time
		print(f"{func}({args}, {kwargs}) took {duration} seconds.")
		return response
	return wrapper


class MONC_Parser:
    def __init__(self, directory: str, dim_bins: list|tuple|set = None):
        self.source_dir = os.path.abspath(directory)
        
        # Check directory is valid
        if not os.path.exists(self.source_dir):
            raise OSError(f'Directory {self.source_dir} not found')
        
        # List of all NC files in target directory.
        # Don't scan subdirectories, as we want to treat each directory
        # as a collection.
        # If there's a risk the list of files will be very large, better
        # to make this a generator, using glob.iglob instead of glob.glob.
        self.nc_files = glob(
                pathname='*.nc',
                root_dir=directory,
                recursive=False
        )

        # Can't find shared string in filenames at this point, as we
        # don't know which files are MONC outputs.
        # self.stem_name = nc_basename(self.nc_files)

        # Target directory for processed files should be sibling of source directory,
        # so it doesn't interfere with scanning for files to be processed.
        # Ideally allow user to specify where processed files should go.
        self.target_dir =  os.path.join(
            os.path.dirname(directory),
            f'{os.path.basename(directory)}+processed'
        )
        if not os.path.exists(self.target_dir):
            os.makedirs(self.target_dir)

        # Separate MONC outputs from other NC files:
        self.non_monc = []
        i = 0
        while i < len(self.nc_files):
            if not is_monc_file(
                os.path.join(self.source_dir, self.nc_files[i])
            ):
                self.non_monc.append(self.nc_files.pop(i))
            else:
                i += 1

        if dim_bins is not None:
            if not isinstance(dim_bins, (list, tuple, set)):
                raise TypeError('dim_bins must be a sequence containing category names based on number of dimensions (not counting time dimension).')
            self.dim_bins = dim_bins
            # Set up dictionary to contain filepaths classified by number of dimensions
            self.by_dim = {}
            # Set up dictionary to contain filepaths of output files.
            self.stem = {name: None for name in set(self.dim_bins)}
            try:
                self.sort_by_dims()
            except ValueError as e:
                raise e

            
        # # Set up dictionary to contain filepaths, categorized by number of dimensions in data.
        # self.source = {v: [] for v in set(DIM_GROUPS)}
        # #     '0+1d': [],
        # #     '2d': [],
        # #     '3d': []
        # # }
        
        self._time_units = None
        # self._reference_datetime = []

        if CONFIG['missing_data']['time_units']['use_input_file']:
            # Attempt to obtain missing time unit data from input file(s)
            time_units_from_input = self.get_input_data(CONFIG['missing_data']['time_units']['input_file_variable'])
            if len(time_units_from_input.keys()) > 1:
                # TODO: Work out how to pick one set of data
                time_units_from_input = list(time_units_from_input.values())[0]
            else:
                time_units_from_input = list(time_units_from_input.values())[0]
            if 'units' not in time_units_from_input:
                # Use fall-back derivation of time units.
                pass
            if 'calendar' not in time_units_from_input:
                time_units_from_input['calendar'] = None
            try:
                self.set_time_units(units=time_units_from_input['units'],
                                    calendar=time_units_from_input['calendar'])
            except Exception as e:
                raise e
        elif CONFIG['missing_data']['time_units']['calendar'] is not None:
            try:
                self.set_time_units(calendar=CONFIG['missing_data']['time_units']['calendar'])
            except Exception as e:
                raise e
        
    def get_time_units(self):
        return None if self._time_units is None else self._time_units.formatted()
    
    def get_calendar(self):
        try:
            return self._time_units.calendar
        except AttributeError:
            return None
        
    def get_reference_date(self):
        try:
            return self._time_units.reftime
        except AttributeError:
            return None
    
    def set_time_units(self, 
                   units: Units|str = None, 
                   since: str|datetime = None, 
                   calendar: str|Units = None):
        if not isinstance(units, Units):
            if since is not None:
                if isinstance(since, datetime):
                    since = since.isoformat()
                elif not isinstance(since, str):
                    raise TypeError('since parameter must be string or datetime object.')
                units = f'{units} since {since}'
            if calendar is not None:
                if isinstance(calendar, Units):
                    calendar = calendar.calendar
                elif not isinstance(calendar, str):
                    raise TypeError('calendar parameter must be string or cfunits.Unit object.')
            else:
                calendar = CONFIG['missing_data']['time_units']['calendar']
            try:
                units = Units(units=units, calendar=calendar)
            except Exception as e:
                raise e
        if not units.isvalid:
            raise ValueError('Units are not valid according to CF convention.')
        self._time_units = units

    def sort_by_dims(self):
        for filename, filepath in [(f, os.path.join(self.source_dir, f)) for f in self.nc_files]:
            with xr.open_dataset(filepath, decode_times=False) as ds:
                # Need to set decode_times=False to pull out time units as defined.
                
                # Find number of dimensions, which determines how file is merged/split
                n_dims = get_n_dims(ds) - 1
                if n_dims in self.by_dim:
                    self.by_dim[n_dims].append(filepath)
                else:
                    self.by_dim[n_dims] = [filepath]

                # Find longest common string in files of same dimension classification (as per dim_bins).
                self.stem[self.dim_bins[n_dims]] = nc_basename(filename, self.stem[self.dim_bins[n_dims]])

        if max(list(self.by_dim.keys())) > len(self.dim_bins):
                raise ValueError('Not enough dimension categories specified.')
                        
    def get_input_data(self, var: str) -> dict:
        '''Attempt to find specified input data in possible input file(s)'''
        input_data = {}
        for f in self.non_monc:
            with xr.open_dataset(os.path.join(self.source_dir, f), decode_times=False) as ds:
                try:
                    input_data[f] = process_input(ds, var)
                except TypeError as e:
                    raise e
                except KeyError:
                    # Variable not found in prospective input file
                    continue
        return input_data

    def process(self, n_dims: int):
        dim_group = self.dim_bins[n_dims]
        if n_dims == 3:
            # split into individual time points
            nc_split(self.by_dim[3])
        else:
            # merge time points
            merged = nc_merge(self.by_dim[n_dims], globals_to_variables=CONFIG['global_to_variable'])


def nc_merge(datasets: list, 
             globals_to_variables: list|set|tuple) -> xr.Dataset:
    # TODO: validate datasets & globals_to_variables
    
    for ds in datasets:
        # Identify time variable
        time_dim = [d for d in ds.dims.keys() if 'time' in d][0]

        # Convert necessary global attributes to variables
        for g in globals_to_variables:
            name = g.replace(' ', '_')  # Replace this with a full "safe_string" function
            data = utils.type_from_str(ds.attrs[g])
            globals[g] = xr.DataArray(
                name=name,
                data=data,
                coords={time_dim: [data]},
                dims=[time_dim]
            )
        ds.assign({name: array for name, array in globals.items()})


def nc_split(ds: xr.Dataset) -> tuple:
    # TODO: validate ds

    pass


def nc_iterator(directory: str) -> None:
    '''
    Raises OSError if directory not found.
    '''
    
    try:
        nc_dir = MONC_Parser(directory)
    except OSError as e:
        raise e

    # # Iterate over all nc files in directory and any subdirectories
    # for nc_file in iglob(pathname='**/*.nc', root_dir=directory, recursive=True):
    for f, p in [(p, os.path.join(nc_dir.source_dir, p)) for p in nc_dir.nc_files]:
        
        # Does the file contain MONC output?
        with xr.open_dataset(p, decode_times=False) as ds:
            # Need to set decode_times=False to pull out time units as defined.
            if is_monc(ds):
                print(f'{f} is MONC output. Processing...')
                
                # Find number of dimensions, which determines how file is merged/split
                n_dims = get_n_dims(ds) - 1
                nc_dir.source[DIM_GROUPS[n_dims]].append(f)
                nc_dir.stem[DIM_GROUPS[n_dims]] = nc_basename(
                    nc_dir.stem[DIM_GROUPS[n_dims]], f
                )
                print(f'{n_dims}d filename stem: {nc_dir.stem[DIM_GROUPS[n_dims]]}')
                # if n_dims < 3:
                #     # Merge with other "0d" and "1d" files (time alone or time+altitude)
                #     nc_dir.source['0+1d'].append(f)
                #     if nc_dir.target['0+1d'] is None:
                #         # # Derive name for target file from current filename
                #         # nc_dir.target['0+1d'] = nc_basename(f)
                        
                #         # Temporarily set target filename to match current filename,
                #         # as cannot determine the common stem until a 2nd file is
                #         # accessed.
                #         nc_dir.target['0+1d'] = f

                #         # Clone current dataset, to become first part of merged dataset.
                #         # merged = ds.copy(deep=True)  # deep option makes a full copy of all data
                        
                #         # Write merged dataset as new NC file.

                #     else:
                #         nc_dir.target['0+1d'] = nc_basename(nc_dir.target['0+1d'], f)
                        
                #         # merge_time_series(new_ds=ds, merge_file=nc_dir.target['0+1d'])

                #     print(f'Merged file name: {nc_dir.target["0+1d"]}')

                # elif n_dims < 4:
                #     # Merge with other "2d" files (time + x + y)
                #     pass
                #     # merge_time_series(new_ds=ds, merged_ds=merged_2d)
                # else:
                #     # Split "3d" (time + 3 spatial dimensions) into individual time points.
                #     # single_timepoints = split_time_series(ds)
                #     # for tp in single_timepoints:
                #     #     # write new NetCDF file to target directory
                #     #     pass
                #     pass

                # process_monc(ds)

                reference_datetime = {k.decode('utf-8').lower().replace('rad_start_', ''): v.decode('utf-8') for [k, v] in ds[CONFIG['options_database']['variable']].data if k.decode('utf-8') in CONFIG['options_database']['reference_datetime']}
                for k, v in reference_datetime.items():
                    try:
                        reference_datetime[k] = int(float(v)) if int(float(v)) == float(v) else float(v)  # int(v) won't work e.g. for v = '2020.0'
                    except ValueError:
                        reference_datetime = None
                reference_datetime = datetime(year=reference_datetime['year'], month=1, day=1) + timedelta(days=(reference_datetime['day'] - 1), seconds=reference_datetime['time'])
                if reference_datetime not in nc_dir.reference_datetime:
                    nc_dir.reference_datetime.append(reference_datetime)

            else:
                print(f'{f} is not MONC output.')
                # Check whether file contains useful input data
                monc_input = process_input(ds)
                monc_input['input file'] = f
                

    if CONFIG['time_units']['use_input_file']:
        time_unit = Units(monc_input['time units'])
        calendar = Units(calendar=(monc_input['time calendar'] or CONFIG['time_units']['calendar'] or 'standard'))
        if (calendar.isvalid and calendar.isreftime):
            calendar = calendar.calendar
        else:
            calendar = 'standard'
        if time_unit.isreftime:
            time_unit = time_unit.formatted()
        else:
            # Look for origin date (reference datetime) in rad_start
            # entries in options_database.
            # Problem is that this would ideally be done on first pass
            # of files, when we're only categorising by dimension.
            # It suggests we need to pull time units from options_database of
            # each file during this first pass.
            pass
        nc_dir.time_units = Units(units=time_unit,
                                  calendar=calendar)
    
    if nc_dir.time_units is None:
        # Set from other config file options
        if 'options_database' in str(CONFIG['time_units']['ref_datetime']):
            if len(nc_dir.reference_datetime) == 1:
                time_unit = f"{VOCABULARY['time']['correct units']} since {nc_dir.reference_datetime[0].isoformat()}"
            else:
                #  Prompt user to select/enter reference datetime
                pass
        else:
            time_unit = Units(f"{CONFIG['time_units']['base_unit']} since {CONFIG['time_units']['ref_datetime']}").formatted()
        calendar = Units(calendar=(CONFIG['time_units']['calendar'] or 'standard'))
        if calendar.isvalid and calendar.isreftime:
            calendar = calendar.calendar
        else:
            calendar = 'standard'
        nc_dir.time_units = Units(units=time_unit,
                                  calendar=calendar)

    print(
        f"""
        Finished classifying files.
        0+1d: {nc_dir.source['0+1d']};
        2d: {nc_dir.source['2d']};
        3d: {nc_dir.source['3d']}.
        Time units: {nc_dir.time_units}.
        """)
    
    # Next, process files according to their dimensionality
    for dimensionality, files in nc_dir.source.items():
        for filename in files:
            with xr.open_dataset(
                os.path.join(nc_dir.source_dir, filename), decode_times=False
                ) as ds:

                # Decide first whether need to merge or split time series.
                # Doing this first avoids unnecessary duplication of 
                # modifications.

                # Split/merge time series
                if DIM_GROUPS.index(dimensionality) >= 3:
                    # Split into individual time points

                    # First, update & add any necessary data.

                    time_vars = [v for v in ds.dims.keys() if 'time' in v.lower()]
                    time_var = time_vars[0] if len(time_vars) == 1 else 'time'  # The latter default will come from config or vocab.

                    # write NC file for each time-point:
                    for time_point in ds[time_var].data:
                        if int(time_point) == time_point:
                            time_point = str(int(time_point))
                        else:
                            time_point = str(time_point).replace('.', ',')
                        output_file = nc_dir.stem[dimensionality] + time_point + '.nc'

                else:
                    # Merge with other files of same dimensionality
                    output_file = os.path.join(nc_dir.target_dir, nc_dir.stem[dimensionality] + '.nc')
                    if os.path.exists(output_file):
                        # Open dataset from output file if it exists already.
                        target_ds = xr.open_dataset(
                            os.path.join(nc_dir.source_dir, filename), 
                            decode_times=False
                        )

                        # Identify time variable, along which to merge:
                        time_vars = [v for v in ds.dims.keys() if 'time' in v.lower()]
                        time_var = time_vars[0] if len(time_vars) == 1 else 'time'  # The latter default will come from config or vocab.

                        # Merge the two datasets
                        target_ds = xr.concat([ds, target_ds], dim=time_var, combine_attrs="drop_conflicts")
                        # target_ds = merge_time_series(to_merge=ds, merge_into=target_ds)
                        print(target_ds.dims)
                    else:
                        target_ds = ds.copy(deep=False)  # This still points to original ds data in memory, as nothing will change before saving.
                    # Save target dataset, then close.
                    target_ds.to_netcdf(output_file)
                    target_ds.close()


                # Get required info from options_database
                odb_req = CONFIG['options_database']
                options_db = {
                    k.decode('utf-8'): utils.type_from_str(v.decode('utf-8'))
                            for [k, v] in ds[odb_req['variable']].data
                                if (k.decode('utf-8') in odb_req['required'] or
                                    k.decode('utf-8') in odb_req['options_to_attrs'])
                }

                # Store any options database key:value pairs listed in
                # config.yml - options_database: options_to_attrs
                # as global attributes.
                for k in odb_req['options_to_attrs']:
                    ds.attrs[k] = options_db[k]

                # Find any dimensions not present as coordinates.
                # Note: this won't (yet) add new dimensions required by changed
                # variable definitions.
                for dim, size in ds.dims.items():
                    if dim not in ds.coords and dim not in odb_req['dimensions']:
                        # Make coordinate based on missing dimension.
                        if dim not in CONFIG['add_dimensions']:
                            raise KeyError(f'Missing coordinate variable, {dim}, not found in config.yml: add_dimensions.')

                        # First, set up array of coordinate points.
                        if isinstance(CONFIG['add_dimensions'][dim]['spacing'], str):
                            # Use spacing from options database
                            try:
                                delta = options_db[CONFIG['add_dimensions'][dim]['spacing']]
                            except KeyError:
                                raise KeyError(f'In config.yml, add_dimensions: {dim}: spacing must be either a numeric value or a key from the options database. {CONFIG["add_dimensions"][dim]["spacing"]} not found in options database.')
                        else:
                            # Use spacing provided
                            try:
                                delta = float(CONFIG['add_dimensions'][dim]['spacing'])
                            except ValueError:
                                raise TypeError('In config.yml, add_dimensions: <dim>: spacing must be either a numeric value or a key from the options database.')
                        
                        first = delta/2 if 'mid' in CONFIG['add_dimensions'][dim]['position'] else 0

                        attributes = {}
                        for k, v in CONFIG['add_dimensions'][dim]['attributes'].items():
                            if k.lower() == 'units':
                                attributes['units'] = cfunits.Units(v).formatted()
                            else:
                                attributes[k.lower().replace(' ', '_')] = v
                        

                        ds = coord_from_dim(ds=ds, dim=dim, data=coord_points(n=size, delta=delta, first=first), attrs=attributes)

                # Update variables for CF compliance and any corrections
                # specified in config.yml.



def parse_options_db(ds: xr.Dataset) -> dict:
    # Get required info from options_database
    required = set([var_list for var_list in (
                    CONFIG['options_database']['required'],
                    CONFIG['options_database']['reference_datetime'],
                    CONFIG['options_database']['options_to_attrs']
    )])
    return {
        k.decode('utf-8'): utils.type_from_str(v.decode('utf-8'))
                for [k, v] in ds[CONFIG['options_database']['variable']].data
                    if k.decode('utf-8') in required
    }


def nc_basename(*args) -> str:
    strings = list(args)
    if None in strings:
        strings.pop(strings.index(None))
    if len(strings) == 0:
        raise TypeError('At least one string must be supplied.')
    s = strings[0]
    if len(strings) > 1:
        for f in strings[1:]:
            match = SequenceMatcher(a=s, b=f).find_longest_match()
            s = s[match.a: match.a + match.size]
    return s


# def process_monc(ds: xr.Dataset) -> xr.Dataset:
def process_monc(parser: MONC_Parser):
    '''
    '''
    


    # Find number of dimensions, which determines how file is merged/split
    n_dims = get_n_dims(ds)
    # print(n_dims, 'dimension(s) found')
    
    if n_dims < 3:
        # Merge with other "0d" and "1d" files (time alone or time+altitude)
        pass
    elif n_dims < 4:
        # Merge with other "2d" files (time + x + y)

        merge_time_series(new_ds=ds, merged_ds=merged_2d)
    else:
        # Split "3d" (time + 3 spatial dimensions) into individual time points.
        single_timepoints = split_time_series(ds)
        for tp in single_timepoints:
            # write new NetCDF file to target directory
            pass
    

def merge_time_series(to_merge: xr.Dataset, merge_into: xr.Dataset = None) -> xr.Dataset:
    # Check both parameters are valid. If no merge_into supplied, clone given
    # dataset into new one.
    if merge_into is None:
        return to_merge.copy(deep=False)  # Shallow copy simply serves as a new pointer to original data in memory
    


def split_time_series(ds: xr.Dataset) -> list:
    pass


def process_input(ds: xr.Dataset, var: str) -> dict:
    '''
    Returns dictionary containing useful input data, if successfully parsed.
    Returns empty dictionary if content does not match expected input data.
    '''
    if not isinstance(var, str):
        raise TypeError('Variable to look for in input file must be supplied as a string.')
    
    if var not in ds.variables:
        raise KeyError(f'{var} not in dataset variables.')

    return ds[var].attrs

    # (units, calendar) = (None, None)

    # # Refer to config file, for whether or not to use time units & calendar,
    # # and which variable to check
    # if CONFIG['time_units']['use_input_file']:
    #     try:
    #         units = ds.variables[CONFIG['time_units']['input_file_variable']].attrs['units']
    #         calendar = ds.variables[CONFIG['time_units']['input_file_variable']].attrs['calendar']
    #     except:
    #         pass
            
    # return {'time units': units, 'time calendar': calendar}


def is_monc(ds: xr.Dataset) -> bool:
    '''
    A MONC (diagnostic) output file will have 'MONC time' in its global attributes.
    '''
    return 'MONC time' in ds.attrs


def is_monc_file(filepath: str) -> bool:
    if not os.path.exists(filepath):
        raise OSError(f'Filepath not found: {filepath}')
    try:
        with xr.open_dataset(filepath, decode_times=False) as ds:
            return is_monc(ds)
    except Exception as e:
        raise e


def get_n_dims(ds: xr.Dataset) -> int:
    return max([len(ds[v].dims) for v in ds.variables if v != CONFIG['options_database']['variable']])


def test():
    # print('Testing NC directory parser:')
    # fail_dir = '/spam/eggs/beans/andspam'
    # pass_dir = os.path.join(os.path.dirname(app_dir), 'test_data')
    # for d in [fail_dir, pass_dir]:
    #     print(f'trying to fetch nc files from {d}...')
    #     try:
    #         # print(get_monc_files(d))#
    #         nc_iterator(d)
    #     except OSError as e:
    #         print(e)
    cf_dir = MONC_Parser(os.path.join(os.path.dirname(app_dir), 'test_data'), dim_bins=DIM_GROUPS)
    process_monc(cf_dir)
    

    # print('Testing nc_basename:')
    # a = 'd20200128_diagnostic_0d_172800.nc'
    # b = 'd20200128_diagnostic_1d_43200.nc'
    # print(f'Common part of {a} + {b}:\n  {nc_basename(a, b)}')
    # a = 'd20200128_diagnostic_0d_172800.nc'
    # b = None
    # print(f'Common part of {a} + {b}:\n  {nc_basename(a, b)}')
    # a = None
    # b = 'd20200128_diagnostic_0d_172800.nc'
    # print(f'Common part of {a} + {b}:\n  {nc_basename(a, b)}')


if __name__== '__main__': test()