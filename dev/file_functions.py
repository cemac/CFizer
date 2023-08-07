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


class NC_Dir:
    def __init__(self, directory: str):
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
        self.target_dir =  os.path.join(
            os.path.dirname(directory),
            f'{os.path.basename(directory)}+processed'
        )
        if not os.path.exists(self.target_dir):
            os.makedirs(self.target_dir)
        
        # Set up dictionary to contain filepaths, categorized by number of dimensions in data.
        self.source = {v: [] for v in set(DIM_GROUPS)}
        #     '0+1d': [],
        #     '2d': [],
        #     '3d': []
        # }
        
        # Set up dictionary to contain filepaths of output files.
        # 3d doesn't need to be tracked, as each output file is produced from only one input file.
        # However, for flexibility and extensibility, keep a list of all
        # target files for each dimensionality.
        self.stem = {v: None for v in set(DIM_GROUPS)}
        # {
        #     '0+1d': None,
        #     '2d': None
        # }
        
        self._time_units = None
        self.reference_datetime = []
        
    @property
    def time_units(self):
        return None if self._time_units is None else self._time_units.formatted()
    
    @time_units.setter
    def time_units(self, units: Units):
        if not isinstance(units, Units):
            raise TypeError('units must be of type cfunits.Units')
        # try:
        #     utils.decode_time_units(units)
        # except ValueError as e:
        #     raise e
        if not units.isvalid:
            raise ValueError('Units are not valid according to CF convention.')
        
        self._time_units = units


def nc_iterator(directory: str) -> None:
    '''
    Raises OSError if directory not found.
    '''
    
    try:
        nc_dir = NC_Dir(directory)
    except OSError as e:
        raise e

    # # Check directory is valid
    # if not os.path.exists(directory):
    #     raise OSError(f'Directory {directory} not found')
    
    # nc_index = {
    #     '0+1d': [],
    #     '2d': [],
    #     '3d': []
    # }

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



                


def nc_basename(*args) -> str:
    strings = list(args)
    if None in strings:
        strings.pop(strings.index(None))
    if len(strings) == 0:
        raise TypeError('At least one string must be supplied.')
    if len(strings) == 1:
        return strings[0]
    # filenames = [os.path.basename(a) for a in args]
    s = strings[0]
    for f in strings[1:]:
        match = SequenceMatcher(a=s, b=f).find_longest_match()
        s = s[match.a: match.a + match.size]
    return s


def process_monc(ds: xr.Dataset) -> xr.Dataset:
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


def process_input(ds: xr.Dataset) -> dict:
    '''
    Returns dictionary containing useful input data, if successfully parsed.
    Returns empty dictionary if content does not match expected input data.
    '''
    (units, calendar) = (None, None)

    # Refer to config file, for whether or not to use time units & calendar,
    # and which variable to check
    if CONFIG['time_units']['use_input_file']:
        try:
            units = ds.variables[CONFIG['time_units']['input_file_variable']].attrs['units']
            calendar = ds.variables[CONFIG['time_units']['input_file_variable']].attrs['calendar']
        except:
            pass
            
    return {'time units': units, 'time calendar': calendar}


def is_monc(ds: xr.Dataset) -> bool:
    '''
    A MONC (diagnostic) output file will have 'MONC time' in its global attributes.
    '''
    return 'MONC time' in ds.attrs


def get_n_dims(ds: xr.Dataset) -> int:
    return max([len(ds[v].dims) for v in ds.variables if v != CONFIG['options_database']['variable']])


def get_monc_files(directory: str) -> dict:
    '''
    Raises OSError if directory not found.
    '''
    # Check directory is valid
    if not os.path.exists(directory):
        raise OSError(f'directory {directory} not found')

    # Get list of all nc files in directory and any subdirectories
    nc_files = glob(pathname='**/*.nc', root_dir=directory, recursive=True)
    # print(len(nc_files), 'files found')
    # print(nc_files)
    monc_files = []
    other_files = []
    
    # discard any file not matching pattern for MONC output
    for f in nc_files:
        filepath = os.path.join(directory, f)
        print(filepath)
        with xr.open_dataset(filepath) as ds:
            monc_files.append(filepath) if 'MONC time' in ds.attrs else other_files.append(filepath)
                
    return {'monc': monc_files, 'other': other_files}


def test():
    print('Testing NC directory parser:')
    fail_dir = '/spam/eggs/beans/andspam'
    pass_dir = os.path.join(os.path.dirname(app_dir), 'test_data')
    for d in [fail_dir, pass_dir]:
        print(f'trying to fetch nc files from {d}...')
        try:
            # print(get_monc_files(d))#
            nc_iterator(d)
        except OSError as e:
            print(e)

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