import xarray as xr
from units import TimeUnits, format_units
from setup import *
from utils import type_from_str, generate_coords #, performance_time
from cfunits import Units
from dask.delayed import Delayed
from time import perf_counter
# from cfizer import OPTIONS_DATABASE, CONFIG


def get_n_dims(dataset: xr.Dataset) -> int:
    return max([
        len(dataset[v].dims) for v in dataset.variables 
        if v != OPTIONS_DATABASE['variable']
        ])


def is_monc(dataset: xr.Dataset) -> bool:
    return MONC_ID_ATTR in set(dataset.attrs).union(set(dataset.variables))


# @performance_time
def split_ds(dataset: xr.Dataset, 
             shared: dict, 
             var: str = 'time') -> tuple[list[xr.Dataset], str]:
    start_time = perf_counter()
    log = []
    if var not in dataset.dims:
        raise AttributeError(f'split_ds: {var} not in dataset dimensions.')
    if shared['verbose']: 
        log.append(
            f'Process {os.getpid()}: splitting dataset '
            f'{dataset.attrs["title"]} by {var}.'
        )
        # print(
        #     f'Process {os.getpid()}: splitting dataset '
        #     f'{dataset.attrs["title"]} by {var}.'
        # )
    base_title = dataset.attrs['title'].strip('_ +,.&') + '_'
    grouped = {point: ds for (point, ds) in dataset.groupby(var)}
    for point, ds in grouped.items():
        # Append relevant coordinate value to title
        grouped[point].attrs['title'] = base_title + str(int(point))
        if shared['verbose']: 
            log.append(
                f"Process {os.getpid()}: Created new dataset with title, "
                f"{grouped[point].attrs['title']}"
            )
            # print(
            #     f"Process {os.getpid()}: Created new dataset with title, "
            #     f"{grouped[point].attrs['title']}"
            # )
    split = list(grouped.values())
    log.append(
        f"Process {os.getpid()}: split_ds took {perf_counter() - start_time} "
        f"seconds."
    )
    return (split, log)


# @performance_time
def ds_to_nc(ds: xr.Dataset, 
             filepath: str, 
             encodings: dict = None, 
             compress: bool = False,
             shared: dict = None) -> None:
    
    if compress:
        ds.to_netcdf(
            path=filepath, 
            encoding=encodings,
            engine='netcdf4',
            format='netCDF4',
        )
    else:
        ds.to_netcdf(
            path=filepath, 
            encoding=encodings
        )
    

# @performance_time
def ds_to_nc_dask(ds: xr.Dataset, 
                  filepath: str, 
                  encodings: dict = None, 
                  compress: bool = False,
                  shared: dict = None) -> Delayed:
    if compress:
        return ds.to_netcdf(
            path=filepath, 
            encoding=encodings, 
            engine='netcdf4',
            format='netCDF4',
            compute=False
            )
    else:
        return ds.to_netcdf(
            path=filepath, 
            encoding=encodings, 
            compute=False
        )
    # writer.compute()


# @performance_time
def perform_write(writer: Delayed,
                  shared: dict = None) -> None:
    writer.compute()


class MoncDs:
    def __init__(self, 
                 dataset: xr.Dataset,
                 time_variable:str=None,
                 n_dims:int=None,
                 title:str=None
                 ) -> None:  # , group:DsGroup=None,
        
        self.log = []
        self.warnings = []
        self.ds = dataset
        # self.time_units = shared['time_units']
        # self.vocab = shared['vocabulary']

        if time_variable:
            self.time_var = time_variable
        # elif group:
        #     self.time_var = group.time_var
        else:
            raise AttributeError(
                "DsGroup: time variable must be specified via the "
                "time_variable parameter."
            )
        
        if n_dims:
             self.n_dims = n_dims
        # elif group:
        #     self.n_dims = group.n_dims
        else:
            # Find number of dimensions from dataset
            self.n_dims = get_n_dims(self.ds) - 1  # Subtract 1 for time.
        
        if title:
            self.ds.attrs['title'] = title
        # elif group: 
        #     self.ds.attrs['title'] = group.stem.strip('_ ') + group.name if group.name not in group.stem else group.stem.strip('_ ')
        elif 'title' not in self.ds.attrs:
            raise AttributeError(
                "MoncDs: dataset title is missing."
            )
        
    def cfize(self, shared: dict):
        '''
        vocabulary can't be global, because this creates conflicts in parallel
        processing.
        '''

        if shared['verbose']: 
            start_time = perf_counter()
            self.log.append(
                f"Process {os.getpid()}: cfize running on MONC dataset "
                f"{self.ds.attrs['title']} with {self.n_dims} spatial "
                f"dimensions."
            )
            # print(
            #     f"Process {os.getpid()}: cfize running on MONC dataset "
            #     f"{self.ds.attrs['title']} with {self.n_dims} spatial "
            #     f"dimensions."
            # )

        # Pull required options_database info into dictionary
        self.get_options(
            fields=CONFIG['options_to_attrs']
        )
        
        # Convert any required options_database pairs to attributes (with 
        # appropriate np.dtypes)
        self.ds.attrs.update(self.options_db)

        # Add global attributes recommended by CF (some values set at 
        # initialisation)
        self.add_cf_attrs()
        
        # Add any missing coordinates
        self.missing_coords()
        
        # update all (non-options-database) variables in either shared['vocabulary'][dim] or 
        # dataset:
        var_list = set(self.ds.variables)
        var_list.discard(OPTIONS_DATABASE['variable'])
        [var_list.discard(v) for v in CONFIG['global_to_variable'].keys()]
        self.cfize_variables(variables=var_list, shared=shared)

        if shared['verbose']: 
            self.log.append(
                f"Process {os.getpid()}: MoncDs.cfize took "
                f"{perf_counter() - start_time} seconds."
            )
            # print(
            #     f"Process {os.getpid()}: MoncDs.cfize took "
            #     f"{perf_counter() - start_time} seconds."
            # )

        # Return processed dataset
        return self.ds

    def missing_coords(self, dim:str=None):
        '''
        Looks for any dimensions that currently don't also exist as coordinate 
        variables.
        dim: str (optional)
            dimension from which coordinate should be created.
        '''
        # Create list of each dim not in coords;
        # Find dtype of existing (spatial) coords (need to deal with conflicts?)
        dim_type = np.dtype('float64')  # Default value, in case none is found
        missing = []
        # If dim specified,
        if dim and dim in self.ds.dims:
            # convert to a list of length one
            missing = [dim]
            if dim in self.ds.coords:
                # If coordinate variable already exists, match its type and
                # add any new attributes.
                dim_type = self.ds.coords[dim].dtype
        else:
            for d in self.ds.dims:
                # Set type to match existing spatial coordinates (assume for now all 
                # are the same)
                if d in self.ds.coords:
                    if d != self.time_var:
                        dim_type = self.ds.coords[d].dtype
                elif d not in OPTIONS_DATABASE['dimensions']:
                    missing.append(d)
            
        
        # For each coord in list:
        for d in missing:
            # Look for any required parameters & attributes in
                # 1. CONFIG
                # 2. VOCAB
            attributes = {}

            # Look in CONFIG for parameters required to set up variable:
            # grid spacing & position (edge/centre).
            if d in CONFIG['new_coordinate_variables'].keys():
                config = CONFIG['new_coordinate_variables'][d]
                try:
                    spacing = config['spacing'] if isinstance(config['spacing'], (int, float)) else self.options_db[config['spacing']]
                    midpoint = 'cent' in config['position'] or 'mid' in config['position']
                    for k, v in config['attributes'].items():
                        if k.lower() == 'units' and Units(v).isvalid:
                            attributes['units'] = format_units(v) 
                        else:
                            attributes[k.lower().replace(' ', '_')] = v if k != 'standard_name' else v.replace(' ', '_')
                except KeyError:
                    raise ConfigError(
                        f"Processing {self.n_dims}d files: Parameters needed "
                        f"to set up {d} as coordinate variable are missing "
                        f"from config file, or required key not found in "
                        f"options_database. "
                        f"For each new coordinate variable, spacing, position "
                        f"and attributes must be specified."
                    )
            else:
                raise ConfigError(
                        f"Processing {self.n_dims}d files: Parameters needed "
                        f"to set up {d} as coordinate variable are missing "
                        f"from config file. "
                        f"For each new coordinate variable, spacing, position "
                        f"and attributes must be specified. Spacing must be "
                        f"either a number or a parameter in the "
                        f"options_database."
                )
            
            # Skip looking for attributes in vocabulary, as these will be 
            # updated in cfize_variables anyway.
            
            # Generate coordinate points as np.ndarray of required dtype
            points = generate_coords(number=self.ds.dims[d],
                                    spacing=spacing,
                                    midpoint=midpoint,
                                    data_type=dim_type)
            
            # Assign dim, attrs & coord points to new xarray.Variable
            new_var = xr.Variable(dims=d, attrs=attributes, data=points)

            # Assign new variable to dataset -> overwrite existing.
            self.ds = self.ds.assign(variables={d: new_var})

    def update_units(self,
                     var: str,
                     shared: dict,
                     updates:dict=None):
        '''Add/update units for specified variable.'''
        
        time_units = shared['time_units']
        if var in CONFIG['new_coordinate_variables']:
            if not updates:
                updates = CONFIG['new_coordinate_variables'][var]['attributes']
            else:
                # Get any info missing from updates out of 
                # CONFIG['new_coordinate_variables'][var]
                [
                    updates.update({k: v}) 
                    for k, v in CONFIG[
                        'new_coordinate_variables'][var]['attributes'].items() 
                    if k not in updates
                ]
                
        if not updates or 'units' not in updates:
            # If no new unit, check existing is present & CF compliant:
            if 'units' in self.ds[var].attrs:
                # If standard_name, look up
                if 'standard_name' in self.ds[var].attrs:
                    # TODO: look up standard names dictionary
                    pass
                # Otherwise, check with cfunits.Units().isvalid
                else:
                    if Units(self.ds[var].attrs['units']).isvalid:
                        if (
                            (var == self.time_var or 
                            ('axis' in self.ds[var].attrs and 
                             self.ds[var].attrs['axis'] == 'T')
                            ) and not shared['time_units']
                        ):
                            if not Units(self.ds[var].attrs['units']).isreftime:
                                raise VocabError(
                                    f"update_units: "
                                    f"{self.ds[var].attrs['units']} is an "
                                    f"invalid unit for {var}: the CF format "
                                    f"for time coordinate variables is "
                                    f"'<time unit> since <reference date[time]"
                                    f"'. No updated unit found in vocabulary."
                                )
                        self.ds[var].attrs['units'] = format_units(
                            self.ds[var].attrs['units']
                        )
                    elif self.ds[var].attrs['units'].lower() == 'fraction':
                        self.ds[var].attrs['units'] = '1'
                    else:
                        # pass
                        # print(
                        #     f"{self.ds[var].attrs['units']}: invalid units for "
                        #     f"variable {var} in dataset "
                        #     f"with title {self.ds.attrs['title']}, and no "
                        #     f"new units specified in vocabulary under "
                        #     f"dimension {self.n_dims}."
                        # )
                        raise VocabError(
                            f"update_units: {self.ds[var].attrs['units']}: "
                            f"invalid units for variable {var} in dataset "
                            f"with title {self.ds.attrs['title']}, and no "
                            f"new units specified in vocabulary under "
                            f"dimension {self.n_dims}."
                        )
            else:
                # print(
                #     f"No units specified for {var} in vocabulary "
                #     f"({self.n_dims}d), and no existing units found."
                # )
                # return
                raise VocabError(
                    f"update_units: CF conventions require that units be "
                    f"specified for each variable. No units specified in "
                    f"vocabulary ({self.n_dims}d) for {var}, and no existing "
                    f"units found."
                )
        else:
            # Units specified in vocabulary.
            # if variable is a time coordinate:
            if var == (
                self.time_var or 
                ('axis' in self.ds[var].attrs and 
                 self.ds[var].attrs['axis'] == 'T') or
                ('axis' in updates and updates['axis'] == 'T')
            ):
                # If supplied, process new units via TimeUnits class.
                # Build any missing attributes from time_units parameter.
                if 'calendar' in updates :
                    calendar = updates['calendar'] 
                elif 'calendar' in self.ds[var].attrs and Units(
                    calendar=self.ds[var].attrs['calendar']).isvalid:
                    calendar = self.ds[var].attrs['calendar']
                else:
                    calendar = time_units.calendar
                self.ds[var].attrs['calendar'] = calendar
                if TimeUnits(units=updates['units']).isreftime:
                    u = updates['units']
                else:
                    if time_units:
                        u = f"{updates['units']} since {time_units.reftime.isoformat(sep=' ')}"
                    elif 'units' in self.ds[var].attrs and Units(self.ds[var].attrs['units']).isreftime:
                        u = f"{updates['units']} since {Units(self.ds[var].attrs['units']).reftime.isoformat(sep=' ')}"
                    else:
                        raise VocabError(
                                    f"update_units: "
                                    f"{updates['units']} is an "
                                    f"invalid unit for {var}: the CF format "
                                    f"for time coordinate variables is "
                                    f"'<time unit> since <reference date[time]"
                                    f"'. No reference date[time] found or "
                                    f"specifed."
                        )
                new_units = TimeUnits(units=u, 
                                    calendar=calendar)
                # Assign T to axis.
                self.ds[var].attrs['axis'] = 'T'

            else:
                # Check any new unit specified is valid
                if Units(updates['units']).isvalid:
                    new_units = Units(updates['units'])
                else:
                    raise VocabError(
                        f"update_units: Invalid units specified in vocabulary "
                        f"for ({self.n_dims}d) variable {var}."
                    )
                
                # If spatial coordinate variable, add/update any specified or 
                # implied axis.
                if var in self.ds.dims:
                    if 'axis' in updates:
                        self.ds[var].attrs['axis'] = updates['axis']
                    elif var in CONFIG['new_coordinate_variables'] and 'axis' in CONFIG['new_coordinate_variables'][var]['attributes']:
                        self.ds[var].attrs['axis'] = CONFIG[
                            'new_coordinate_variables'][var][
                                'attributes']['axis']
                    elif var[0].lower() in {'x', 'y', 'z'}:
                        if not shared['quiet']:
                            self.log.append(
                                f"MoncDs.update_units: assumed axis "
                                f"{var[0].upper()} for variable {var}."
                            )
                            # print(
                            #     f"MoncDs.update_units: assumed axis "
                            #     f"{var[0].upper()} for variable {var}."
                            # )
                        self.ds[var].attrs['axis'] = var[0].upper()
                    elif 'axis' not in self.ds[var].attrs:
                        # Prompt user if required.
                        raise ConfigError(
                            f"update_units: Axis should be specified for "
                            f"spatial coordinate "
                            f"variables. No axis attribute found in dataset, "
                            f"vocabulary or configuration file for variable "
                            f"{var}."
                        )

                    # Add positive attribute if required/specified (not expected).

            # If any new unit equivalent but not equal to any existing, run 
            # conversion.
            if 'units' in self.ds[var].attrs:
                old_units = Units(self.ds[var].attrs['units'])
                if new_units.equivalent(old_units) and not new_units.equals(old_units):
                    Units.conform(
                        x=self.ds[var].data, 
                        from_units=old_units, 
                        to_units=new_units, 
                        inplace=True
                    )

            # Apply new units if present, or update format of existing.
            # Don't use cfunits.Units.formatted() method, because its ordering
            # is weird; instead, call a function that recognises '/' as 
            # inverting any exponent (including implied 1), interprets '^' as
            # exponent operator
            self.ds[var].attrs['units'] = format_units(new_units)

    def cfize_variables(self, 
                        variables: list|set|tuple,
                        shared: dict):
        
        if shared['verbose']: 
            self.log.append(
                f"Process {os.getpid()}: cfizing variables on MONC dataset "
                f"{self.ds.attrs['title']} with {self.n_dims} spatial "
                f"dimensions."
            )
            # print(
            #     f"Process {os.getpid()}: cfizing variables on MONC dataset "
            #     f"{self.ds.attrs['title']} with {self.n_dims} spatial "
            #     f"dimensions."
            # )
        
        # Process any coordinate variables first, as changes in these may affect
        # processing of data variables.
        coord_vars = [v for v in variables if v in self.ds.dims]
        variables = coord_vars + [v for v in variables if v not in coord_vars]

        # For each variable in args:
        for var in variables:
            if var not in self.ds.variables:
                raise AttributeError(
                    f"cfize_variables: Variable {var} not found in dataset."
                )
                # ENHANCEMENT: If variable not in dataset, as is, look for 
                # wildcards or use it as a stem. E.g. time_series* should identify 
                # time_series_300_1800 as a match.

            # Find variable's updates in shared['vocabulary'][dim][variable].
            try:
                updates = shared['vocabulary'][self.n_dims][var]
            except KeyError:
                # print(
                #     f"cfize_variables: {var} not found in vocabulary "
                #     f"(dimension {self.n_dims})."
                # )
                raise VocabError(
                    f"cfize_variables: {var} not found in vocabulary "
                    f"(dimension {self.n_dims})."
                )
                # # check/update any units already present.
                # self.update_units(var=var, time_units=shared['time_units'])
                # continue

            # update variable's dimensions if required
            if 'dimension_changes' in updates:
                # Check specified dim(s) exist(s)
                if not all([
                    k in self.ds[var].dims
                    for k in updates['dimension_changes'].keys()
                ]):
                    raise VocabError(
                        f"cfize_variables: {k}, specified in vocabulary - "
                        f"dimension_changes, is not a dimension of {var}."
                    )
                # swap_dims
                self.ds[var] = self.ds[var].swap_dims(updates['dimension_changes'])

                # If new dim not in coords, call missing_coords function, 
                # perhaps with new coord name as argument
                for d in set(self.ds[var].dims) - set(self.ds[var].dims).intersection(set(self.ds.coords)):
                    self.missing_coords(dim=d)
                
            # Add/update names:
            # standard_name &/or long_name
            if 'standard_name' in updates:
                # TODO: check standard_name is in list of standard names.
                # This should be done in setting up vocabulary.yml, but best to
                # check here nonetheless.

                self.ds[var].attrs['standard_name'] = updates['standard_name']

            if 'long_name' in updates:
                self.ds[var].attrs['long_name'] = updates['long_name']

            # check at least one present
            if 'standard_name' not in self.ds[var].attrs and 'long_name' not in self.ds[var].attrs:
                # self.warnings.append(VocabError(
                #     f"{var}: CF requires at least one of standard_name and "
                #     f"long_name to be assigned to each variable."
                # ))
                raise VocabError(
                    f"cfize_variables: {var}: CF requires at least one of "
                    f"standard_name and long_name to be assigned to each "
                    f"variable."
                )
            
            # Add/update units
            self.update_units(var=var, 
                              updates=updates, 
                              shared=shared)
            
            # Update variable name if specified.
            if 'updated_name' in updates:
                self.ds = self.ds.rename({var: updates['updated_name']})
                '''
                TODO: If updated variable is the time variable, need to pass 
                back updated name, or probe for this change after return.
                '''
                if var == self.time_var:
                    self.time_var = updates['updated_name']
            
            # Convert any perturbation variables to absolute values.
            if 'perturbation_to_absolute' in updates and updates['perturbation_to_absolute']:
                if 'reference_variable' not in updates:
                    raise VocabError(
                        f"cfize_variables: {variable}: If "
                        f"perturbation_to_absolute is True, "
                        f"reference_variable must contain the name of the "
                        f"variable containing reference value(s).")
                
                # Get reference variable data array from (updated) vocabulary.
                ref_array = updates['ref_array']
                
                # Add reference DataArray to perturbation variable DataArray, & 
                # re-assign to variable as new array with existing attributes, 
                # dimensions, etc.
                absolute = self.ds[var] + ref_array
                # Adding two xarray.DataArrays doesn't preserve attributes,
                # so need to create new data array combining new data with
                # existing attributes.
                # print("Trying to build new data array from sum of perturbation & reference.")
                # print("absolute:", absolute)
                # print(absolute.data)
                # print("existing:", self.ds[var])
                self.ds[var] = xr.DataArray(
                    data=absolute.data,
                    coords=self.ds[var].coords,
                    dims=self.ds[var].dims,
                    name=self.ds[var].name,
                    attrs=self.ds[var].attrs
                )

    def add_cf_attrs(self, **kwargs):
            # TODO: <version> and <url> and any other required data will be assigned during packaging.

            defaults = {}  # {attr: None for attr in CF_ATTRIBUTES}
            # defaults['title'] = kwargs['title'] if 'title' in kwargs
            defaults['history'] = f'{datetime.now().isoformat(timespec="minutes", sep=" ")}: output files processed using CFizer version <version>, <url>.'
            defaults['conventions'] = f"CF-{CF_VERSION}"
            for attr in CF_ATTRIBUTES:
                # Only overwrite existing value if a new value was specified in
                # config file.
                if attr in self.ds.attrs :
                    if attr in CONFIG and CONFIG[attr] is not None:
                        self.ds.attrs[attr] = CONFIG[attr]
                elif (attr in CONFIG and CONFIG[attr]) or (attr in defaults and defaults[attr]):
                    self.ds.attrs[attr] = CONFIG[attr] if attr in CONFIG else defaults[attr]
                else:
                    # print(
                    #     f"{attr} must be included as a global attribute, "
                    #     f"specified in config.yml."
                    # )
                    self.warnings.append(ConfigWarning(
                        f"add_cf_attrs: CF-1.10, section 2.6.2, recommends "
                        f"{attr} be included as a global attribute. This can "
                        f"be specified in config.yml."
                    ))

    def get_options(self, fields: list|set|tuple = None) -> dict:
        '''
        Note: the xarray.Dataset.assign_attrs() method isn't suitable here, as it
        creates a new dataset, rather than updating the existing one.
        '''
        # Get required fields from argument; import all if none supplied.
        if fields is None:
            if self.ds[OPTIONS_DATABASE['variable']].dtype == 'S1':
                # This is used if concat_characters=False options is used in
                # xarray.open_dataset. It seems to create some problems down
                # the line, so is not advised.
                options = {
                    ''.join([c.decode('utf-8') for c in k]): 
                    type_from_str(''.join([c.decode('utf-8') for c in v]))
                    for k, v in self.ds[OPTIONS_DATABASE['variable']].data
                }  # If dataset opened with concat_characters=False
            else:
                options = {
                    k.decode('utf-8'): type_from_str(v.decode('utf-8'))
                    for [k, v] in self.ds[OPTIONS_DATABASE['variable']].data}
        else:
            if self.ds[OPTIONS_DATABASE['variable']].dtype == 'S1':
                options = {
                    ''.join([c.decode('utf-8') for c in k]): 
                    type_from_str(''.join([c.decode('utf-8') for c in v]))
                    for k, v in self.ds[OPTIONS_DATABASE['variable']].data
                    if ''.join([c.decode('utf-8') for c in k]) in fields
                }  # If dataset opened with concat_characters=False
            else:
                options = {
                    k.decode('utf-8'): type_from_str(v.decode('utf-8'))
                    for [k, v] in self.ds[OPTIONS_DATABASE['variable']].data
                    if k.decode('utf-8') in fields}
        
        # Drop any options that are reset by DEPHY, if used.
        if (all([opt in options for opt in DEPHY_OPTIONS]) and
            all([options[opt] for opt in DEPHY_OPTIONS])):
            for a in DROP_FOR_DEPHY:
                if a in options:
                    options.pop(a)

        self.options_db = options
    