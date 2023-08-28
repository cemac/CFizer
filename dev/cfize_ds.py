import xarray as xr
from units import TimeUnits
from time import sleep
import os.path as op
from setup import *
from utils import type_from_str, format_units, generate_coords
from cfunits import Units
from groups import DsGroup


def get_n_dims(dataset: xr.Dataset) -> int:
    return max([
        len(dataset[v].dims) for v in dataset.variables 
        if v != OPTIONS_DATABASE['variable']
        ])


class MoncDs:
    def __init__(self, 
                 dataset: xr.Dataset,
                 time_units: TimeUnits,
                 time_variable:str=None,
                 n_dims:int=None,
                 title:str=None,
                 group:DsGroup=None,
                 ) -> None:
        
        self.ds = dataset
        self.time_units = time_units

        if time_variable:
            self.time_var = time_variable
        elif group:
            self.time_var = group.time_var
        else:
            raise AttributeError(
                "DsGroup: time variable must be specified via either of the "
                "time_variable or group parameters."
            )
        
        if n_dims:
             self.n_dims = n_dims
        elif group:
            self.n_dims = group.n_dims
        else:
            # Find number of dimensions from dataset
            self.n_dims = get_n_dims(self.ds) - 1  # Subtract 1 for time.
        
        if title:
            self.ds.attrs['title'] = title
        elif group: 
            self.ds.attrs['title'] = group.stem.strip('_ ') + group.name if group.name not in group.stem else group.stem.strip('_ ')
        elif 'title' not in self.ds.attrs:
            raise AttributeError(
                "MoncDs: dataset title is missing. This can be defined via the "
                "title or group parameters."
            )
        
    def cfize(self):
        # Pull required options_database info into dictionary
        self.options_db = get_options(
            dataset=self.ds, 
            fields=CONFIG['options_to_attrs']
        )

        # Convert any required options_database pairs to attributes (with 
        # appropriate np.dtypes)
        self.ds.attrs.update(self.options_db)

        # Add global attributes recommended by CF (some values set at 
        # initialisation)
        add_cf_attrs(dataset=self.ds)

        # Add any missing coordinates
        self.missing_coords()

        # update all (non-options-database) variables in either VOCAB[dim] or 
        # dataset:
        var_list = set(self.ds.variables)
        var_list.discard(OPTIONS_DATABASE['variable'])
        [var_list.discard(v) for v in CONFIG['global_to_variable'].keys()]
        self.cfize_variables(variables=var_list)

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
                    raise KeyError(
                        f"Processing {self.n_dims}d files: Parameters needed "
                        f"to set up {d} as coordinate variable are missing "
                        f"from config file, or required key not found in "
                        f"options_database. "
                        f"For each new coordinate variable, spacing, position "
                        f"and attributes must be specified."
                    )
            else:
                raise KeyError(
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
                     updates:dict=None):
        '''Add/update units for specified variable.'''
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
                # for k, v in CONFIG['new_coordinate_variables'][var]['attributes'].items():
                #     if k not in updates:
                #         updates[k] = v
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
                        self.ds[var].attrs['units'] = format_units(
                            self.ds[var].attrs['units']
                        )
                    elif self.ds[var].attrs['units'].lower() == 'fraction':
                        self.ds[var].attrs['units'] = '1'
                    else:
                        print(
                            f"{self.ds[var].attrs['units']}: invalid units for "
                            f"variable {var} in dataset "
                            f"with title {self.ds.attrs['title']}, and no "
                            f"new units specified in vocabulary under "
                            f"dimension {self.n_dims}."
                        )
                        # raise ValueError(
                        #     f"{self.ds[var].attrs['units']}: invalid units for "
                        #     f"variable {var} in dataset "
                        #     f"with title {self.ds.attrs['title']}, and no "
                        #     f"new units specified in vocabulary under "
                        #     f"dimension {self.n_dims}."
                        # )
            else:
                print(
                    f"No units specified for {var} in vocabulary "
                    f"({self.n_dims}d), and no existing units found."
                )
                return
                # raise AttributeError(
                #     f"CF conventions require that units be specified for each "
                #     f"variable. No units specified in vocabulary "
                #     f"({self.n_dims}d) for {var}, and no existing units found."
                # )
        else:
            # Units specified in vocabulary.
            # if variable is a time coordinate:
            if var == self.time_var:
                # If supplied, process new units via TimeUnits class.
                # Build any missing attributes from time_units parameter.
                if 'calendar' in updates :
                    calendar = updates['calendar'] 
                elif 'calendar' in self.ds[var].attrs and Units(
                    calendar=self.ds[var].attrs['calendar']).isvalid:
                    calendar = self.ds[var].attrs['calendar']
                else:
                    calendar = self.time_units.calendar
                self.ds[var].attrs['calendar'] = calendar
                if not TimeUnits(units=updates['units']).isreftime:
                    u = f"{updates['units']} since {self.time_units.reftime.isoformat()}"
                new_units = TimeUnits(units=u, 
                                    calendar=calendar)
                # Assign T to axis.
                self.ds[var].attrs['axis'] = 'T'

            else:
                # Check any new unit specified is valid
                if Units(updates['units']).isvalid:
                    new_units = Units(updates['units'])
                else:
                    raise ValueError(
                        f"Invalid units specified for ({self.n_dims}d) "
                        f"variable {var}."
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
                    elif 'axis' not in self.ds[var].attrs:
                        # Prompt user if required.
                        raise AttributeError(
                            f"Axis should be specified for spatial coordinate "
                            f"variables. No axis attribute found in dataset or "
                            f"configuration file for "
                            f"variable {var}."
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
                        variables: list|set|tuple):
        
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

            # Find variable's updates in VOCAB[dim][variable].
            try:
                updates = vocabulary[self.n_dims][var]
            except KeyError:
                print(
                    f"cfize_variables: {var} not found in vocabulary "
                    f"(dimension {self.n_dims})."
                )
                # check/update any units already present.
                self.update_units(var=var)
                continue

            # update variable's dimensions if required
            if 'dimension_changes' in updates:
                # Check specified dim(s) exist(s)
                if not all([
                    k in self.ds[var].dims
                    for k in updates['dimension_changes'].keys()
                ]):
                    raise KeyError(
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
                raise KeyError("CF requires at least one of standard_name and "
                            "long_name to be assigned to each variable.")
            
            # Add/update units
            self.update_units(var=var, updates=updates)
            # if 'units' not in updates:
            #     # If no new unit, check existing is present & CF compliant:
            #     if 'units' in self.ds[var].attrs:
            #         # If standard_name, look up
            #         if 'standard_name' in self.ds[var].attrs:
            #             # TODO: look up standard names dictionary
            #             pass
            #         # Otherwise, check with cfunits.Units().isvalid
            #         else:
            #             if Units(self.ds[var].attrs['units']).isvalid:
            #                 self.ds[var].attrs['units'] = format_units(
            #                     self.ds[var].attrs['units']
            #                 )
            #             else:
            #                 raise ValueError(
            #                     f"Invalid units for variable {var} in dataset "
            #                     f"with title {self.ds.attrs['title']}, and not "
            #                     f"new units specified in vocabulary under "
            #                     f"dimension {self.n_dims}."
            #                 )
            # else:
            #     # Units specified in vocabulary.
            #     # if variable is a time coordinate:
            #     if var == self.time_var:
            #         # If supplied, process new units via TimeUnits class.
            #         # Build any missing attributes from time_units parameter.
            #         if 'calendar' in updates :
            #             calendar = updates['calendar'] 
            #         elif 'calendar' in self.ds[var].attrs and Units(
            #             calendar=self.ds[var].attrs['calendar']).isvalid:
            #             calendar = self.ds[var].attrs['calendar']
            #         else:
            #             calendar = self.time_units.calendar
            #         self.ds[var].attrs['calendar'] = calendar
            #         if not TimeUnits(units=updates['units']).isreftime:
            #             u = f"{updates['units']} since {self.time_units.reftime.isoformat()}"
            #         new_units = TimeUnits(units=u, 
            #                             calendar=calendar)
            #         # Assign T to axis.
            #         self.ds[var].attrs['axis'] = 'T'

            #     else:
            #         # Check any new unit specified is valid
            #         if Units(updates['units']).isvalid:
            #             new_units = Units(updates['units'])
            #         else:
            #             raise ValueError(
            #                 f"Invalid units specified in vocabulary (dim "
            #                 f"{self.n_dims}) for variable {var}."
            #             )

            #     # If any new unit equivalent but not equal to any existing, run 
            #     # conversion.
            #     if 'units' in self.ds[var].attrs:
            #         old_units = Units(self.ds[var].attrs['units'])
            #         if new_units.equivalent(old_units) and not new_units.equals(old_units):
            #             Units.conform(
            #                 x=self.ds[var].data, 
            #                 from_units=old_units, 
            #                 to_units=new_units, 
            #                 inplace=True
            #             )

            #     # Apply new units if present, or update format of existing.
            #     # Don't use cfunits.Units.formatted() method, because its ordering
            #     # is weird; instead, call a function that recognises '/' as 
            #     # inverting any exponent (including implied 1), interprets '^' as
            #     # exponent operator
            #     self.ds[var].attrs['units'] = format_units(new_units)

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
                    raise KeyError(
                        f"{variable}: If perturbation_to_absolute is True, "
                        f"reference_variable must contain the name of the "
                        f"variable containing reference value(s).")
                # Look up reference variable's dimension
                ref_var = updates['reference_variable']
                if ref_var not in self.ds.variables:
                    # Look up NC file containing reference variable in other 
                    # dimension groups

                    # # Open dataset of corresponding group
                    # while not ref_group.processed:  # op.exists(ref_ds_path):
                    #        sleep(1)
                    
                    # '''
                    # This assumes required variable is in a merged dataset. Will 
                    # need more rigorous mapping if this cannot be assumed.
                    # '''
    
                    # NB: can't wait for reference_vars[ref_var]['filepath']
                    # to appear here, as it will be updated on a different
                    # process, with its own copy of the dictionary.
                    
                    with xr.open_dataset(
                        reference_vars[ref_var]['filepath'],
                        decode_times=False
                    ) as ref_ds:
                        # Pull out reference variable as DataArray.
                        ref_array = ref_ds[ref_var]
                else:
                    # Pull out reference variable as DataArray.
                    ref_array = self.ds[ref_var]

                # Add reference DataArray to perturbation variable DataArray, & 
                # re-assign to variable as new array with existing attributes, 
                # dimensions, etc.
                absolute = self.ds[var] + ref_array
                # Adding two xarray.DataArrays doesn't preserve attributes,
                # so need to create new data array combining new data with
                # existing attributes.
                self.ds[var] = xr.DataArray(
                    data=absolute.data,
                    coords=self.ds[var].coords,
                    dims=self.ds[var].dims,
                    name=self.ds[var].name,
                    attrs=self.ds[var].attrs,
                    indexes=self.ds[var].indexes
                )

        # Return updated dataset.


# def cfize_dataset(dataset: xr.Dataset,
#                   n_dims: int,
#                   title: str,
#                   time_var: str,
#                   time_units: TimeUnits) -> xr.Dataset:
    
#     # Pull required options_database info into dictionary
#     options_db = get_options(dataset=dataset, fields=CONFIG['options_to_attrs'])

#     # Convert any required options_database pairs to attributes (with 
#     # appropriate np.dtypes)
#     dataset.attrs.update(options_db)

#     # Add global attributes recommended by CF (some values set at 
#     # initialisation)
#     add_cf_attrs(dataset=dataset)
    
#     # Add any missing coordinates
#     dataset = missing_coords(dataset=dataset,
#                              n_dims=n_dims,
#                              time_var=time_var,
#                              options=options_db)

#     # update all variables in either VOCAB[dim] or dataset:
#     var_list = set(dataset.variables)
#     var_list.discard(OPTIONS_DATABASE['variable'])
#     [var_list.discard(v) for v in CONFIG['global_to_variable'].keys()]
#     dataset = cfize_variables(dataset=dataset,
#                               n_dims=n_dims,
#                               time_var=time_var,
#                               options=options_db,
#                               time_units=time_units,
#                               variables=var_list)

#     # Return updated dataset.
#     pass


def get_options(dataset: xr.Dataset, fields: list|set|tuple = None) -> dict:
    '''
    Note: the xarray.Dataset.assign_attrs() method isn't suitable here, as it
    creates a new dataset, rather than updating the existing one.
    '''
    # Get required fields from argument; import all if none supplied.
    if fields is None:
        if dataset[OPTIONS_DATABASE['variable']].dtype == 'S1':
            # This is used if concat_characters=False options is used in
            # xarray.open_dataset. It seems to create some problems down
            # the line, so is not advised.
            options = {
                ''.join([c.decode('utf-8') for c in k]): 
                type_from_str(''.join([c.decode('utf-8') for c in v]))
                for k, v in dataset[OPTIONS_DATABASE['variable']].data
            }  # If dataset opened with concat_characters=False
        else:
            options = {
                k.decode('utf-8'): type_from_str(v.decode('utf-8'))
                for [k, v] in dataset[OPTIONS_DATABASE['variable']].data}
    else:
        if dataset[OPTIONS_DATABASE['variable']].dtype == 'S1':
            options = {
                ''.join([c.decode('utf-8') for c in k]): 
                type_from_str(''.join([c.decode('utf-8') for c in v]))
                for k, v in dataset[OPTIONS_DATABASE['variable']].data
                if ''.join([c.decode('utf-8') for c in k]) in fields
            }  # If dataset opened with concat_characters=False
        else:
            options = {
                k.decode('utf-8'): type_from_str(v.decode('utf-8'))
                for [k, v] in dataset[OPTIONS_DATABASE['variable']].data
                if k.decode('utf-8') in fields}
    
    # Drop any options that are reset by DEPHY, if used.
    if (all([opt in options for opt in DEPHY_OPTIONS]) and
        all([options[opt] for opt in DEPHY_OPTIONS])):
        for a in DROP_FOR_DEPHY:
            if a in options:
                options.pop(a)

    return options
    

def add_cf_attrs(dataset: xr.Dataset, **kwargs):
        # TODO: <version> and <url> and any other required data will be assigned during packaging.
        defaults = {}  # {attr: None for attr in CF_ATTRIBUTES}
        # defaults['title'] = kwargs['title'] if 'title' in kwargs
        defaults['history'] = f'{datetime.now().isoformat(timespec="minutes", sep=" ")}: output files processed using CFizer version <version>, <url>.'
        defaults['conventions'] = f"CF-{CF_VERSION}"
        for attr in CF_ATTRIBUTES:
            # Only overwrite existing value if a new value was specified in
            # config file.
            if attr in dataset.attrs :
                if attr in CONFIG and CONFIG[attr] is not None:
                    dataset.attrs[attr] = CONFIG[attr]
            elif (attr in CONFIG and CONFIG[attr]) or (attr in defaults and defaults[attr]):
                dataset.attrs[attr] = CONFIG[attr] if attr in CONFIG else defaults[attr]
            else:
                print(
                    f"{attr} must be included as a global attribute, specified "
                    f"in config.yml."
                )
                # raise AttributeError(
                #     f"{attr} must be included as a global attribute, specified "
                #     f"in config.yml."
                # )


# def missing_coords(dataset: xr.Dataset, 
#                    n_dims: int, 
#                    time_var: str,
#                    options: dict = None,
#                    dim: str = None) -> xr.Dataset:
#     '''
#     Looks for any dimensions that currently don't also exist as coordinate 
#     variables.
#     dataset:        dataset to be processed.
#     n_dims:         number of spatial dimensions in dataset.
#     time_var:       time coordinate variable of dataset.
#     options:        dictionary containing key:value pairs required to populate 
#                     missing coordinates.
#     dim (optional): dimension from which coordinate should be created.
#     '''
#     # Create list of each dim not in coords;
#     # Find dtype of existing (spatial) coords (need to deal with conflicts?)
#     dim_type = np.dtype('float64')  # Default value, in case none is found
#     missing = []
#     for d in dataset.dims:
#         # Set type to match existing spatial coordinates (assume for now all 
#         # are the same)
#         if d in dataset.coords:
#             if d != time_var:
#                 dim_type = dataset.coords[d].dtype
#         elif d not in OPTIONS_DATABASE['dimensions']:
#             missing.append(d)
        
#     # If dim specified,
#     if dim and dim in dataset.dims:
#         # convert to a list of length one
#         missing = [dim]
#         if dim in dataset.coords:
#             # If coordinate variable already exists, match its type and
#             # add any new attributes.
#             dim_type = dataset.coords[dim].dtype
    
#     # else:
#     #     # Add any coordinates listed in config file that weren't already 
#     #     # found in missing coordinates
#     #     [missing.append(d) for d in CONFIG['new_coordinate_variables'].keys() if d not in missing]  # Remove this: it causes coordinates to be added when they are not present in any variables' dimensions
        
#     # For each coord in list:
#     for d in missing:
#         # Look for any required parameters & attributes in
#             # 1. CONFIG
#             # 2. VOCAB
#         attributes = {}

#         # Look in CONFIG for parameters required to set up variable:
#         # grid spacing & position (edge/centre).
#         if d in CONFIG['new_coordinate_variables'].keys():
#             config = CONFIG['new_coordinate_variables'][d]
#             try:
#                 spacing = config['spacing'] if isinstance(config['spacing'], (int, float)) else options[config['spacing']]
#                 midpoint = 'cent' in config['position'] or 'mid' in config['position']
#                 for k, v in config['attributes'].items():
#                     if k.lower() == 'units' and Units(v).isvalid:
#                         attributes['units'] = format_units(v) 
#                     else:
#                         attributes[k.lower().replace(' ', '_')] = v if k != 'standard_name' else v.replace(' ', '_')
#             except KeyError:
#                 raise KeyError(
#                     f"Processing {n_dims}d files: Parameters needed to set up "
#                     f"{d} as coordinate variable are missing from config file, "
#                     f"or required key not found in options_database. "
#                     f"For each new coordinate variable, spacing, position and "
#                     f"attributes must be specified."
#                 )
#         else:
#             raise KeyError(
#                     f"Processing {n_dims}d files: Parameters needed to set up "
#                     f"{d} as coordinate variable are missing from config file. "
#                     f"For each new coordinate variable, spacing, position and "
#                     f"attributes must be specified. Spacing must be either a "
#                     f"number or a parameter in the options_database."
#             )
        
#         # Skip looking for attributes in vocabulary, as these will be updated
#         # in cfize_variables anyway.
#         # # Look in vocabulary for any attributes.
#         # # Allow these values to override any from CONFIG, if they are valid.
#         # if d in vocabulary[n_dims].keys():
#         #     for k, v in vocabulary[n_dims][d].items():
#         #         if k.lower() == 'axis':
#         #             attributes.update({k.lower():v.upper()})
#         #         elif k.lower() == 'units' and Units(v).isvalid:
#         #             attributes.update({k.lower(): format_units(v)})
#         #         elif k.lower() == 'standard_name' or k.lower() == 'standard name':
#         #             attributes.update({'standard_name': v.replace(' ', '_')})
#         #         elif k.lower() == 'long_name' or k.lower() == 'long name':
#         #             attributes.update({'long_name': v})
        
#         # Generate coordinate points as np.ndarray of required dtype
#         points = generate_coords(number=dataset.dims[d],
#                                  spacing=spacing,
#                                  midpoint=midpoint,
#                                  data_type=dim_type)
        
#         # Assign dim, attrs & coord points to new xarray.Variable
#         new_var = xr.Variable(dims=d, attrs=attributes, data=points)

#         # Assign new variable to dataset -> overwrite existing.
#         dataset = dataset.assign(variables={dim: new_var})

#     # Return updated dataset
#     return dataset


# def cfize_variables(dataset: xr.Dataset, 
#                     n_dims: int, 
#                     time_var: str, 
#                     options: dict, 
#                     time_units: TimeUnits,
#                     variables: list|set|tuple) -> xr.Dataset:
    
#     # For each variable in args:
#     for var in variables:
#         if var not in dataset.variables:
#             raise AttributeError(
#                 f"cfize_variables: Variable {var} not found in dataset."
#             )
#             # ENHANCEMENT: If variable not in dataset, as is, look for 
#             # wildcards or use it as a stem. E.g. time_series* should identify 
#             # time_series_300_1800 as a match.

#         # Find variable's updates in VOCAB[dim][variable].
#         try:
#             updates = vocabulary[n_dims][var]
#         except KeyError:
#             # TODO: once vocab complete, make this an exception.
#             print(
#                 f"cfize_variables: {var} not found in vocabulary (dimension "
#                 f"{n_dims})."
#             )
#             # TODO: update units anyway
#             continue

#         # update variable's dimensions if required
#         if 'dimension_changes' in updates:
#             # Check specified dim(s) exist(s)
#             if not all([
#                 k in dataset[var].dims
#                 for k in updates['dimension_changes'].keys()
#             ]):
#                 raise KeyError(
#                     f"cfize_variables: {k}, specified in vocabulary - "
#                     f"dimension_changes, is not a dimension of {var}."
#                 )
#             # swap_dims
#             dataset[var] = dataset[var].swap_dims(updates['dimension_changes'])

#             # If new dim not in coords, call missing_coords function, perhaps 
#             # with new coord name as argument
#             for d in set(dataset[var].dims) - set(dataset[var].dims).intersection(set(dataset.coords)):
#                 dataset = missing_coords(
#                     dataset=dataset,
#                     n_dims=n_dims,
#                     time_var=time_var,
#                     options=options,
#                     dim=d
#                 )
            
#         # Add/update names:
#         # standard_name &/or long_name
#         if 'standard_name' in updates:
#             # TODO: check standard_name is in list of standard names.
#             # This should be done in setting up vocabulary.yml, but best to
#             # check here nonetheless.

#             dataset[var].attrs['standard_name'] = updates['standard_name']

#         if 'long_name' in updates:
#             dataset[var].attrs['long_name'] = updates['long_name']

#         # check at least one present
#         if 'standard_name' not in dataset[var].attrs and 'long_name' not in dataset[var].attrs:
#             raise KeyError("CF requires at least one of standard_name and "
#                            "long_name to be assigned to each variable.")
        
#         # Add/update units
#         if 'units' not in updates:
#             # If no new unit, check existing is present & CF compliant:
#             if 'units' in dataset[var].attrs:
#                 # If standard_name, look up
#                 if 'standard_name' in dataset[var].attrs:
#                     # TODO: look up standard names dictionary
#                     pass
#                 # Otherwise, check with cfunits.Units().isvalid
#                 else:
#                     if Units(dataset[var].attrs['units']).isvalid:
#                         dataset[var].attrs['units'] = format_units(
#                             dataset[var].attrs['units']
#                         )
#                     else:
#                         raise ValueError(
#                             f"Invalid units for variable {var} in dataset with "
#                             f"title {dataset.attrs['title']}, and not new "
#                             f"units specified in vocabulary under dimension "
#                             f"{n_dims}."
#                         )
#         else:
#             # Units specified in vocabulary.
#             # if variable is a time coordinate:
#             if var == time_var:
#                 # If supplied, process new units via TimeUnits class.
#                 # Build any missing attributes from time_units parameter.
#                 if 'calendar' in updates :
#                     calendar = updates['calendar'] 
#                 elif 'calendar' in dataset[var].attrs and Units(
#                     calendar=dataset[var].attrs['calendar']).isvalid:
#                     calendar = dataset[var].attrs['calendar']
#                 else:
#                     calendar = time_units.calendar
#                 dataset[var].attrs['calendar'] = calendar
#                 if not TimeUnits(units=updates['units']).isreftime:
#                     u = f"{updates['units']} since {time_units.reftime.isoformat()}"
#                 new_units = TimeUnits(units=u, 
#                                       calendar=calendar)
#                 # Assign T to axis.
#                 dataset[var].attrs['axis'] = 'T'

#             else:
#                 # Check any new unit specified is valid
#                 if Units(updates['units']).isvalid:
#                     new_units = Units(updates['units'])
#                 else:
#                     raise ValueError(
#                         f"Invalid units specified in vocabulary (dim {n_dims}) "
#                         f"for variable {var}."
#                     )

#             # If any new unit equivalent but not equal to any existing, run 
#             # conversion.
#             if 'units' in dataset[var].attrs:
#                 old_units = Units(dataset[var].attrs['units'])
#                 if new_units.equivalent(old_units) and not new_units.equals(old_units):
#                     Units.conform(
#                         x=dataset[var].data, 
#                         from_units=old_units, 
#                         to_units=new_units, 
#                         inplace=True
#                     )

#             # Apply new units if present, or update format of existing.
#             # Don't use cfunits.Units.formatted() method, because its ordering
#             # is weird; instead, call a function that recognises '/' as 
#             # inverting any exponent (including implied 1), interprets '^' as
#             # exponent operator
#             dataset[var].attrs['units'] = format_units(new_units)

#         # If spatial coordinate variable, add/update any specified or implied 
#         # axis.

#             # Prompt user if required.

#             # Add positive attribute if required/specified (not expected).

#         # Update variable name if specified.
            
#             # '''
#             # If updated variable is the time variable, need to pass back 
#             # updated name, or probe for this change after return.
#             # '''
        
#         # Convert any perturbation variables to absolute values.

#             # Look up reference variable's dimension
#             # ref_group = DIM_GROUPS[reference_vars[var]]

#             # # Open dataset of corresponding group
#             # while not ref_group.processed:  # op.exists(ref_ds_path):
#             #        sleep(1)
            
#             # '''
#             # This assumes required variable is in a merged dataset. Will need 
#             # more rigorous mapping if this cannot be assumed.
#             # '''

#             # Pull out reference variable as DataArray.

#             # Add reference DataArray to perturbation variable DataArray, & 
#             # re-assign to variable as new array with existing attributes etc.
#             # absolute = self.ds[variable] + ref_var
#             # self.ds[variable] = xr.DataArray(
#             #     data=absolute.data,
#             #     coords=self.ds[variable].coords,
#             #     dims=self.ds[variable].dims,
#             #     name=self.ds[variable].name,
#             #     attrs=self.ds[variable].attrs,
#             #     indexes=self.ds[variable].indexes
#             # )

#     # Return updated dataset.

    
