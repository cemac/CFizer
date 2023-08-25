import xarray as xr
from algorithm import TimeUnits
from time import sleep
import os.path as op
from setup import *
from utils import type_from_str, format_units, generate_coords
from cfunits import Units


def cfize_dataset(dataset: xr.Dataset,
                  n_dims: int,
                  title: str,
                  time_var: str,
                  time_units: TimeUnits) -> xr.Dataset:
    
    # Pull required options_database info into dictionary
    options_db = get_options(dataset=dataset, fields=CONFIG['options_to_attrs'])

    # Convert any required options_database pairs to attributes (with 
    # appropriate np.dtypes)
    dataset.attrs.update(options_db)

    # Add global attributes recommended by CF (some values set at 
    # initialisation)
    add_cf_attrs(dataset=dataset)
    
    # Add any missing coordinates
    dataset = missing_coords(dataset=dataset,
                             n_dims=n_dims,
                             time_var=time_var,
                             options=options_db)

    # update all variables in either VOCAB[dim] or dataset:
    var_list = set(dataset.variables)
    var_list.discard(OPTIONS_DATABASE['variable'])
    [var_list.discard(v) for v in CONFIG['global_to_variable'].keys()]
    dataset = cfize_variables(dataset=dataset,
                              n_dims=n_dims,
                              time_var=time_var,
                              options=options_db,
                              variables=var_list)

    # Return updated dataset.
    pass


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
        defaults = {attr: None for attr in CF_ATTRIBUTES}
        # defaults['title'] = kwargs['title'] if 'title' in kwargs
        defaults['history'] = f'{datetime.now().isoformat(timespec="minutes")}: output files processed using CFizer version <version>, <url>.'
        defaults['conventions'] = f"CF-{CF_VERSION}"
        for attr in CF_ATTRIBUTES:
            # Only overwrite existing value if a new value was specified in
            # config file.
            if attr in dataset.attrs :
                if attr in CONFIG and CONFIG[attr] is not None:
                    dataset.attrs[attr] = CONFIG[attr]
            else:
                dataset.attrs[attr] = CONFIG[attr] if attr in CONFIG else defaults[attr]


def missing_coords(dataset: xr.Dataset, 
                   n_dims: int, 
                   time_var: str,
                   options: dict = None,
                   dim: str = None) -> xr.Dataset:
    '''
    Looks for any dimensions that currently don't also exist as coordinate 
    variables.
    dataset:        dataset to be processed.
    n_dims:         number of spatial dimensions in dataset.
    time_var:       time coordinate variable of dataset.
    options:        dictionary containing key:value pairs required to populate 
                    missing coordinates.
    dim (optional): dimension from which coordinate should be created.
    '''
    # Create list of each dim not in coords;
    # Find dtype of existing (spatial) coords (need to deal with conflicts?)
    dim_type = np.dtype('float64')  # Default value, in case none is found
    missing = []
    for d in dataset.dims:
        # Set type to match existing spatial coordinates (assume for now all 
        # are the same)
        if d in dataset.coords:
            if d != time_var:
                dim_type = dataset.coords[d].dtype
        elif d not in OPTIONS_DATABASE['dimensions']:
            missing.append(d)
        
    # If dim specified,
    if dim and dim in dataset.dims:
        # convert to a list of length one
        missing = [dim]
        if dim in dataset.coords:
            # If coordinate variable already exists, match its type and
            # add any new attributes.
            dim_type = dataset.coords[dim].dtype
    
    # else:
    #     # Add any coordinates listed in config file that weren't already 
    #     # found in missing coordinates
    #     [missing.append(d) for d in CONFIG['new_coordinate_variables'].keys() if d not in missing]  # Remove this: it causes coordinates to be added when they are not present in any variables' dimensions
        
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
                spacing = config['spacing'] if isinstance(config['spacing'], (int, float)) else options[config['spacing']]
                midpoint = 'cent' in config['position'] or 'mid' in config['position']
                for k, v in config['attributes'].items():
                    if k.lower() == 'units' and Units(v).isvalid:
                        attributes['units'] = format_units(v) 
                    else:
                        attributes[k.lower().replace(' ', '_')] = v if k != 'standard_name' else v.replace(' ', '_')
            except KeyError:
                raise KeyError(
                    f"Processing {n_dims}d files: Parameters needed to set up "
                    f"{d} as coordinate variable are missing from config file, "
                    f"or required key not found in options_database. "
                    f"For each new coordinate variable, spacing, position and "
                    f"attributes must be specified."
                )
        else:
            raise KeyError(
                    f"Processing {n_dims}d files: Parameters needed to set up "
                    f"{d} as coordinate variable are missing from config file. "
                    f"For each new coordinate variable, spacing, position and "
                    f"attributes must be specified. Spacing must be either a "
                    f"number or a parameter in the options_database."
            )
        
        # Skip looking for attributes in VOCABULARY, as these will be updated
        # in cfize_variables anyway.
        # # Look in VOCABULARY for any attributes.
        # # Allow these values to override any from CONFIG, if they are valid.
        # if d in VOCABULARY[n_dims].keys():
        #     for k, v in VOCABULARY[n_dims][d].items():
        #         if k.lower() == 'axis':
        #             attributes.update({k.lower():v.upper()})
        #         elif k.lower() == 'units' and Units(v).isvalid:
        #             attributes.update({k.lower(): format_units(v)})
        #         elif k.lower() == 'standard_name' or k.lower() == 'standard name':
        #             attributes.update({'standard_name': v.replace(' ', '_')})
        #         elif k.lower() == 'long_name' or k.lower() == 'long name':
        #             attributes.update({'long_name': v})
        
        # Generate coordinate points as np.ndarray of required dtype
        points = generate_coords(number=dataset.dims[d],
                                 spacing=spacing,
                                 midpoint=midpoint,
                                 data_type=dim_type)
        
        # Assign dim, attrs & coord points to new xarray.Variable
        new_var = xr.Variable(dims=d, attrs=attributes, data=points)

        # Assign new variable to dataset -> overwrite existing.
        dataset = dataset.assign(variables={dim: new_var})

    # Return updated dataset
    return dataset


def cfize_variables(dataset: xr.Dataset, 
                    n_dims: int, 
                    time_var: str, 
                    options: dict, 
                    variables: list|set|tuple) -> xr.Dataset:
    
    # For each variable in args:
    for var in variables:
        if var not in dataset.variables:
            raise AttributeError(
                f"cfize_variables: Variable {var} not found in dataset."
            )
            # ENHANCEMENT: If variable not in dataset, as is, look for 
            # wildcards or use it as a stem. E.g. time_series* should identify 
            # time_series_300_1800 as a match.

        # TODO: if var is time_var, need to use wildcards to match to vocab
        if var == time_var:
            pass
    
        # Find variable's updates in VOCAB[dim][variable].
        try:
            updates = VOCABULARY[n_dims][var]
        except KeyError:
            # TODO: once vocab complete, make this an exception.
            print(
                f"cfize_variables: {var} not found in vocabulary (dimension "
                f"{n_dims})."
            )
            continue

        # update variable's dimensions if required
        if 'dimension_changes' in updates:
            # Check specified dim(s) exist(s)
            if not all([
                k in dataset[var].dims
                for k in updates['dimension_changes'].keys()
            ]):
                raise KeyError(
                    f"cfize_variables: {k}, specified in vocabulary - "
                    f"dimension_changes, is not a dimension of {var}."
                )
            # swap_dims
            dataset[var] = dataset[var].swap_dims(updates['dimension_changes'])

            # If new dim not in coords, call missing_coords function, perhaps 
            # with new coord name as argument
            for d in set(dataset[var].dims) - set(dataset[var].dims).intersection(set(dataset.coords)):
                dataset = missing_coords(
                    dataset=dataset,
                    n_dims=n_dims,
                    time_var=time_var,
                    options=options,
                    dim=d
                )
            
        # Add/update names:
        # standard_name &/or long_name
        if 'standard_name' in updates:
            # TODO: check standard_name is in list of standard names.
            # This should be done in setting up vocabulary.yml, but best to
            # check here nonetheless.

            dataset[var].attrs['standard_name'] = updates['standard_name']

        if 'long_name' in updates:
            dataset[var].attrs['long_name'] = updates['long_name']

        # check at least one present
        if 'standard_name' not in dataset[var].attrs and 'long_name' not in dataset[var].attrs:
            raise KeyError("CF requires at least one of standard_name and "
                           "long_name to be assigned to each variable.")
        
        # Add/update units

            # If no new unit, check existing is present & CF compliant:

                # If standard_name, look up

                # Otherwise, check with cfunits.Units().isvalid

            # if variable is a time coordinate:

                # If supplied, process new units via TimeUnits class

                    # Build any missing attributes from time_units parameter

                # Assign T to axis.

            # else:
                # Check any new unit specified is valid

            # If any new unit equivalent but not equal to any existing, run 
            # conversion.

            # Apply new units if present, or update format of existing.
            # Don't use cfunits.Units.formatted() method, because its ordering
            # is weird; instead, call a function that recognises '/' as 
            # inverting any exponent (including implied 1), interprets '^' as
            # exponent operator
            # dataset[var].attrs['units'] = format_units(units)

        # If spatial coordinate variable, add/update any specified or implied 
        # axis.

            # Prompt user if required.

            # Add positive attribute if required/specified (not expected).

        # Update variable name if specified.
            
            # '''
            # If updated variable is the time variable, need to pass back 
            # updated name, or probe for this change after return.
            # '''
        
        # Convert any perturbation variables to absolute values.

            # Look up reference variable's dimension
            # ref_group = DIM_GROUPS[reference_vars[var]]

            # # Open dataset of corresponding group
            # while not ref_group.processed:  # op.exists(ref_ds_path):
            #        sleep(1)
            
            # '''
            # This assumes required variable is in a merged dataset. Will need 
            # more rigorous mapping if this cannot be assumed.
            # '''

            # Pull out reference variable as DataArray.

            # Add reference DataArray to perturbation variable DataArray, & 
            # re-assign to variable as new array with existing attributes etc.
            # absolute = self.ds[variable] + ref_var
            # self.ds[variable] = xr.DataArray(
            #     data=absolute.data,
            #     coords=self.ds[variable].coords,
            #     dims=self.ds[variable].dims,
            #     name=self.ds[variable].name,
            #     attrs=self.ds[variable].attrs,
            #     indexes=self.ds[variable].indexes
            # )

    # Return updated dataset.

    
