import xarray as xr
from algorithm import TimeUnits
from time import sleep
import os.path as op
from setup import *
from utils import format_units


def cfize_dataset(dataset: xr.Dataset,
                  dimension: int,
                  title: str,
                  time_units: TimeUnits) -> xr.Dataset:
    
    # Pull required options_database info into dictionary

    # Add global attributes recommended by CF (some values set at 
    # initialisation)

    # Convert any required options to attributes (with appropriate np.dtypes)

    # Add any missing coordinates
        # dataset = missing_coords(dataset)

    # update all variables in either VOCAB[dim] or dataset:
        # dataset = cfize_variables(dataset, *var_list) to 

    # Return updated dataset.
    pass


def missing_coords(dataset: xr.Dataset, dim: str = None) -> xr.Dataset:
    # Find dtype of existing (spatial) coords (need to deal with conflicts?)

    # If dim specified,
        # convert to a list of length one

    # else:
        # Create list of each dim not in coords, or of specified new coords:

    # For each coord in list:
        # Look for any required attributes
            # 1. VOCAB
            # 2. CONFIG
        
        # Get spacing options variable & edge/centre from CONFIG.

        # Get grid spacing

        # Generate coordinate points as np.ndarray of required dtype

        # Assign dim, attrs & coord points to new xarray.Variable

        # Assign new variable to dataset -> overwrite existing.

    # Return updated dataset
    
    pass


def cfize_variables(dataset: xr.Dataset, *args: str) -> xr.Dataset:
    pass
    # For each variable in args:

        # If variable not in dataset, as is, look for wildcards or use it as a 
        # stem. E.g. time_series* should identify time_series_300_1800 as a
        # match.
    
        # Find variable's updates in VOCAB[dim][variable].

        # update variable's dimensions if required

            # Check specified dim(s) exist(s)

            # swap_dims

            # If new dim not in coords, call missing_coords function, perhaps 
            # with new coord name as argument

        # Add/update names:
            # standard_name &/or long_name

            # check at least one present
        
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

    
