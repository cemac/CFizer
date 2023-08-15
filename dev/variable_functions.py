import xarray as xr
from cfunits import Units


def cf_compliant(variable: xr.DataArray) -> bool:
    '''
    Presence of axis attribute not enforced for coordinates, but advisory given if absent.
    Time does not require axis=T, assuming standard_name is appropriate, so this is not checked here.
    Does not check for positive attribute, which is needed when direction of coordinate variable is ambiguous.
    '''
    # Check name is defined
    if not ('standard_name' in variable.attrs or 'long_name' in variable.attrs):
        return False
    # Check valid units
    if 'units' not in variable.attrs or not Units(units=variable.attrs['units']).isvalid:
        return False
    # If coordinate variable, check for axis
    if variable.name in variable.coords and variable.name[0].lower() in {
         'x', 'y', 'z'}:
            if 'axis' not in variable.attrs:
                print(f'ADVISORY: `{variable.name}`: CF Conventions recommend adding an axis attribute for spatial coordinates.')
    return True


def cfize_variable(variable: xr.DataArray) -> xr.DataArray:
    # Find variable in vocabulary

    # Apply changes to variable where necessary

    # Check units are valid and apply standard formatting (e.g. kg/kg becomes 1)
    
    # return new version of variable's data array.
    pass


def test():
    pass

if __name__ == "__main__": test()
