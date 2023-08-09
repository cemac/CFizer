import xarray as xr
from cfunits import Units


def cfize_variable(variable: xr.DataArray) -> xr.DataArray:
    # Find variable in vocabulary

    # Apply changes to variable where necessary

    # Check units are valid and apply standard formatting (e.g. kg/kg becomes 1)
    
    # return new version of variable's data array.
    pass


def coord_points(number: int, spacing: float|int, midpoint: bool = False) -> list:
    pass
