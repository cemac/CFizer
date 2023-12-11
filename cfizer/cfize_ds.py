import xarray as xr
from cfizer.units import TimeUnits, format_units
from cfizer.startup import *
from cfizer.utils import type_from_str, generate_coords  # , performance_time
from cfunits import Units
from dask.delayed import Delayed
from time import perf_counter, strftime, localtime
from datetime import datetime
import numpy as np
from typing import Union
import os


def get_n_dims(dataset: xr.Dataset) -> int:
    if not g:
        raise NameError(
            f"Process {os.getpid()}: get_time_var: Global dictionary g not "
            f"available."
        )

    return max(
        [
            len(dataset[v].dims)
            for v in dataset.variables
            if v != g["CONFIG"]["options_database"]["variable"]
        ]
    )


def is_monc(dataset: xr.Dataset) -> bool:
    if not g:
        raise NameError(
            f"Process {os.getpid()}: get_time_var: Global dictionary g not "
            f"available."
        )

    return g["CONFIG"]["monc_id_attribute"] in set(dataset.attrs).union(
        set(dataset.variables)
    )


def split_ds(
    dataset: xr.Dataset, shared: dict, var: str = "time"
) -> tuple[list[xr.Dataset], str]:
    
    start_time = perf_counter()
    log = []
    if var not in dataset.dims:
        raise AttributeError(f"split_ds: {var} not in dataset dimensions.")
    if shared["verbose"]:
        log.append(
            f"{strftime('%H:%M:%S', localtime())} Process {os.getpid()}: "
            f"splitting dataset {dataset.attrs['title']} by {var}."
        )
    base_title = dataset.attrs["title"].strip("_ +,.&") + "_"
    # Separate along the specified axis using xarray.Dataset.groupby method.
    # Use deep copy option to produce independent new datasets.
    grouped = {point: ds.copy(deep=True) for (point, ds) in dataset.groupby(var)}
    for point in grouped.keys():
        # Append relevant coordinate value to title
        grouped[point].attrs["title"] = base_title + str(int(point))
        if shared["verbose"]:
            log.append(
                f"{strftime('%H:%M:%S', localtime())} Process {os.getpid()}: "
                f"Created new dataset with "
                f"title, {grouped[point].attrs['title']}"
            )
    split = list(grouped.values())
    log.append(
        f"         Process {os.getpid()}: split_ds took "
        f"{perf_counter() - start_time} seconds."
    )
    return (split, log)


def ds_to_nc(
    ds: xr.Dataset,
    filepath: str,
    encodings: dict = None,
    compress: bool = False,
    shared: dict = None,
) -> None:
    if compress:
        ds.to_netcdf(
            path=filepath,
            encoding=encodings,
            engine="netcdf4",
            format="netCDF4",
        )
    else:
        ds.to_netcdf(path=filepath, encoding=encodings)


def ds_to_nc_dask(
    ds: xr.Dataset,
    filepath: str,
    encodings: dict = None,
    compress: bool = False,
    shared: dict = None,
) -> Delayed:
    if compress:
        return ds.to_netcdf(
            path=filepath,
            encoding=encodings,
            engine="netcdf4",
            format="netCDF4",
            compute=False,
        )
    else:
        return ds.to_netcdf(path=filepath, encoding=encodings, compute=False)


def perform_write(writer: Delayed, shared: dict = None) -> None:
    writer.compute()


class MoncDs:
    def __init__(
        self,
        dataset: xr.Dataset,
        time_variable: str = None,
        n_dims: int = None,
        title: str = None,
    ) -> None:
        self.log = []
        self.warnings = []
        self.ds = dataset
        self.options_db = None  # Not essential, but considered good practice

        if time_variable:
            self.time_var = time_variable
        else:
            raise AttributeError(
                "MoncDs: time variable must be specified via the "
                "time_variable parameter."
            )

        if n_dims is not None:
            self.n_dims = n_dims
        else:
            # Find number of dimensions from dataset
            self.n_dims = get_n_dims(self.ds) - 1  # Subtract 1 for time.

        if title:
            self.ds.attrs["title"] = title
        elif "title" not in self.ds.attrs:
            raise AttributeError("MoncDs: dataset title is missing.")

    def cfize(self, shared: dict):
        """
        Run CF compliance functions on dataset.

        shared: takes place of a global dictionary, to give access to 
        vocabulary. Can't use globals directly, because it creates conflicts in 
        parallel processing.
        """

        if shared["verbose"]:
            start_time = perf_counter()
            self.log.append(
                f"{strftime('%H:%M:%S', localtime())} Process {os.getpid()}: "
                f"cfize running on MONC dataset {self.ds.attrs['title']} with "
                f"{self.n_dims} spatial dimensions."
            )

        # Pull required options_database info into dictionary
        self.get_options(shared=shared)  # Sets self.options_db

        # Convert any required options_database pairs to attributes (with
        # appropriate np.dtypes)
        self.ds.attrs.update(self.options_db)

        # Add global attributes recommended by CF (some values set at
        # initialisation)
        self.add_cf_attrs(shared=shared)

        # Add any missing coordinates
        self.missing_coords(shared=shared)

        # update all (non-options-database) variables in dataset (this will
        # look for required data in shared['vocabulary'][dim], but will only
        # use what corresponds to variables present):
        var_list = set(self.ds.variables)
        var_list.discard(shared["CONFIG"]["options_database"]["variable"])
        [
            var_list.discard(v)
            for v in [
                k.replace(" ", "_").replace("-", "_")
                for k in shared["CONFIG"]["global_to_variable"].keys()
            ]
        ]  # Replacement needed to account for updated name.
        self.cfize_variables(variables=var_list, shared=shared)

        if shared["verbose"]:
            self.log.append(
                f"         Process {os.getpid()}: MoncDs.cfize took "
                f"{perf_counter() - start_time} seconds."
            )

        # Return processed dataset
        return self.ds

    def missing_coords(self, shared: dict, dim: str = None):
        """
        Looks for any dimensions that currently don't also exist as coordinate
        variables, adding if necessary.
        dim: str (optional)
            dimension from which coordinate should be created.
        """
        # Create list of each dim not in coords;
        # Find dtype of existing (spatial) coords (need to deal with conflicts?)
        dim_type = np.dtype("float64")  # Default value, in case none is found
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
                # Set type to match existing spatial coordinates (assume for 
                # now all are the same)
                if d in self.ds.coords:
                    if d != self.time_var:
                        dim_type = self.ds.coords[d].dtype
                # If dimension is not an existing coordinate variable, and not
                # one associated with the options database "variable", add to
                # the list of missing coordinates.
                elif d not in shared["CONFIG"]["options_database"]["dimensions"]:
                    missing.append(d)

        # For each coord in list:
        for d in missing:
            # Look for any required parameters & attributes in
            # shared["CONFIG"]. Although some may be defined in vocabulary,
            # this won't have the grid spacing and other needed data.
            attributes = {}

            # Look in g['CONFIG'] for parameters required to set up variable:
            # grid spacing & position (edge/centre).
            if d in shared["CONFIG"]["new_coordinate_variables"].keys():
                config = shared["CONFIG"]["new_coordinate_variables"][d]
                try:
                    spacing = (
                        config["spacing"]
                        if isinstance(config["spacing"], (int, float))
                        else self.options_db[config["spacing"]]
                    )
                    midpoint = (
                        "cent" in config["position"] or "mid" in config["position"]
                    )  # Accepts either spelling of 'centre', or mid, middle, mid-point, etc. Anything else treated as 'edge'.
                    if "attributes" in config:
                        for k, v in config["attributes"].items():
                            if k.lower() == "units" and Units(v).isvalid:
                                attributes["units"] = format_units(v)
                            else:
                                attributes[k.lower().replace(" ", "_")] = (
                                    v if k != "standard_name" else v.replace(" ", "_")
                                )
                        else:
                            # If attributes in vocabulary, skip for now, as 
                            # these will be updated in cfize_variables anyway.
                            if d not in shared["vocabulary"][self.n_dims]:
                                raise ConfigError(
                                    f"Processing {self.n_dims}d files: "
                                    f"dimension {d} has no corresponding "
                                    f"coordinate variable, and CF metadata for "
                                    f"this variable is missing from both "
                                    f"config.yml and vocabulary.yml."
                                )
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

            # Generate coordinate points as np.ndarray of required dtype
            points = generate_coords(
                number=self.ds.dims[d],
                spacing=spacing,
                midpoint=midpoint,
                data_type=dim_type,
            )

            # Assign dim, attrs & coord points to new xarray.Variable
            new_var = xr.Variable(dims=d, attrs=attributes, data=points)

            # Assign new variable to dataset -> overwrite existing.
            self.ds = self.ds.assign(variables={d: new_var})

    def update_units(self, var: str, shared: dict, updates: dict = None):
        """Add/update/check units for specified variable."""

        time_units = shared["time_units"]
        # Missing unit data may be found outside vocabulary if var is a new
        # coordinate variable.
        if var in shared["CONFIG"]["new_coordinate_variables"]:
            if not updates:
                updates = shared["CONFIG"]["new_coordinate_variables"][var][
                    "attributes"
                ]
            else:
                # Get any info missing from updates out of
                # shared['CONFIG']['new_coordinate_variables'][var]
                [
                    updates.update({k: v})
                    for k, v in shared["CONFIG"]["new_coordinate_variables"][var][
                        "attributes"
                    ].items()
                    if k not in updates
                ]

        if not updates or "units" not in updates:
            # If no new units spec, check existing is present & CF compliant:
            if "units" in self.ds[var].attrs:
                # If standard_name, look up
                if "standard_name" in self.ds[var].attrs:
                    # TODO: look up standard names dictionary
                    pass
                # Otherwise, check with cfunits.Units().isvalid
                else:
                    if Units(self.ds[var].attrs["units"]).isvalid:
                        # If unit is valid, check whether it's a time-related
                        # unit, in which case it needs to have a reference time.
                        if (
                            var == self.time_var
                            or (
                                "axis" in self.ds[var].attrs
                                and self.ds[var].attrs["axis"] == "T"
                            )
                        ) and not shared["time_units"]:
                            if not Units(self.ds[var].attrs["units"]).isreftime:
                                raise VocabError(
                                    f"update_units: "
                                    f"{self.ds[var].attrs['units']} is an "
                                    f"invalid unit for {var}: the CF format "
                                    f"for time coordinate variables is "
                                    f"'<time unit> since <reference date[time]"
                                    f"'. No updated unit found in vocabulary."
                                )
                        # If unit is valid, format it to look like units in CF
                        # standard names table.
                        self.ds[var].attrs["units"] = format_units(
                            self.ds[var].attrs["units"]
                        )
                    elif self.ds[var].attrs["units"].lower() == "fraction":
                        self.ds[var].attrs["units"] = "1"
                    else:
                        raise VocabError(
                            f"update_units: {self.ds[var].attrs['units']}: "
                            f"invalid units for variable {var} in dataset "
                            f"with title {self.ds.attrs['title']}, and no "
                            f"new units specified in vocabulary under "
                            f"dimension {self.n_dims}."
                        )
            else:
                raise VocabError(
                    f"update_units: CF conventions require that units be "
                    f"specified for each variable. No units specified in "
                    f"vocabulary ({self.n_dims}d) for {var}, and no existing "
                    f"units found."
                )
        else:
            # Units specified in vocabulary.
            # Check any new unit specified is valid
            if not Units(updates["units"]).isvalid:
                raise VocabError(
                    f"update_units: Invalid units specified in vocabulary "
                    f"for ({self.n_dims}d) variable {var}."
                )

            # if variable is a time coordinate:
            if var == (
                self.time_var
                or ("axis" in self.ds[var].attrs and self.ds[var].attrs["axis"] == "T")
                or ("axis" in updates and updates["axis"] == "T")
            ):
                # If supplied, process new units via TimeUnits class.
                # Build any missing attributes from time_units parameter.
                if "calendar" in updates:
                    calendar = updates["calendar"]
                elif (
                    "calendar" in self.ds[var].attrs
                    and Units(calendar=self.ds[var].attrs["calendar"]).isvalid
                ):
                    calendar = self.ds[var].attrs["calendar"]
                else:
                    calendar = time_units.calendar
                self.ds[var].attrs["calendar"] = calendar
                if TimeUnits(units=updates["units"]).isreftime:
                    u = updates["units"]
                else:
                    if time_units:
                        u = f"{updates['units']} since {time_units.reftime.isoformat(sep=' ')}"
                    elif (
                        "units" in self.ds[var].attrs
                        and Units(self.ds[var].attrs["units"]).isreftime
                    ):
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
                new_units = TimeUnits(units=u, calendar=calendar)
                # Assign T to axis.
                self.ds[var].attrs["axis"] = "T"

            else:
                # Units for dimension other than time
                new_units = Units(updates["units"])

                # If spatial coordinate variable, add/update any specified or
                # implied axis. Log assumption if implied.
                if var in self.ds.dims:
                    if "axis" in updates:
                        self.ds[var].attrs["axis"] = updates["axis"]
                    elif (
                        var in shared["CONFIG"]["new_coordinate_variables"]
                        and "axis"
                        in shared["CONFIG"]["new_coordinate_variables"][var][
                            "attributes"
                        ]
                    ):
                        self.ds[var].attrs["axis"] = shared["CONFIG"][
                            "new_coordinate_variables"
                        ][var]["attributes"]["axis"]
                    elif var[0].lower() in {"x", "y", "z"}:
                        if not shared["quiet"]:
                            self.log.append(
                                f"{strftime('%H:%M:%S', localtime())} MoncDs."
                                f"update_units: assumed axis "
                                f"{var[0].upper()} for variable {var}."
                            )
                        self.ds[var].attrs["axis"] = var[0].upper()
                    elif "axis" not in self.ds[var].attrs:
                        # TODO: Prompt user if required.
                        raise ConfigError(
                            f"update_units: Axis should be specified for "
                            f"spatial coordinate variables. No axis attribute "
                            f"found in dataset, vocabulary or configuration "
                            f"file for variable {var}."
                        )

                    # TODO: Add positive attribute if required/specified 
                    # (not expected).

            # If any new unit equivalent but not equal to the existing, convert.
            if "units" in self.ds[var].attrs:
                old_units = Units(self.ds[var].attrs["units"])
                if new_units.equivalent(old_units) and not new_units.equals(old_units):
                    Units.conform(
                        x=self.ds[var].data,
                        from_units=old_units,
                        to_units=new_units,
                        inplace=True,
                    )

            # Apply new units if present, or update format of existing.
            # Don't use cfunits.Units.formatted() method, because its ordering
            # is weird; instead, call a function that recognises '/' as
            # inverting any exponent (including implied 1), interprets '^' as
            # exponent operator
            self.ds[var].attrs["units"] = format_units(new_units)

    def cfize_variables(self, variables: Union[list, set, tuple], shared: dict):
        if shared["verbose"]:
            self.log.append(
                f"{strftime('%H:%M:%S', localtime())} Process {os.getpid()}: "
                f"cfizing variables on MONC dataset {self.ds.attrs['title']} "
                f"with {self.n_dims} spatial dimensions."
            )

        # Process any coordinate variables first, as changes in these may affect
        # processing of data variables.
        coord_vars = [v for v in variables if v in self.ds.dims]
        variables = coord_vars + [v for v in variables if v not in coord_vars]

        # For each variable in args:
        for var in variables:
            if var not in self.ds.variables:
                raise AttributeError(
                    f"cfize_variables: Variable {var} not found in dataset."
                )  # This won't happen in current version, as variables argument *is* from self.ds.variables.
                # ENHANCEMENT: If variable not in dataset, as is, look for
                # wildcards or use it as a stem. E.g. time_series* and 
                # time_series should both identify time_series_300_1800 as a 
                # match.

            # Find variable's updates in shared['vocabulary'][dim][variable].
            try:
                updates = shared["vocabulary"][self.n_dims][var]
            except KeyError:
                # If not in vocabulary, check variable is already CF compliant.
                if not(
                    "standard_name" in self.ds[var].attrs or "long_name" in self.ds[var].attrs
                    ) or "units" not in self.ds[var].attrs:
                    raise VocabError(
                        f"cfize_variables: {var} not found in vocabulary "
                        f"(dimension {self.n_dims})."
                    )
                # Subsequent method calls will evaluate validity of units and, 
                # if present, standard_name.
                updates = {}  # No keys means no updates - only checking.
            
            # update variable's dimensions if required
            if "dimension_changes" in updates:
                # Check specified dim(s) exist(s)
                if not all(
                    [
                        k in self.ds[var].dims
                        for k in updates["dimension_changes"].keys()
                    ]
                ):
                    raise VocabError(
                        f"cfize_variables: At least one dimension specified in "
                        f"vocabulary - dimension_changes, is not a dimension "
                        f"of {var}."
                    )
                # Swap dimensions using xarray.DataArray.swap_dims method.
                self.ds[var] = self.ds[var].swap_dims(updates["dimension_changes"])

                # If new dim not in coords, call missing_coords function,
                # with new coord name as argument
                for d in set(self.ds[var].dims) - set(self.ds[var].dims).intersection(
                    set(self.ds.coords)
                ):
                    self.missing_coords(shared=shared, dim=d)

            # Add/update names:
            # standard_name &/or long_name
            if "standard_name" in updates:
                self.ds[var].attrs["standard_name"] = updates["standard_name"]

            if "standard_name" in self.ds[var].attrs:
                # TODO: check standard_name is in list of standard names.
                # This should be done in setting up vocabulary.yml, but best to
                # check here nonetheless. This also then works for cases where
                # standard_name is not in vocabulary entry, but is already
                # present in the variable's metadata.
                pass

            if "long_name" in updates:
                self.ds[var].attrs["long_name"] = updates["long_name"]

            # check at least one present
            if (
                "standard_name" not in self.ds[var].attrs
                and "long_name" not in self.ds[var].attrs
            ):
                raise VocabError(
                    f"cfize_variables: {var}: CF requires at least one of "
                    f"standard_name and long_name to be assigned to each "
                    f"variable."
                )

            # Add/update units
            self.update_units(var=var, updates=updates, shared=shared)

            # Update variable name if specified.
            if "updated_name" in updates:
                self.ds = self.ds.rename({var: updates["updated_name"]})
                """
                If updated variable is the time variable, need to pass 
                back updated name, or probe for this change after return.
                """
                if var == self.time_var:
                    self.time_var = updates["updated_name"]

            # Convert any perturbation variables to absolute values.
            if (
                "perturbation_to_absolute" in updates
                and updates["perturbation_to_absolute"]
            ):
                if "reference_variable" not in updates:
                    raise VocabError(
                        f"cfize_variables: {var}: If "
                        f"perturbation_to_absolute is True, "
                        f"reference_variable must contain the name of the "
                        f"variable containing reference value(s)."
                    )

                # Get reference variable data array from (updated) vocabulary.
                # This will crash if reference data are not available, and crash
                # needs to be handled by function calling MoncDs.cfize method.
                ref_array = updates["ref_array"]

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
                )

    def add_cf_attrs(self, shared: dict, **kwargs):
        defaults = {}

        # History attribute should be an audit trail, and so include
        # whatever processing is done to each file, and updated just
        # prior to saving NetCDF file. Default is to at least report CFizer
        # version.
        defaults[
            "history"
        ] = f'{datetime.now().isoformat(timespec="minutes", sep=" ")}: MONC output files converted by CFizer version {VERSION}, https://github.com/cemac/CFizer.'
        defaults["Conventions"] = f"CF-{CF_VERSION}"
        for attr in CF_ATTRIBUTES:
            # Only overwrite existing value if a new value was specified in
            # config file.
            if attr in self.ds.attrs:
                if attr in shared["CONFIG"] and shared["CONFIG"][attr] is not None:
                    self.ds.attrs[attr] = shared["CONFIG"][attr]
            elif (attr in shared["CONFIG"] and shared["CONFIG"][attr]) or (
                attr in defaults and defaults[attr]
            ):
                self.ds.attrs[attr] = (
                    shared["CONFIG"][attr]
                    if attr in shared["CONFIG"]
                    else defaults[attr]
                )
            else:
                self.warnings.append(
                    ConfigWarning(
                        f"add_cf_attrs: CF-1.10, section 2.6.2, recommends {attr} be included as a global attribute. This can be specified in config.yml."
                    )
                )

    def get_options(self, shared: dict) -> dict:
        """
        Note: the xarray.Dataset.assign_attrs() method isn't suitable here, as 
        it creates a new dataset, rather than updating the existing one.
        """

        fields = shared["CONFIG"]["options_to_attrs"]
        # Get required fields from argument; import all if none supplied.
        if fields is None:
            if self.ds[shared["CONFIG"]["options_database"]["variable"]].dtype == "S1":
                # This is used if concat_characters=False options is used in
                # xarray.open_dataset. It seems to create some problems down
                # the line, so is not advised.
                options = {
                    "".join([c.decode("utf-8") for c in k]): type_from_str(
                        "".join([c.decode("utf-8") for c in v])
                    )
                    for k, v in self.ds[
                        shared["CONFIG"]["options_database"]["variable"]
                    ].data
                }  # If dataset opened with concat_characters=False
            else:
                options = {
                    k.decode("utf-8"): type_from_str(v.decode("utf-8"))
                    for [k, v] in self.ds[
                        shared["CONFIG"]["options_database"]["variable"]
                    ].data
                }
        else:
            if self.ds[shared["CONFIG"]["options_database"]["variable"]].dtype == "S1":
                options = {
                    "".join([c.decode("utf-8") for c in k]): type_from_str(
                        "".join([c.decode("utf-8") for c in v])
                    )
                    for k, v in self.ds[
                        shared["CONFIG"]["options_database"]["variable"]
                    ].data
                    if "".join([c.decode("utf-8") for c in k]) in fields
                }  # If dataset opened with concat_characters=False
            else:
                options = {
                    k.decode("utf-8"): type_from_str(v.decode("utf-8"))
                    for [k, v] in self.ds[
                        shared["CONFIG"]["options_database"]["variable"]
                    ].data
                    if k.decode("utf-8") in fields
                }

        # Drop any options that are reset by DEPHY, if used.
        if all(
            [
                opt in options and options[opt] == val
                for opt, val in shared["CONFIG"]["dephy_true_if"].items()
            ]
        ):
            for a in shared["CONFIG"]["drop_for_dephy"]:
                if a in options:
                    options.pop(a)

        self.options_db = options
