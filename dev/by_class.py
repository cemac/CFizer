import xarray as xr
from setup import *
import os
from glob import glob
from cfunits import Units
from dataset_functions import *
from variable_functions import *
from difflib import SequenceMatcher
from utils import type_from_str, generate_coords, stem_str


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
                 directory: str = None,
                 target: str = None) -> None:
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
        self.sort_nc()
        # dim_groups = {
        #     name: DsGroup(name=name, action=action) 
        #     for name, action in self.DIM_ACTIONS.items()
        # }
        self.by_dim = {n_dim: DsGroup(
            name=group, action=self.DIM_ACTIONS[group]
        ) for n_dim, group in self.DIM_GROUPS}
        # {
        #     n_dims: DimGroup(group=dim_groups[group])
        #     for n_dims, group in self.DIM_GROUPS.items()
        # }
        self.time_units = None
        self.target_dir = target or f'{os.path.dirname(self.directory)}+processed'
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
            with xr.open_dataset(filepath, decode_times=False) as ds:
                if is_monc(ds):
                    self.monc_files.append(filepath)
                    # Find number of (non-options-database) dimensions, and add
                    # filepath to the appropriate dataset group.
                    n_dims = get_n_dims(ds) - 1
                    self.by_dim[n_dims].add(filepath=filepath)
                else:
                    self.input_files.append(filepath)

    def merge_by_group(self):
        to_merge = {}
        for name, action in self.DIM_ACTIONS.items():
            if action == 'merge':
                merge_group = [group for group in self.by_dim.values() if group.name == name]
                if len(merge_group) > 1:
                    to_merge[name] = merge_group
                else:
                    to_merge.pop(name)
        merged = {}
        for name, merge_list in to_merge.items():
            merged[name] = merge_ds_groups(merge_list)

    def datasets(self):
        # Generator
        return (xr.open_dataset(path, decode_times=False) for path in self.monc_files)
    
    def cfize(self):
        for dim, group in self.by_dim.items():
            if group.action == 'merge':
                # extract global attributes for each file, and assign each as a `data_var` with dimension of time.
                for ds in group.datasets():
                    ds.attr_to_var(CONFIG['global_to_variable'])

                # Merge all datasets in group to create single time-series dataset
                group.process()

                # if dim > 0:
                #     # Check whether dimension group below has same name, & therefore should be merged.
                #     if group.name == self.by_dim[dim - 1].name:
                #         # merge output dataset with that of dimension below
                #         pass
            else:
                # 
                pass
        
        # Merge any groups with the same name: these comprise different
        # dimensions that should be merged, i.e. merging the 0d and 1d
        # time series into 0+1d.
        # Collect the merged groups, and any groups not needing to be merged,
        # into self.processed dict.
        for group in self.processed.keys():
            groups_to_merge = [g for g in self.by_dim.values() if g == group]
            if len(groups_to_merge) > 1:
                # Create new DsGroup containing groups to be merged
                self.processed[group] = DsGroup(groups=groups_to_merge)
                # Perform the merge
                self.processed[group].process()
            else:
                self.processed[group] = groups_to_merge[0]

        
        
        
            
# class SubGroup:
#     def __init__(self, *args) -> None:
#         if len(args) <= 1 or not all(
#             [isinstance(arg, DsGroup) for arg in args]):
#             raise TypeError('Expecting a sequence of DsGroup objects.')
#         if not all([len(group.processed) == 1 for group in args]):
#             raise TypeError("Each DsGroup.processed attribute must contain a single dataset.")
        
#         self.ds = [group.processed for group in args]
#         # elif all([isinstance(arg, xr.MoncDataset) for arg in args]):
#         #     self.ds = [arg.ds for arg in args]
#         # elif all([isinstance(arg, xr.Dataset) for arg in args]):
#         #     self.ds = args
#         # else:
#         #     raise TypeError('Valid arguments must be DsGroup types.')
        
#         self.stem = stem_str([group.stem for group in args])
#         self.processed = None
        
#     def merge_dims(self):
#         self.processed = None




class DsGroup:
    '''
    Each file in a DsGroup has the same variables, coordinates & dimensions.
    '''
    def __init__(self, 
                 name: str = '', 
                 action: str = None, 
                 filepaths: list|set|tuple = [],
                 groups: list|set|tuple = None) -> None:
        
        # TODO: validate parameters

        if groups and all([isinstance(g, DsGroup) for g in groups]):
            # Check each group has only one dataset in its processed attribute.
            if not all([isinstance(group.processed, xr.Dataset) for group in groups]):
                raise TypeError('To merge DsGroup objects, the processed attribute of each must contain a single xarray.Dataset object.')
            # merge existing groups
            self.name = stem_str([group.name for group in groups])
            self.action = 'merge_groups'
            self.process = self.merge_groups
            self.stem = stem_str([group.stem for group in groups])
            self.filepaths = []
            for group in groups:
                self.filepaths += group.filepaths
            self.to_merge = [group.processed for group in groups]
            self.processed = None
            # self.process()
        else:
            self.name = name
            self.action = action.lower()
            if self.action == 'split':
                self.process = self.split_times
            elif self.action == 'merge':
                self.process = self.merge_times
            # elif self.action = 'merge_groups':
            #     self.process = self.merge_groups
            else:
                self.process = self.keep_time_series
            self.filepaths = filepaths
            self.processed = []
            self.stem = stem_str(self.filepaths)
            # Although it would be best to verify that all member datasets have the same time variable, for now, just use the first one.
            self.time_var = MoncDataset(filepath=self.filepaths[0]).time_var
            # TODO: check the dataset object created for self.time_var is not kept in memory after this __init__ function.

    def datasets(self):
        '''
        This returns a generator that yields the datasets contained in each
        filepath.
        '''
        # return (xr.open_dataset(path, decode_times=False) for path in self.filepaths)
        return (MoncDataset(path) for path in self.filepaths)

    def add(self, filepath):
        # TODO: validate filepath
        self.filepaths.append(filepath)
        self.stem = stem_str(filepath, self.stem)

    # def stem_str(self, *args):
    #     if not args:
    #         raise ValueError('At least one filename/path must be provided.')
    #     filenames = [os.path.splitext(os.path.basename(path))[0] for path in args]  # If each path contains only a filename (no extension) this will still work.
    #     stem = filenames[0]
    #     if len(args) > 1:
    #         for f in filenames[1:]:
    #             match = SequenceMatcher(a=stem, b=f).find_longest_match()
    #             stem = stem[match.a: match.a + match.size]
    #     return stem

    def keep_time_series(self):
        # copy all datasets/files into self.processed
        self.processed = [ds for ds in self.datasets()]

    def merge_times(self):
        '''
        Using the option, data_vars='minimal', and merging only on
        coords=[self.time_var] ensures the options_database variable is not
        replicated for every time-point.
        combine_attrs='drop_conflicts' can be used because any necessary global
        attributes can be preserved by MoncDataset.attr_to_var().
        '''
        self.processed = xr.combine_by_coords(
            [monc_ds.ds for monc_ds in self.datasets()],
            combine_attrs='drop_conflicts',
            join='exact',
            coords=[self.time_var],
            data_vars='minimal'
        )

    def split_times(self):
        split = None
        self.processed += split

    def merge_groups(self):
        self.processed = merge_dimensions(*self.to_merge)
        

# class DimGroup:
#     def __init__(self, group: DsGroup) -> None:
#         self.group = group

#     def add_file(self, n_dims: int, filepath: str):
#         if n_dims in self.groups:
#             self.groups[n_dims].add(filepath)
#         else:
#             self.groups[n_dims] = DsGroup(
#                 name=self.DIM_GROUPS[n_dims]['group'],
#                 action=self.DIM_GROUPS[n_dims]['action'],
#                 filepaths=[filepath]
#             )

#     def merge_dims(self, *args):
#         to_merge = {}
#         # Find dims to merge
#         for group, action in self.DIM_ACTIONS.items():
#             to_merge[group] = [n_dims for n_dims, name in self.DIM_GROUPS.items() if name == group and action == 'merge']
#             if len(to_merge[group]) <= 1:
#                 to_merge.pop(group)
#         print(to_merge)


class TimeUnits(Units):
    def __init__(self, units=None, calendar=None, formatted=False, names=False, definition=False, _ut_unit=None):
        super().__init__(units, calendar, formatted, names, definition, _ut_unit)
        if not self.isvalid:
            raise ValueError('Units not valid according to cfunits module.')
    
    def time_unit(self) -> str:
        return self.units.split(' since ')[0]
    
    def since(self) -> str:
        return self.reftime.isoformat()
    
    def cf(self) -> str:
        return self.formatted()
    

class OptionsDb:
    _dephy_options = ('dephy_file', 'dephy_forcings_enabled')
    _drop_for_dephy = ("longitude", "latitude", "z0")

    def __init__(self, dataset: xr.Dataset, fields: list|set|tuple = None) -> None:
        
        if not isinstance(dataset, xr.Dataset):
            raise TypeError('dataset argument must be of type xarray.Dataset')
        
        self._ds = dataset

        # Get required fields from argument; import all if none supplied.
        if fields is None:
            [self.__setattr__(
                k.decode('utf-8'), type_from_str(v.decode('utf-8'))
                ) for [k, v] in self._ds.options_database.data]
        else:
            [self.__setattr__(
                k.decode('utf-8'), type_from_str(v.decode('utf-8'))
                ) for [k, v] in self._ds.options_database.data
                if k.decode('utf-8') in fields]
        
        if (all([self.has_attr(opt) for opt in self._dephy_options]) &
            all([self.__getattribute__(opt) for opt in self._dephy_options])):
            for a in self._drop_for_dephy:
                if self.has_attr(a):
                    self.__delattr__(a)
        
    def has_attr(self, attr: str) -> bool:
        try:
            self.__getattribute__(attr)
        except AttributeError:
            return False
        else:
            return True
        
    def to_attrs(self) -> xr.Dataset:
        for k in [k for k in dir(self) if (k[0] != '_')]:
            self._ds.attrs[k] = self.__getattribute__(k)
        return self._ds
        
    # def dx(self):
    #     return self.dxx
    
    # def dy(self):
    #     return self.dyy
    
    def get(self, attr: str):
        return self.__getattribute__(attr)
    

class MoncDataset:
    ''''
    Notes:
        This has been set up as a standalone class, rather than inheriting from xr.Dataset, to avoid some complexities of inheriting from xarray data classes.
    '''
    do_not_duplicate = ('MONC_timestep', )

    def __init__(self, filepath: str = None, dataset: xr.Dataset = None):
        if dataset:
            if not isinstance(dataset, xr.Dataset):
                raise TypeError
            self.ds = dataset if is_monc(dataset) else None
        elif filepath:
            if not isinstance(filepath, str):
                raise TypeError('filepath must be a string')
            try:
                self.ds = xr.open_dataset(filepath, decode_times=False)
            except Exception as e:
                raise e
            else:
                if not is_monc(self.ds):
                    self.ds = None
        else:
            raise ValueError('Either a filepath (string) or xarray.Dataset must be passed in.')
        # data_vars = [MoncVariable(data_array=var) for var in ds.data_vars.values()]
        if self.ds:
            self.options = OptionsDb(self.ds, fields=CONFIG['options_to_attrs'])
            time_vars = [d for d in self.ds.dims if 'time' in d]
            if len(time_vars) != 1
                # Throw error or prompt user to select which dimension to use
                raise AttributeError(f'More than one time dimension found: {time_vars}')
            self.time_var = time_vars[0]
    
    def add_cf_attrs(self):
        pass

    def attr_to_var(self, attr: str|list|set|tuple = None):
        
        if attr is None:
            # convert all attributes to variables
            attr = set(self.ds.attrs.keys())
        
        if isinstance(attr, str):
            attr = {attr}
        
        last_time_point = self.ds[self.time_var].data[-1]
        globals = {}
        # if isinstance(attr, (list, set, tuple)) and all(
        #     [isinstance(a, str) for a in attr]):
        if all([isinstance(a, str) for a in attr]):
            for g in attr:
                name = g.replace(' ', '_')  # TODO: Replace this with a full "safe_string" function
                if name in self.do_not_duplicate:
                    data = type_from_str(self.ds.attrs[g])
                    coords = (last_time_point,)
                else:
                    data = [type_from_str(self.ds.attrs[g])]*len(self.ds[self.time_var].data)
                    coords = self.ds[self.time_var]
                globals[g] = xr.DataArray(
                    name=name,
                    data=data,
                    coords={self.time_var: coords},
                    dims=(self.time_var,)
                )
        else:
            raise TypeError('attr must be a string, a collection of strings, or None.')
        self.ds = self.ds.assign(
            {name: array for name, array in globals.items()}
        )

    # def var_from_global(self, attr: str) -> xr.DataArray:
    #         name = attr.replace(' ', '_')  # TODO: Replace this with a full "safe_string" function
    #         if name in self.do_not_duplicate:
    #             data = type_from_str(self.ds.attrs[attr])
    #             coords = (last_time_point,)
    #         else:
    #             data = [type_from_str(self.ds.attrs[attr])]*len(self.ds[self.time_var].data)
    #             coords = self.ds[self.time_var]
    #         return xr.DataArray(
    #             name=name,
    #             data=data,
    #             coords={self.time_var: coords},
    #             dims=(self.time_var,)
    #         )

    def missing_coords(self):
        # Look for any dimensions that currently don't also exist as coordinate variables
        missing = [dim for dim in self.ds.dims if (
            dim not in self.ds.coords &
            dim not in OPTIONS_DATABASE['dimensions'])]
        
        # Look for any coordinates listed in config file that weren't already found in missing coordinates
        [missing.append(dim) for dim in CONFIG['new_coordinate_variables'].keys() if dim not in missing]
        
        for dim in missing:
            attributes = {}
            if dim in CONFIG['new_coordinate_variables']:
                # Get required information from config file, if available
                config = CONFIG['new_coordinate_variables'][dim]
                spacing = config['spacing']
                midpoint = 'cent' in config['position'] or 'mid' in config['position']
                for k, v in config['attributes'].items():
                    if k.lower() == 'units':
                        attributes['units'] = cfunits.Units(v).formatted()
                    else:
                        attributes[k.lower().replace(' ', '_')] = v if k != 'standard_name' else v.replace(' ', '_')
            else:
                # If not in config file, assume it's x, y, xu or yu, with associated default attributes etc.
                if dim[0] == 'x':
                    spacing = self.options.get('dxx')
                    attributes = {
                        'axis': 'X',
                        'units': 'm'
                    }
                elif dim[0] == 'y':
                    spacing = self.options.get('dyy')
                    attributes = {
                        'axis': 'Y',
                        'units': 'm'
                    }
                else:
                    # TODO: prompt user for grid spacing
                    raise AttributeError(f'Unknown grid spacing for {dim}.')
                attributes['units'] = 'm'
                midpoint = dim[-1] not in ('u', 'v')
                attributes['long_name'] = f'{dim[0]}-coordinate in Cartesian system (cell-{'centers' if midpoint else 'edges'})'
                
            # Generate data points for coordinate variable
            points = generate_coords(number=self.ds.dims[dim],
                                     spacing=spacing,
                                     midpoint=midpoint)
            
            # create new xarray variable based on the new coordinate
            new_var = xr.Variable(dims=dim, attrs=attributes, data=points)

            # Add new coordinate variable to dataset
            self.ds = ds.assign(variables={dim: new_var})
        
    def cf_var(self, variable: str = None) -> None:
        # Check variable exists in self.ds.

        # Find variable in vocabulary

        # Apply changes to variable where necessary

        # Check units are valid and apply standard formatting (e.g. kg/kg becomes 1)
        
        # return new version of variable's data array.
        pass


# class MoncVariable(xr.DataArray):
#     __slots__ = ()

#     def __init__(self, data: Any = dtypes.NA, coords: Sequence[Tuple] | Mapping[Any, Any] | None = None, dims: Hashable | Sequence[Hashable] | None = None, name: Hashable = None, attrs: Mapping = None, indexes: Dict = None, fastpath: bool = False, data_array: xr.DataArray = None):
        
#         if data_array is not None:
#             super().__init__(data=data_array.data, coords=data_array.coords, dims=data_array.dims, name=data_array.name, attrs=data_array.attrs, indexes=data_array.indexes)
#         else:
#             super().__init__(data=data, coords=coords, dims=dims, name=name, attrs=attrs, indexes=indexes, fastpath=fastpath)


def test():
    # data_dir = os.path.join(os.path.dirname(app_dir), 'test_data')
    
    # print("Testing MoncDataset instantiation: passing in filepath")
    # filename = 'd20200128_diagnostic_1d_3600.nc'
    # print(f"Is {filename} MONC output? {MoncDataset(filepath=os.path.join(data_dir, filename)).is_monc()}")

    # print("Testing MoncDataset instantiation: passing in xarray.Dataset")
    # filename = 'd20200128_diagnostic_1d_3600.nc'
    # ds = xr.open_dataset(os.path.join(data_dir, filename), decode_times=False)
    # print(f"Is {filename} MONC output? {MoncDataset(dataset=ds).is_monc()}")

    


if __name__ == '__main__': test()