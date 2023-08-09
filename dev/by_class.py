import xarray as xr
from setup import *
import os
from glob import glob
from cfunits import Units
from dataset_functions import *
from variable_functions import *
from difflib import SequenceMatcher


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
        self.processed = []

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
                    n_dims = get_n_dims(ds) - 1
                    self.by_dim[n_dims].add(filepath=filepath)
                else:
                    self.input_files.append(filepath)

    def merge_groups(self):
        to_merge = {}
        for name, action in self.DIM_ACTIONS.items():
            if action == 'merge':
                merge_group = [group for group in self.by_dim.values() if group.name == name]
                if len(merge_group) > 1:
                    to_merge[name] = merge_group
        merged = {}
        for name, merge_list in to_merge.items():
            merged[name] = merge_dims(merge_list)

    def datasets(self):
        # Generator
        return (xr.open_dataset(path, decode_times=False) for path in self.monc_files)
    

class DsGroup:
    '''
    Each file in a DsGroup has the same variables, coordinates & dimensions.
    '''
    def __init__(self, 
                 name: str, 
                 action: str = None, 
                 filepaths: list|set|tuple = []) -> None:
        
        # TODO: validate parameters

        self.name = name
        if action.lower() == 'split':
            self.action = self.split_times
        elif action.lower() == 'merge':
            self.action = self.merge_times
        else:
            self.action = None
        self.filepaths = filepaths
        self.processed = []
        self.stem = self.name_stem(filepaths)

    def datasets(self):
        '''
        This returns a generator that yields the datasets contained in each
        filepath.
        '''
        return (xr.open_dataset(path, decode_times=False) for path in self.filepaths)
    
    def add(self, filepath):
        # TODO: validate filepath
        self.filepaths.append(filepath)
        self.stem = self.name_stem(filepath, self.stem)

    def name_stem(self, *args):
        if not args:
            raise ValueError('At least one filename/path must be provided.')
        filenames = [os.path.splitext(os.path.basename(path))[0] for path in args]  # If each path contains only a filename (no extension) this will still work.
        self.stem = filenames[0]
        if len(args) > 1:
            for f in filenames[1:]:
                match = SequenceMatcher(a=self.stem, b=f).find_longest_match()
                self.stem = self.stem[match.a: match.a + match.size]

    def process(self):
        if self.action is not None:
            self.action()

    def merge_times(self):
        merged = None
        self.processed = [merged]

    def split_times(self):
        split = None
        self.processed += split


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
    

# class MoncDataset(xr.Dataset):
#     __slots__ = ()  # __slots__ must be declared in subclasses of any xarray classes.
    
#     def __init__(self, filepath: str = None, dataset: xr.Dataset = None):
#         # This is quite a messy way to set up, but does allow a filepath to be
#         # passed in for dataset instantiation, rather than having to open it
#         # as an xr.Dataset first. Calling super().__init__ is a convenient way
#         # of setting up all inherited attributes & methods.

#         if dataset is not None:
#             ds = dataset
#         elif filepath is not None:
#             ds = xr.open_dataset(filepath, decode_times=False)
#         else:
#             raise ValueError('Either a filepath (string) or xarray.Dataset must be passed in.')
#         data_vars = [MoncVariable(data_array=var) for var in ds.data_vars.values()]
#         super().__init__(ds.data_vars, ds.coords, ds.attrs)
        

#     def is_monc(self) -> bool:
#         return MONC_ID_ATTR in self.attrs

#     def add_cf_attrs(self):
#         pass

#     def attr_to_var(self, attr: str = None):
#         if attr is None:
#             # convert all attributes to variables, using this function recursively.
#             for attr in self.attrs.keys():
#                 self.attr_to_var(attr=attr)
#         # TODO

#     def missing_coords(self):
#         pass


# class MoncVariable(xr.DataArray):
#     __slots__ = ()

#     def __init__(self, data: Any = dtypes.NA, coords: Sequence[Tuple] | Mapping[Any, Any] | None = None, dims: Hashable | Sequence[Hashable] | None = None, name: Hashable = None, attrs: Mapping = None, indexes: Dict = None, fastpath: bool = False, data_array: xr.DataArray = None):
        
#         if data_array is not None:
#             super().__init__(data=data_array.data, coords=data_array.coords, dims=data_array.dims, name=data_array.name, attrs=data_array.attrs, indexes=data_array.indexes)
#         else:
#             super().__init__(data=data, coords=coords, dims=dims, name=name, attrs=attrs, indexes=indexes, fastpath=fastpath)


def test():
    data_dir = os.path.join(os.path.dirname(app_dir), 'test_data')
    
    # print("Testing MoncDataset instantiation: passing in filepath")
    # filename = 'd20200128_diagnostic_1d_3600.nc'
    # print(f"Is {filename} MONC output? {MoncDataset(filepath=os.path.join(data_dir, filename)).is_monc()}")

    # print("Testing MoncDataset instantiation: passing in xarray.Dataset")
    # filename = 'd20200128_diagnostic_1d_3600.nc'
    # ds = xr.open_dataset(os.path.join(data_dir, filename), decode_times=False)
    # print(f"Is {filename} MONC output? {MoncDataset(dataset=ds).is_monc()}")


if __name__ == '__main__': test()