from difflib import SequenceMatcher
import os.path as op
import os
import xarray as xr
from glob import iglob
from setup import *
import re


def stem_str(*args: str):
    '''
    Finds common segment of a list of filepaths/names. Also works on any 
    strings.
    '''
    if not args:
        raise ValueError('stem_str: At least one string must be provided.')
    if None in args:
        return stem_str(*[arg for arg in args if arg is not None])
    
    filenames = [op.splitext(op.basename(path))[0] for path in args]
    # If each path contains only a filename (no extension), or a string other 
    # than a filepath, this will still work.
    
    stem = filenames[0]
    if len(args) > 1:
        for f in filenames[1:]:
            match = SequenceMatcher(a=stem, b=f).find_longest_match()
            stem = stem[match.a: match.a + match.size]
    
    return stem


def get_time_var(*args: int|xr.Dataset) -> str:
    '''
    
    '''
    # If argument is an integer, assume it's the number of spatial dimensions
    if isinstance(args[0], int):
        if args[0] not in vocabulary:
            raise KeyError(
                'get_time_var: Vocabulary file has no specifications for '
                f'{args[0]} spatial dimensions.')
        time_vars = [v for v in vocabulary[args[0]].keys() 
                     if 'time' in v.lower()]
        if not time_vars:
            raise KeyError('get_time_var: No time variable found in '
                           f'vocabulary for {args[0]}d datasets.')
    elif isinstance(args[0], xr.Dataset):
        time_vars = [d for d in args[0].dims if 'time' in d.lower()]
        if not time_vars:
            raise KeyError(f'No time variable found in dataset.')
    else:
        raise AttributeError(
            'get_time_var: Either a number of spatial dimensions or an xarray.'
            'Dataset must be passed in as a parameter.')
    if len(time_vars) > 1:
        as_string = ''
        for i, s in enumerate(time_vars):
            as_string += f'\n[{i}] {s}'
        while True:
            sel = input(f'Please select time variable for {args[0]}d datasets. Enter number: {as_string}')
            if isinstance(sel, int) and sel >= 0 and sel < len(time_vars):
                break
            else:
                print('Invalid selection.')
        return time_vars[sel]
    return time_vars[0]


def ds_feed(file_list: list|set|tuple|str):
    '''
    This generator replaces the DsGroup.datasets generator, because generators 
    can't be passed (as a "pickle") to a new Process as an attribute, method or 
    function within an object.
    '''
    # ENHANCEMENT: test that all paths supplied are valid & are NC files.
    if isinstance(file_list, str):
        file_list = [file_list]
    return (xr.open_dataset(path, decode_times=False) for path in file_list)
    

def variable_glob(ds: xr.Dataset, var_glob: str) -> list:
    '''
    Takes a glob-style string and finds matching variable(s) in a dataset.
    '''
    match_str = re.compile(var_glob.replace('*', '.*').replace('?','.'))
    return [match_str.fullmatch(v).string for v in set(ds.variables) if match_str.fullmatch(v)]


class DsGroup:
    '''
    Each file in a DsGroup has the same variables, coordinates & dimensions.
    '''
    def __init__(self, 
                 name: str = '', 
                 n_dims: int = None,
                 time_variable: str = '',
                 action: str = None, 
                 filepaths: list|set|tuple = None,
                 groups: list|set|tuple = None) -> None:
        '''
        n_dims: Number of *spatial* dimensions (X, Y, Z), ignoring variations 
                such as z & zn.
        '''
        # If groups parameter contains existing DsGroup objects:
        if groups and all([isinstance(g, DsGroup) for g in groups]):
            # check each DsGroup has only one dataset in its processed 
            # attribute.
            if not all([isinstance(g.processed, xr.Dataset) for g in groups]):
                raise TypeError('DsGroup: To merge DsGroup objects, the '
                                'processed attribute of each must contain a '
                                'single xarray.Dataset object.')
            self.name = stem_str(*[group.name for group in groups])
            self.n_dims = sum([g.n_dims for g in groups])  # This is a little dubious, given there's no guarantee dimension groups don't overlap, but works for MONC outputs where 0d & 1d are grouped.
            self.time_var = time_variable or groups[0].time_var
            if not all([g.time_var == self.time_var for g in groups]):
                raise ValueError('DsGroup: When attempting to combine groups, '
                                 'found mismatch in time variables.')
            self.action = action or 'merge_groups'
            if 'merge' in self.action:
                self.process = self.merge_groups
            self.stem = stem_str(*[g.stem for g in groups])
            self.filepaths = None
            self.to_process = [g.processed for g in groups]  # Dataset objects
            self.processed = None  # Dataset object, when completed
        else:
            self.n_dims = n_dims    # TODO: Could set to parse files to infer 
                                    # if not given.
            self.name = name if name else f'{str(n_dims)}d'
            # If no action is specified, will only be processed for CF 
            # compliance, but files won't be merged/split.
            self.action = action.lower() if action else None
            self.filepaths = filepaths if filepaths else []
            # Cannot set default valueas an empty list in arguments, or it 
            # assigns SAME empty list to every instance.

            self.processed = None  # Filepath(s) of output
            self.stem = stem_str(*self.filepaths) if self.filepaths else None

            # Although it would be best to verify that all member datasets have 
            # the same time variable, for now, just use the first one.
            try:
                self.time_var = time_variable or get_time_var(self.n_dims)
            except KeyError:
                if self.filepaths:
                    # Although it would be best to verify that all member 
                    # datasets have the same time variable, for now, just use 
                    # the first one.
                    try:
                        self.time_var = get_time_var(
                            xr.open_dataset(
                            self.filepaths[0], decode_times=False
                            ))
                    except:
                        self.time_var = None    # Search in datasets as files 
                                                # are added.

                    # TODO: verify all datasets in group have same time 
                    # variable & time points.
                    # if all([tv==time_var[0] for tv in time_var.values()]):
                    #     time_var = time_var[0]
                    # else:
                    #     #TODO: disambiguate.
                    #     raise ValueError(
                    #         f"Inconsistent time variables across dataset group {group.name}."
                    #         )
                else: 
                    self.time_var = None    # Search in datasets as files are 
                                            # added.

    def datasets(self):
        return ds_feed(self.filepaths)
    
    def add(self, filepath: str) -> None:
        
        if not op.exists(filepath):
            raise OSError(f'DsGroup.add: Filepath {filepath} not found.')
        
        if op.isdir(filepath):
            for path in iglob(op.join(filepath, '*.nc'), 
                              os.pathsep, 
                              recursive=False):
                self.add(path)

        # Add filepath to collection
        self.filepaths.append(filepath)
        
        # Update name stem common to all files in group.
        self.stem = stem_str(self.stem, filepath)

        # print(f"Added {op.basename(filepath)} to {self.name}; stem: {self.stem}")

        # If don't have a time variable yet, attempt to find in new file.
        if not self.time_var:
            try:
                self.time_var = get_time_var(
                    xr.open_dataset(
                        filepath, 
                        decode_times=False
                    )
                )
            except KeyError:
                print(f'DsGroup.add: Time variable not found in {filepath}.')
        # If time variable contains wildcards, seek matching variable in ds. If
        # found, update vocabulary accordingly.
        elif len({'*', '?'}.intersection(self.time_var)) > 0:
            exact_time_var = variable_glob(
                xr.open_dataset(
                    filepath, 
                    decode_times=False
                ), self.time_var
            )
            if len(exact_time_var) != 1:
                # TODO: disambiguation
                raise AttributeError(
                    'Multiple variables found in dataset that match time '
                    f'variable pattern, {self.time_var}.')
            exact_time_var = exact_time_var[0]
            vocabulary[self.n_dims][exact_time_var] = vocabulary[self.n_dims][self.time_var]
            vocabulary[self.n_dims].pop(self.time_var)
            self.time_var = exact_time_var
        else:
            # TODO: confirm time variable in new dataset matches self.time_var
            pass

    def merge_times(self) -> None:
        pass

    def merge_groups(self) -> None:
        pass


