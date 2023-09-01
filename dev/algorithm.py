#!/usr/bin/env python

# Initialise, incl. from VOCAB & CONFIG files.
from setup import *
import xarray as xr
import os
import os.path as op
from glob import iglob
from cfunits import Units
import numpy as np
from datetime import datetime, timezone
import typing
# from utils import type_from_str, tf, performance_time
# from dask.diagnostics import ProgressBar
from time import perf_counter, localtime, strftime
import argparse
from multiprocessing import Process, Pool #, set_start_method  # multiprocess not available in Jaspy environment.
from cfize_ds import *
from units import TimeUnits
from groups import DsGroup, globals_to_vars
import re
import sys


def process_large(
        filepath: str, 
        group: DsGroup,
        title: str, 
        shared: dict
    ) -> str:

    if shared['verbose']:
        start_time = perf_counter()

    warnings = []
    errors = []
    log = []
    
    update_globals = {}  # This is not currently used, but is here for consistency with cf_merge

    # time_units = shared['time_units']
    target_dir = shared['target_dir']
    # reference_vars = shared['reference_vars']

    # if global_vars:
    #     global reference_vars, vocabulary
    #     reference_vars = global_vars['reference_vars']
    #     vocabulary = global_vars['vocabulary']

    group.processed = []

    if shared['verbose']:
        log.append(f"{strftime('%H:%M:%S', localtime())} Process {os.getpid()}: process_large running on group "
                   f"{group.name} - file {op.basename(filepath)}. Title passed "
                   f"in: {title}.")
        # print(
        #     f"Process {os.getpid()}: process_large running on group "
        #     f"{group.name} - file {op.basename(filepath)}."
        #     f"Title passed in: {title}."
        # )
    
    with xr.open_dataset(filepath, 
                         decode_times=False
                         ) as dataset:
        
        # apply chunking if large (e.g. ds.nbytes >= 1e9).
        if dataset.nbytes >= 1e9:
            dataset.chunk(chunks=CHUNKING_DIMS)

        # Assign any required global attributes to variables, 
        # with associated attributes. Remove global attribute.
        '''
        Assume these are only correct for last time-point in series.
        '''
        
        last_time_point = dataset[group.time_var].data[-1]  # time_var[i]
        try:
            new_vars = globals_to_vars(
                        ds=dataset, 
                        time_var=group.time_var, 
                        last_time_point=last_time_point
            )  # time_var[i]
        except Exception as e:
            errors.append['process_large: globals_to_vars: ' + str(e)]
            return {
                'warnings': warnings,
                'errors': errors,
                'log': log
            }
        for name, array in new_vars.items():
            dataset.attrs.pop(name)
            dataset = dataset.assign({name: array})  # Overwrite with updated

        # Create MoncDs object for futher processing, and update dataset's title
        try:
            monc_ds = MoncDs(
                dataset=dataset,
                time_variable=group.time_var,
                n_dims=group.n_dims,
                title=title
            )
        except AttributeError as e:
            errors.append('process_large: MoncDs: ' + str(e))
            return {
                'warnings': warnings,
                'errors': errors,
                'log': log
            }
        #Call CF compliance function:
            # Adds any missing global attributes required by CF convention
            # Derives any required global attributes from options database
        try:
            monc_ds.cfize(shared=shared)
        except (ConfigError or VocabError or AttributeError) as e:
            log += monc_ds.log
            errors.append('process_large: MoncDs.cfize: ' + str(e))
            warnings += monc_ds.warnings
            return {
                'warnings': warnings,
                'errors': errors,
                'log': log
            }
        log += monc_ds.log
        warnings += monc_ds.warnings
        

        if group.action == 'split':
            # Split dataset by time-point, yielding multiple 
            # new datasets. Append time to titles.
            (processed, split_log) = split_ds(
                dataset=monc_ds.ds,
                shared=shared, 
                var=monc_ds.time_var
            )
            log += split_log
            log.append(
                f"split_ds returned datasets with titles " +
                ", ".join(ds.attrs['title'] for ds in processed)
            )
            # processed only ever needs to hold latest collection of datasets.
        else:
            processed = monc_ds.ds  # Will this persist after context of ds ends?

    writers = []
    # For each new dataset:
    for i, ds in enumerate(processed):
        # Update required global attributes (MONC time, previous diagnostic write at).
        for attr, (func, attr_type) in split_attrs.items():
            if attr in ds.attrs:
                ds.attrs[attr] = attr_type(func(i, processed))
            else:
                errors.append(AttributeError(
                    f"setup.py: split_attrs specifies derivation of {attr} "
                    f"in split datasets, but attribute was not found."
                ))
                return {
                    'warnings': warnings,
                    'errors': errors,
                    'log': log
                }

        # set filepath
        filepath = op.join(
            target_dir, f"{ds.attrs['title']}.nc"
            )
        # set encodings
        encodings = {
            k:{
                'dtype': v.dtype,
                '_FillValue': None,
                COMPRESSION[1]: COMPRESSION[0],
                'complevel': COMPRESSION[2]
            } for k, v in ds.variables.items()
        } if COMPRESSION[0] else {
            k:{
                'dtype': v.dtype,
                '_FillValue': None
            } for k, v in ds.variables.items()
        }

        # Export to NetCDF (set as single-command function, 
        # so can wrap in performance_time).
        # Testing dask version
        if shared['verbose']: 
            log.append(
                f"{strftime('%H:%M:%S', localtime())} Process {os.getpid()}: "
                f"preparing delayed writer for {filepath}."
            )
            # print(
            #     f"Process {os.getpid()}: preparing delayed writer for "
            #     f"{filepath}."
            # )
        try:
            if shared['verbose']: w_start = perf_counter()
            writers.append(ds_to_nc_dask(ds=ds,
                                    filepath=filepath,
                                    encodings=encodings,
                                    compress=COMPRESSION[0],
                                    shared=shared
                                    ))
            if shared['verbose']: log.append(
                f"         Process {os.getpid()}: ds_to_nc_dask took "
                f"{perf_counter() - w_start} seconds."
            )
            ds.close()
            group.processed.append(filepath)
        except Exception as e:
            errors.append('process_large: ds_to_nc_dask: ' + str(e))
            return {
                'warnings': warnings,
                'errors': errors,
                'log': log
            }

    if shared['verbose']: 
        log.append(
            f"{strftime('%H:%M:%S', localtime())} Process {os.getpid()}: "
            f"Computing writes."
            )
        # print(f"Process {os.getpid()}: Computing writes.")
    try:
        [perform_write(writer, shared=shared) for writer in writers]
        for writer in writers:
            if shared['verbose']: w_start = perf_counter()
            perform_write(writer, shared=shared)
            if shared['verbose']: 
                log.append(
                    f"         Process {os.getpid()}: perform_write took "
                    f"{perf_counter() - w_start} seconds."
                )
    except Exception as e:
        errors.append('process_large: perform_write: ' + str(e))
        return {
            'warnings': warnings,
            'errors': errors,
            'log': log
        }
    
    # TODO: Run CF checker on output files
    
    if shared['verbose']:
        log.append(
            f"         Process {os.getpid()}: process_large took "
            f"{perf_counter() - start_time} seconds."
        )
        # print(
        #     f"Process {os.getpid()}: process_large took "
        #     f"{perf_counter() - start_time} seconds."
        # )

    return {
        'update_group':{
            'processed': group.processed, 
            'time_var': group.time_var
        },
        'update_globals': update_globals,
        'warnings': warnings,
        'errors': errors,
        'log': log
    }
    # return (group, shared)


def cf_merge(group: DsGroup, 
             shared: dict) -> str:
    
    errors = []
    warnings = []
    log = []
    if shared['verbose']: 
        log.append(
            f"{strftime('%H:%M:%S', localtime())} Process {os.getpid()}: "
            f"cf_merge running on group {group.name}"
        )
        # print(
        #     f"Process {os.getpid()}: cf_merge running on group {group.name}"
        # )
    
    try:
        rtn = group.merge_times(shared=shared)
    except Exception as e:
        errors.append('cf_merge: DsGroup.merge_times: ' + str(e))
        return {
            'warnings': warnings,
            'errors': errors,
            'log': log
        }
    else:
        if 'error' in rtn and rtn['error']:
            errors.append('cf_merge: DsGroup.merge_times: ' + str(rtn['error']))
            return {
                'warnings': warnings,
                'errors': errors,
                'log': log
            }
        if 'log' in rtn and rtn['log']:
            log += rtn['log']
    try:
        rtn = group.cfize_and_save(shared=shared)
    except Exception as e:
        errors.append('cf_merge: DsGroup.cfize_and_save: ' + str(e))
        return {
            'warnings': warnings,
            'errors': errors,
            'log': log
        }
    else:
        update_globals = rtn['update_globals'] if 'update_globals' in rtn else None
        update_group = rtn['update_group'] if 'update_group' in rtn else {
            'processed': group.processed, 
            'time_var': group.time_var
        }
        if 'errors' in rtn and rtn['errors']:
            errors.append(
                'cf_merge: DsGroup.cfize_and_save: ' + '; '.join(
                    [str(e) for e in rtn['errors']]
                )
            )
        warnings += rtn['warnings'] if 'warnings' in rtn else []
        if 'log' in rtn and rtn['log']:
            log += rtn['log']
    
    return {
        'update_group': update_group,
        'update_globals': update_globals,
        'warnings': warnings,
        'errors': errors,
        'log': log
    }


def process_parallel(
        groups: dict, 
        n_proc: int, 
        shared: dict
):
    '''
    groups:     dimension-specific groups to be processed.
    n_proc:     number of processes in process pool (doesn't include controller)
    shared:     dictionary containing what would be global variables if all
                functions were running on a single process. Global constants
                will still work as normal.
    '''

    target_dir = shared['target_dir']
    time_units = shared['time_units']
    # vocabulary = shared['vocabulary']
    reference_vars = shared['reference_vars']
    warnings = []

    # global reference_vars, vocabulary
    # Verify n_proc + controller doesn't exceed available cores
    # TODO: confirm this works on Jasmin!
    if n_proc + 1 > os.cpu_count(): n_proc = os.cpu_count() - 1

    results = {dim: [] for dim in groups.keys()}  # set up empty dict for results from process pool

    # Set up process pool
    # set_start_method('fork')  # only available on Unix systems. not Mac/Win.
    # fork method allows each new process to inherit the parent's whole state
    # (global variables etc). Otherwise, each spawned process will have only
    # values set at initialisation.
    with Pool(processes=n_proc) as pool:
        
        # Allocate datasets to be processed in increasing order of dimensions
        # Work with 0-2d files first, because 3d processing depends on 0d/1d.
        '''
        ENHANCEMENT: this processing in order of increasing dimension is 
        critical for converting perturbations to absolutes. For non-MONC data, 
        it may not be safe to assume that reference variables will be in
        lower-dimension datasets, so waiting for their completion would hang 
        execution or fail. A safer alternative would be to export any reference 
        variables to a standalone NetCDF file, which can then be accessed when 
        needed. In this case, there should be a flag to indicate when the
        required variable is available.
        '''
        for dim in sorted(list(groups.keys())):
            group = groups[dim]
            if shared['verbose']: print(
                    f"{strftime('%H:%M:%S', localtime())} "
                    f"Process {os.getpid()}: Setting up parallel processing of "
                    f"{dim}:{group.name}."
                )
            
            if group.action == 'split' or dim == 3:
                # The 2nd condition ensures that large files are allocated to
                # the process pool, whether they are to be split or not.
                
                # Derive base title/filename from group's stem & dimension
                title = group.stem + group.name if group.name not in group.stem else group.stem
                
                for_controller = []

                # Wait for 0+1d processing to finish, so reference variables
                # are available for perturbations.
                for d in {v['dim'] for v in reference_vars.values()}:
                    for r in results[d]:
                        result_dict = r.get()
                        if 'log' in result_dict and result_dict['log']:
                            print(''.join(
                                [entry + '\n' for entry in result_dict['log']]
                            ))
                        if 'errors' in result_dict and result_dict['errors']:
                            err_msg = '; '.join([str(e) for e in result_dict['errors']])
                            sys.exit(err_msg)
                        if 'warnings' in result_dict:
                            warnings += result_dict['warnings']
                        if 'update_group' in result_dict:
                            update_group = result_dict['update_group']
                            if groups[d].processed:
                                groups[d].processed = list({
                                    *update_group['processed'], 
                                    *groups[d].processed
                                })
                            else:
                                groups[d].processed = update_group['processed']
                            groups[d].time_var = update_group['time_var']
                        else:
                            update_group = None
                        update_globals = result_dict['update_globals'] if 'update_globals' in result_dict else None
                        if update_globals:
                            if 'vocabulary' in update_globals: 
                                [[shared['vocabulary'][v_dim][update_var].update(update_dict) 
                                  for update_var, update_dict in updates.items()] 
                                  for v_dim, updates in update_globals['vocabulary'].items()]
                                # for v_dim, updates in update_globals['vocabulary'].items():
                                #     for update_var, update_dict in updates.items():
                                #         shared['vocabulary'][v_dim][update_var].update(update_dict)
                        # [reference_vars.update({v: result[v]}) 
                        #  for v in result.keys() if 'filepath' in result[v]]
                    # groups[d].processed = [
                    #     result.get() for result in results[d]
                    # ]  # This should only set groups[d].processed to its
                    #     # current value, but serves as a "completed" flag.
                
                # Round-robin allocation of file processing to controller and 
                # workers:
                for i, filepath in enumerate(group.filepaths):
                    # Allocate one file per round to controller
                    if i%(n_proc+1) == 0:
                        for_controller.append(filepath)
                    else:
                        results[dim].append(
                            pool.apply_async(
                                func=process_large, 
                                kwds={
                                    'filepath':filepath, 
                                    'group':group,
                                    'title':title, 
                                    'shared':shared
                                }
                            )
                        )
                
                # Now work on remaining files on controller process
                # group.processed = []
                # [
                #     [
                #         group.processed.append(ds) 
                #         for ds in process_large(
                #             filepath=filepath,
                #             group=group,
                #             title=title,
                #             shared=shared
                #         )
                #     ] for filepath in for_controller
                # ]
                for filepath in for_controller:
                    r = process_large(
                        filepath=filepath,
                        group=group,
                        title=title,
                        shared=shared
                    )
                    if 'log' in r and r['log']:
                        print(''.join(
                            [entry + '\n' for entry in r['log']]
                        ))
                    if 'errors' in r and r['errors']:
                        sys.exit('; '.join([str(e) for e in r['errors']]))
                    if 'warnings' in r:
                        warnings += r['warnings']
                    update_group = r['update_group'] if 'update_group' in r else None
                    update_globals = r['update_globals'] if 'update_globals' in r else None
                    if groups[dim].processed:
                        groups[dim].processed = list({
                            *groups[dim].processed, 
                            *update_group['processed']
                        })
                        # process_large(
                            #     filepath=filepath,
                            #     group=group,
                            #     title=title,
                            #     shared=shared
                            # )[0]
                    else:
                        groups[dim].processed = update_group['processed']
                        # process_large(
                        #         filepath=filepath,
                        #         group=group,
                        #         title=title,
                        #         shared=shared
                        #     )[0]
                    groups[dim].time_var = update_group['time_var'] 
                
                # Gather completed jobs
                # TODO: is this efficient, or can it be done after merging 
                # groups?
                # If this moves below, each group's processed attribute should
                # be updated in the process_large function, then the group
                # returned.
                # [result.get() for result in results[dim]]
                for result in results[dim]:
                    # [update_group, 
                    #  update_globals] = result.get()
                    r = result.get()
                    if 'log' in r and r['log']:
                        print(''.join(
                            [entry + '\n' for entry in r['log']]
                        ))
                    if 'errors' in r and r['errors']:
                        sys.exit('; '.join([str(e) for e in r['errors']]))
                    if 'warnings' in r:
                        warnings += r['warnings']
                    if 'update_group' in r:
                        update_group = r['update_group']
                        if groups[dim].processed:
                            groups[dim].processed = list(
                                {*groups[dim].processed,
                                *update_group['processed']}
                            )
                        else:
                            groups[dim].processed = update_group['processed']
                        groups[dim].time_var = update_group['time_var']
                    else:
                        update_group = None
                    update_globals = r['update_globals'] if 'update_globals' in r else None
                    
            # If group is to be merged:
            elif group.action == 'merge':
                # results[dim] = pool.apply_async(
                #     func=cf_merge,
                #     kwds={
                #         'group': group,
                #         'time_units': time_units,
                #         'target_dir': target_dir
                #     }
                # )
                results[dim].append(
                    pool.apply_async(
                        func=cf_merge,
                        kwds={
                            'group': group,
                            'shared': shared
                        }
                    )
                )
        
            # Otherwise, process for CF compliance, but leave dataset as a 
            # single, standalone file.
            else:
                # ENHANCEMENT: Could still evaluate by size whether it's 
                # processed locally or passed to the Pool.

                # Round-robin allocation of file processing to controller and 
                # workers:
                results[dim].append(
                    pool.apply_async(
                        func=cfize_all_datasets,
                        kwds={
                            'group': group,
                            'shared': shared
                        }
                    )
                )

        # Check all other processes are complete
        # for result in results:
        #     print(result.get())  # assuming here that and DirectoryParser.process_by_dim and process_large each return values.
        for dim in sorted(list(groups.keys())):
            # Sort because expect smaller to finish first
            if not groups[dim].processed:
                groups[dim].processed = []
                for result in results[dim]:
                    # [update_group, 
                    #  update_globals] = result.get()
                    r = result.get()
                    if 'log' in r and r['log']:
                        print(''.join(
                            [entry + '\n' for entry in r['log']]
                        ))
                    if 'errors' in r and r['errors']:
                        sys.exit('; '.join([str(e) for e in r['errors']]))
                    if 'warnings' in r:
                        warnings += r['warnings']
                    if 'update_group' in r:
                        update_group = r['update_group']
                        groups[dim].processed = list(
                            {*groups[dim].processed,
                             *update_group['processed']}
                        )
                        groups[dim].time_var = update_group['time_var']
                    else:
                        update_group = None
                    update_globals = r['update_globals'] if 'update_globals' in r else None
                    
                # [groups.update({dim: result.get()[0]}) for result in results[dim]]
                # groups[dim].processed = [
                #     result.get() for result in results[dim]
                # ]
    return warnings


def cfize_all_datasets(
        group: DsGroup, 
        shared: dict
    ) -> list:
    '''
    Should print & return names of saved files
    TODO: This function still to be written tested.
    '''
    
    update_globals = {}
    errors = []
    warnings = []
#     print(f"Processing {group.n_dims}:{group.name}.")

#     files = []

#     processing = list(group.datasets())  # Open all datasets in group.
#     for i, dataset in enumerate(processing):
#         with MoncDs(
#             dataset=dataset,
#             time_units=time_units,
#             time_variable=group.time_var,
#             n_dims=group.n_dims,
#             title=group.stem.strip('_ ') + group.name if group.name not in group.stem else group.stem.strip('_ ')
#         ) as ds:
    
#             # Call CF compliance function on dataset.
#             processed = ds.cfize()  # xarray.Dataset

#             # Append title with (last) time-point
#             title = processed.attrs['title'].strip(' _') + '_' + processed[ds.time_var].data.tolist()[-1]

#             # Set filepath and save as processed
#             filepath = op.join(target_dir, f"{title}.nc")

#             # TODO: can't yet provide filepath to reference variables if not merged;
#             # this would need references_vars[v]['filepath'] to contain a list or
#             # dict with e.g. timepoints as keys.
#             # # Check for any reference variables for perturbation variables, and note
#             # # filepath if found:
#             # [
#             #     reference_vars[v].update({'filepath':filepath})
#             #     for v in reference_vars.keys()
#             #     if v in processed.variables
#             # ]

#             # Set encoding
#             # TODO: encoding needs to be set, with each variable's encoding specified. 
#             # Otherwise, _FillValue is applied to all, including coordinates, the 
#             # latter contravening CF Conventions.
#             encodings = {
#                 k:{
#                     'dtype': v.dtype,
#                     '_FillValue': None
#                 } for k, v in processed.variables.items()
#             }  # if k == 'options_database' or k in ds.ds.coords

#             # Save dataset to NetCDF
#             # xarray docs report engine='h5netcdf' may sometimes be 
#             # faster. However, it doesn't natively handle the 
#             # various string length types used here.
#             processed.to_netcdf(
#                 path=filepath, 
#                 encoding=encodings)  # , engine='h5netcdf'

#             # TODO: Run each file through cf-checker?
            
#             # Close dataset
#             processed.close()

#         files.append(filepath)

    return {
        'update_group':{
            'processed': group.processed, 
            'time_var': group.time_var
        },
        'update_globals': update_globals,
        'warnings': warnings,
        'errors': errors
    }


def process_serial(
        groups: dict, 
        shared: dict
):
    '''
    groups:         dimension-specific groups to be processed.
    shared:         dictionary containing what would be global variables if all
                    functions were running on a single process. Global constants
                    will still work as normal.
        target_dir: directory in which to write processed NC files.
        time_units: CF-compliant unit to use for time coordinate.
    '''

    target_dir = shared['target_dir']
    time_units = shared['time_units']
    warnings = []

    # Processe datasets in increasing order of dimensions.
    # Work with 0-2d files first, because 3d processing depends on 0d/1d.
    '''
    ENHANCEMENT: this processing in order of increasing dimension is 
    critical for converting perturbations to absolutes. For non-MONC data, 
    it may not be safe to assume that reference variables will be in
    lower-dimension datasets, so waiting for their completion would hang 
    execution or fail. A safer alternative would be to export any reference 
    variables to a standalone NetCDF file, which can then be accessed when 
    needed. In this case, there should be a flag to indicate when the
    required variable is available.
    '''
    for dim in sorted(list(groups.keys())):
        group = groups[dim]
        if shared['verbose']: print(f"{strftime('%H:%M:%S', localtime())}"
                                    f"Serial processing {dim}:{group.name}.")
        
        if group.action == 'split':
            
            # Derive base title/filename from group's stem & dimension
            title = group.stem + group.name if group.name not in group.stem else group.stem
            
            group.processed = [process_large(
                    filepath=filepath,
                    group=group,
                    title=title,
                    shared=shared
                ) for filepath in group.filepaths]

        # If group is to be merged:
        elif group.action == 'merge':
            try:
                rtn = group.merge_times(shared=shared)
                if 'log' in rtn:
                    log += rtn['log']
                if 'error' in rtn and rtn['error']:
                    sys.exit(log + '\nprocess_serial: DsGroup.merge_times' + str(rtn['error']))

            # group.processed = cf_merge(
            #         group=group,
            #         time_units=time_units,
            #         target_dir=target_dir
            # )
                rtn = group.cfize_and_save(shared=shared)
                if 'errors' in rtn and rtn['errors']:
                    sys.exit('process_serial: DsGroup.cfize_and_save: ' + 
                             '; '.join([str(e) for e in rtn['errors']]))
                if 'warnings' in rtn:
                    warnings += rtn['warnings']
            except Exception as e:
                sys.exit('process_serial (merging): ' + str(e))
    
        # Otherwise, process for CF compliance, but leave dataset as a 
        # single, standalone file.
        else:
            # TODO: to be tested
            group.cf_only(time_units=time_units, target_dir=target_dir)
            # group.processed = cf_only(
            #     group=group,
            #     time_units=time_units,
            #     target_dir=target_dir
            # )
    return warnings


def sort_nc(directory, shared: dict) -> [dict, list]:
    '''
    Sorts into MONC output files and other NC files, the latter listed as
    possible input files.

    It also categorizes the MONC outputs according to their number of 
    dimensions.

    It would be best if all files of a given run were in one directory, 
    rather than split between 0-2d and 3d.
    '''
    if verbose: print(strftime('%H:%M:%S', localtime()), 'Categorising files by dimension')

    # global vocabulary, reference_vars

    # by_dim = {n_dim: DsGroup(
    #     name=group, n_dims=n_dim, action=DIM_ACTIONS[group], shared=shared
    # ) for n_dim, group in DIM_GROUPS.items()}
    by_dim = {}
    for n_dim, group in DIM_GROUPS.items():
        by_dim[n_dim] = DsGroup(
            name=group, n_dims=n_dim, action=DIM_ACTIONS[group], shared=shared
        )
        if by_dim[n_dim].log:
            print(
                [entry + '\n' for entry in by_dim[n_dim].log]
            )
        
    input_files = []

    # For each NC file in source directory (recursive parsing or not?):
    for filepath in iglob(f'{op.join(directory, "*.nc")}', 
                          root_dir=os.pathsep, 
                          recursive=False):
        with xr.open_dataset(filepath, 
                             decode_times=False) as ds:  
                            # concat_characters=False
            # Using concat_characters=False preserves the string dimension 
            # of the options_database variable. However, it makes dealing
            # with the binary -> string conversion painful and messy.

            # Categorise MONC output / other
            if is_monc(ds):
                # Find number of (non-options-database) dimensions.
                n_dims = get_n_dims(ds) - 1  # Subtract 1 for time.
                # Add filepath to relevant dimension group.
                by_dim[n_dims].add(filepath=filepath)
                # print(f"{op.basename(filepath)}: {n_dims} spatial dimensions.")
            else:
                # Categorise as potential input file
                input_files.append(filepath)
                # print(f"{op.basename}: possible input file.")
    
    return [by_dim, input_files]



def time_units_from_input(filepaths: list) -> TimeUnits:
    '''
    Attempt to find time unit data in possible input file(s).
    Refers to FROM_INPUT_FILE['reftime'] from setup.py, for variable name
    under which time units are expected.
    '''

    input_data = {}
    for f in filepaths:
        with xr.open_dataset(f, decode_times=False) as ds:  
            # , concat_characters=False
            try:
                input_data[f] = ds[FROM_INPUT_FILE['reftime']].attrs
            except KeyError:
                # Variable not found in prospective input file; move on to
                # next input file, if any.
                continue
    if len(input_data) != 1:
        # Request user input for tie-breaking.
        if len(input_data) > 1:
            print('Multiple time units found:')
            for i, u in enumerate(input_data.values()):
                print(f'[{i}] {u}')
            while True:
                user_input = input('Please select number or enter units for time, including reference date in ISO format.')
                if user_input.isnumeric():
                    try:
                        input_data = list(input_data.values())[int(user_input)]
                    except:
                        print('Invalid number; try again.')
                        continue
                    else:
                        break
                else:
                    input_data = {'units': user_input}
        else:
            input_data = {
                'units': 
                input('Enter units for time, including reference date in ISO format:')
            }
    else:
        input_data = list(input_data.values())[0]
    while True:
        try:
            base_unit, ref_date = re.split(' since ', input_data['units'])
            try:
                ref_date = datetime.fromisoformat(ref_date)
            except:
                input_data['units'] = input(e + 'Enter units for time, including reference date in ISO format:')
                continue
            if not ref_date.tzinfo:
                ref_date = ref_date.replace(tzinfo=timezone.utc)
            input_data['units'] = f"{base_unit} since {ref_date.isoformat(sep=' ')}"

            # BUG: is there a way to persuade cfunits to use the time zone from
            # the reference date? It uses cftime.datetime type for reftime,
            # which doesn't allow for time zone or offset.
            time_units = TimeUnits(
                units=input_data['units'], 
                calendar=input_data['calendar'] if 'calendar' in input_data else DEFAULT_CALENDAR
            )
        except ValueError or KeyError as e:
            input_data['units'] = input(e + 'Enter units for time, including reference date in ISO format:')
            continue
        else:
            return time_units
            
            
def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'source_dir',
        help='Directory containing MONC output NetCDF files to process.'
    )
    parser.add_argument(
        '--target_dir', '-t', 
        dest='target_dir', 
        help='Directory for processed NetCDF files.',
        required=False
    )
    parser.add_argument(
        '--reference_time', '-r',
        dest='ref_time',
        help='Date or date-time of origin for time units, in ISO format (yyyy-mm-dd[{T| }hh:mm[:ss][{+|-}hh:mm]]).',
        required=False
    )
    parser.add_argument(
        '--calendar', '-c',
        dest='calendar',
        choices={'gregorian', 'standard', 'none', 'proleptic_gregorian',
                 '360_day', 'noleap', '365_day', 'all_leap', '366_day',
                 'julian'},
        help='Calendar to use for time units.',
        required=False
    )
    parser.add_argument(
        '--cpus', '-p',
        dest='n_proc',
        type=int,
        choices=range(1, 17),
        help='Total number of cores to use for parallel processing (including controller process).',
        required=False
    )

    parser.add_argument(
        '--keep-interim', '-i',
        dest='keep_interim',
        action='store_true',
        # type=tf,
        # choices={True, False},
        default=False,
        help='Set if multiple dimension sets are being merged, but you want to keep the single-dimension-set files as well.',
        required=False
    )

    parser.add_argument(
        '--verbose', '-v',
        dest='verbose',
        action='store_true',
        # type=tf,
        # choices=(True, False),
        default=False,
        help='Set to report all progress to stdout.',
        required=False
    )

    parser.add_argument(
        '--quiet', '-q',
        dest='quiet',
        action='store_true',
        # type=tf,
        # choices=(True, False),
        default=False,
        help='Set to suppress warnings being printed to stdout.',
        required=False
    )

    return parser.parse_args()

def main():
    start_time = perf_counter()  # wrapping in performance counter doesn't work, presumably because of multiprocessing.

    # get command line arguments
    args = parse_arguments()
    
    global quiet, verbose, vocabulary
    quiet = args.quiet
    verbose = args.verbose

    warnings = []

    if verbose: print(f"{strftime('%H:%M:%S', localtime())} "
                      f"Main app process id: {os.getpid()}")  # , started {start_time}


    # Validate supplied directory to parse.
    if not op.exists(args.source_dir):
        sys.exit(OSError(f'Source directory {op.abspath(source_dir)} not found.'))
    source_dir = op.abspath(args.source_dir)
    
    time_units = None
    # Validate & set reference time if supplied.
    if args.ref_time:
        try:
            ref_time = datetime.fromisoformat(args.ref_time)
            if not ref_time.tzinfo:
                ref_time = ref_time.replace(tzinfo=timezone.utc)
        except Exception as e:
            sys.exit(f"{args.ref_time} is an invalid reference time. {str(e)}")
    else:
        ref_time = None

    # Validate calendar and if possible set time units.
    calendar = DEFAULT_CALENDAR
    if args.calendar:
        try:
            calendar = Units(calendar=args.calendar).calendar
        except Exception as e:
            sys.exit(
                f"Calendar {args.calendar} specified in arguments is not valid."
            )  # This should not be possible
    if ref_time:
        try:
            time_units = TimeUnits(
                units=f"{default_time_unit} since {ref_time.isoformat(sep=' ')}",
                calendar=calendar)
        except ValueError as e:
            sys.exit(
                f"Could not create valid time units from reference date-time {ref_time} and calendar {calendar}. {str(e)}"
            )
    
    # Get/set directory for output files.
    target_dir = op.abspath(args.target_dir) if args.target_dir else op.abspath(
        op.join(op.dirname(source_dir), 
                f'{op.basename(source_dir)}+processed')
    )
    # If target directory doesn't exist, create it.
    if not op.exists(target_dir):
        try:
            os.makedirs(target_dir)
        except OSError as e:
            sys.exit(
                f"Unable to create target directory, {target_dir}. {str(e)}"
            )
    if not os.access(target_dir, os.W_OK):
        sys.exit(OSError(
            f"Write permission denied for target directory, {target_dir}."
        ))
    if not args.target_dir:
        if not quiet:
            print(
                f"{strftime('%H:%M:%S', localtime())} "
                f"Processed files will be saved to: {target_dir}"
            )

    # Set number of processes to place in process pool.
    # Subtract 1 from number of processes specified, to use one as controller.
    n_proc = args.n_proc - 1 if args.n_proc else 0
    if n_proc > os.cpu_count() - 1:
        sys.exit(OSError(
            f"Not enough cores available for {args.n_proc} processes. Maximum: {os.cpu_count()}."
        ))
        # TODO: check os.cpu_count() works correctly on JASMIN, when mutliple
        # cores are allocated.
    
    shared = {
        'target_dir': target_dir,
        'verbose': verbose,
        'quiet': quiet
    }   # These are global variables (not constants), but must be passed 
        # explicitly to any functions that may run on a separate process.

    # Parse directory & categorise NC files by type and dimensions
    try:
        [group_by_dim, input_files] = sort_nc(source_dir, shared=shared)
    except Exception as e:
        sys.exit('sort_nc: ' + str(e))
    
    # ID time variable of each group from VOCAB.
    
    # Attempt to derive time unit info, if not already supplied at command line.
    # Assume time units will be common to all datasets in directory.
    if not time_units:
        time_units = time_units_from_input(input_files) if input_files else None
        if shared['verbose'] and time_units:
            print(strftime('%H:%M:%S', localtime()), 'Time units found in input file:', (time_units.formatted()))

    # NOTE: by this stage, time_units should contain a valid reference date /
    # datetime and calendar. The base units will be set to the default 
    # 'seconds', but this will be overridden by any units found in time 
    # coordinate variable(s).
            
    shared.update({
        'time_units': time_units,
        'vocabulary': vocabulary,
        'reference_vars': reference_vars
    })   # These are global variables (not constants), but must be passed 
        # explicitly to any functions that may run on a separate process.
    
    # For each dimension group:
    # If n_proc>0, use process pool; otherwise process sequentially.
    try:
        if n_proc > 0:
            warnings += process_parallel(
                groups=group_by_dim, 
                n_proc=n_proc, 
                shared=shared
            )
        else:
            warnings += process_serial(
                groups=group_by_dim, 
                shared=shared
            )
    except (ConfigError or VocabError) as e:
        sys.exit('Processing: ' + str(e))
    except Exception as e:
        sys.exit('Processing: ' + str(e))

    # TODO: move the following to the process_* routines?
    # Find dimension groups to be merged
    groups_to_merge = {}
    for group_name in DIM_ACTIONS.keys():
        groups_to_merge[group_name] = [
            g for g in group_by_dim.values() if g.name == group_name
        ]
        if group_name in groups_to_merge and len(groups_to_merge[group_name]) <= 1:
            groups_to_merge.pop(group_name)
        
    # For each set of dimension groups to be merged:
    for name, groups in groups_to_merge.items():
        # Create new group, comprising the groups to be merged.
        merger = DsGroup(groups=groups, shared=shared)
        if merger.log:
            print(
                ''.join([entry + '/n' for entry in merger.log])
            )
        
        if merger.action == 'merge_groups':
            # Merge each group's resultant single dataset
            try:
                (merged_file, log) = merger.merge_groups(shared=shared)
            except xr.MergeError as e:
                sys.exit("Merge error in DsGroup.merge_groups: " + str(e))
            except Exception as e:
                sys.exit("Unexpected error in DsGroup.merge_groups: " + str(e))
            if log:
                print(''.join([entry + '\n' for entry in log]))
            
        group_by_dim[re.split(merger.stem,merger.name)[1]] = merger

        # Delete interim NC files, unless flagged to do otherwise.
        # TODO: this needs to wait until all 3d are done.
        if not args.keep_interim:
            try:
                [os.remove(f) for f in merger.filepaths]
            except OSError as e:
                if not shared['quiet']:
                    print(strftime('%H:%M:%S', localtime()), "Failed to delete interim NC files. " + str(e))
            else:
                if shared['verbose']:
                    print(
                        f"{strftime('%H:%M:%S', localtime())} Removed interim "
                        f"files:\n        ", 
                        "\n         ".join([
                            op.basename(f) for f in merger.filepaths
                        ])
                    )

    # Output list of actions taken: each list of merged files & what file they were merged into; each split file and list of files it was split into; each file processed without merge/split & what its new version is called.
    print("\nSUMMARY\n=======")
    for k, group in group_by_dim.items():
        print(
            f"Group {k}: source files "
            f"{[op.basename(path) for path in group.filepaths]} --> "
            f"{group.action + ' --> ' if group.action and group.action != 'merge_groups' else ''}cfize --> "
            f"{[op.basename(path) for path in group.processed]}"
        )
    print()
    # Garbage collection if required.

    end_time = perf_counter()
    if verbose: print(f"         Main app process {os.getpid()} took "
                      f"{end_time - start_time} seconds")

    if len(warnings) > 0:
        if shared['quiet']:
            sys.exit('Quiet mode: warnings suppressed.')
        else:
            sys.exit('\n'.join([str(e) for e in warnings]))
    else:
        sys.exit(0)


if __name__ == '__main__': main()