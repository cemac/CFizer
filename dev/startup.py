'''
Name:           cfizer.py
Author:         Cameron Wilson, CEMAC, University of Leeds
Date:           August 2023
CFizer Version: See __version__

'''

# import yaml
# import os
# from datetime import datetime
# import numpy as np
from __init__ import *  # If init takes care of once-only declarations, this module allows them to be imported en masse to other routines, without their being rerun, I think.
# from utils import vocab_from_xls

# COMPRESSION = (True, 'zlib', 2) # True/False, compression type, complevel (1:9)
#     # See https://unidata.github.io/netcdf4-python/efficient-compression-of-netcdf-variables

# DIM_GROUPS = {
#     0: '0+1d', 
#     1: '0+1d', 
#     2: '2d', 
#     3: '3d'
# }

# DIM_ACTIONS = {
#     '0+1d': 'merge',
#     '2d': 'merge', 
#     '3d': 'split'
#     }

# MONC_ID_ATTR = 'MONC timestep'  # It is essential that this be preserved across all modified datasets, either as attribute or variable.

# OPTIONS_DATABASE = {
#     'variable': 'options_database',
#     'dimensions': {'number_options', 'kvp', 'string'}
# }

# FROM_INPUT_FILE = {
#     'reftime': 'time'
# }

# CHUNKING_DIMS = {'z': 2, 'zn':2}  # Ideally chunk sizes should be power of 2.

# DEPHY_OPTIONS = ('dephy_file',)  # 'dephy_forcings_enabled'
# DROP_FOR_DEPHY = ("longitude", "latitude", "z0")

# Global attributes and their Numpy data types.
# do_not_propagate = {'MONC timestep': np.int32}
# split_attrs = {
#     'MONC time': (
#         lambda i, ds_list: ds_list[i]['time'].data.tolist(), 
#         np.float64
#     ),
#     'Previous diagnostic write at': (
#         lambda i, ds_list: ds_list[0].attrs['Previous diagnostic write at'] 
#         if i == 0 else ds_list[i - 1]['time'].data.tolist(), 
#         np.float64
#     )
# }
# group_attrs = {'created': datetime}


