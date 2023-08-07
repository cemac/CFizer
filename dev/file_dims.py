import asyncio
import xarray as xr
from glob import glob
import os
import psutil
from time import perf_counter


def performance_time(func):
	def wrapper(*args, **kwargs):
		start_time = perf_counter()
		response = func(*args, **kwargs)  # run the wrapped function
		end_time = perf_counter()
		duration = end_time - start_time
		print(f"{func}({args}, {kwargs}) took {duration} seconds.")
		return response
	return wrapper


async def async_open_ds(filepath):
     print(f't={perf_counter()}: opening {filepath}...')
     ds = xr.open_dataset(filepath)
     print(f't={perf_counter()}: loaded {filepath}')
     return (filepath, ds)


async def get_n_dims(nc_filepath):
    print('opening', nc_filepath, 't=', perf_counter())
    with xr.open_dataset(nc_filepath) as ds:
        print('parsing dimensions, t=',perf_counter())
        n_dims = max(
            [len(ds[v].dims) for v in ds.variables if v != 'options_database']
        )
    print('done with dataset, t=', perf_counter())
    # n_dims = 2
    return (nc_filepath, n_dims)


async def categorise_by_dims(directory):
    dim_dict = {
        '0+1d': [],
        '2d': [],
        '3d': []
    }
    file_list = glob(
        pathname=f'{directory}/*.nc', recursive=False
    )
    tasks = [get_n_dims(f) for f in file_list]
    
    for future in asyncio.as_completed(tasks):
        data = await future
        if data[1] < 3:
            dim_dict['0+1d'].append(data[0])
        elif data[1] < 4:
            dim_dict['2d'].append(data[0])
        else:
            dim_dict['3d'].append(data[0])
    print(dim_dict)


async def nc_by_dims(directory):
    dim_dict = {
        '0+1d': [],
        '2d': [],
        '3d': []
    }
    file_list = glob(
        pathname=f'{directory}/*.nc', recursive=False
    )
    print(f't={perf_counter()}: setting up tasks...')
    loading_nc = [async_open_ds(f) for f in file_list]
    for future in asyncio.as_completed(loading_nc):
        print(future)
        (f, ds) = await future
        print(f't={perf_counter()}: parsing dimensions of {f}...')
        n_dims = max(
            [len(ds[v].dims) for v in ds.variables if v != 'options_database']
        )
        ds.close()
        print(f't={perf_counter()}: done with dataset')
        if n_dims < 3:
            dim_dict['0+1d'].append(f)
        elif n_dims < 4:
            dim_dict['2d'].append(f)
        else:
            dim_dict['3d'].append(f)
    print(f't={perf_counter()}: Done.')
    print(dim_dict)


if __name__ == '__main__':
    test_dir = os.path.join(os.path.dirname(os.getcwd()), 'test_data')
    # loop = asyncio.get_event_loop()
    # asyncio.run(categorise_by_dims(test_dir))
    asyncio.run(nc_by_dims(test_dir))
    # file_categories = loop.run_until_complete(categorise_by_dims(directory=test_dir))
    # print(file_categories)