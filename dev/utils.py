import re
from time import perf_counter
from datetime import timedelta, datetime
from cfunits import Units


def performance_time(func):
	def wrapper(*args, **kwargs):
		start_time = perf_counter()
		response = func(*args, **kwargs)  # run the wrapped function
		end_time = perf_counter()
		duration = end_time - start_time
		print(f"{func}({args}, {kwargs}) took {duration} seconds.")
		return response
	return wrapper


def vocab_from_xls(filepath: str) -> dict:
    '''
    This will for Excel file to convert to vocabulary dictionary, if it 
    contains the right columns (fields).
    It will attempt to fill any missing data based on CF conventions, including:
    - inferring cell_method from presence of terms such as "mean" and "max" in variable name.
    - looking up standard_name, if supplied.
    It will check any units given using `cfunits`, and check against CF standard names definitions.

    The intention is to compile a full lookup table, such that the main 
    application does not need to infer anything or look anywhere else for the 
    required data.
    '''
    import openpyxl as xl

    VOCABULARY = {}
    # 
    return VOCABULARY


def options_db_to_global_attrs():
    pass


def is_sequence(object):
    return isinstance(object, (
        list,
        set,
        tuple,
        range
    ))


def is_sequence_of(object, type: type):
    if is_sequence(object):
        return all([isinstance(element, type) for element in object])
    else:
        return False


def type_from_str(string: str):
    if not isinstance(string, str):
        # Can't work with it, so return as is.
        return string
    if string.lower() == 'true':
        return True
    elif string.lower() == 'false':
        return False
    else:
        try:
            float(string)
            # Note: string.isnumeric() returns false for numbers containing 
            # decimal points or scientific notation.
        except ValueError:
            # not a number
            return string  # leave as a string
        else:
            if '.' not in string:
                try:
                    int(string)
                except ValueError:
                    # not an integer;
                    # should never reach this option.
                    return float(string)
                else:
                    return int(string)
            else:
                return float(string)  # floating point, whether integer or not


def decode_time_units(units: str) -> tuple:
    try:
        (unit, ref_datetime) = units.lower().split(' since ')
    except AttributeError or ValueError:
        raise ValueError('''Argument must be a string in the format:
                        <time-unit> since <time-origin>,
                         with <time-origin> in ISO format.''')
    
    try:
        dt = datetime.fromisoformat(ref_datetime)
    except:
        raise ValueError('Origin date-time must be in ISO format.')
    
    return (unit, dt)


def timedelta_from_units(number: float|int, unit: str) -> timedelta:
    
    try:
        float(number)
    except ValueError:
        raise TypeError('First argument must be numeric')

    time_units = dict(
        days=0,
        seconds=0, 
        microseconds=0, 
        milliseconds=0, 
        minutes=0, 
        hours=0, 
        weeks=0
    )

    if unit not in time_units:
        raise KeyError('Invalid unit supplied')
    
    time_units[unit] = number

    return timedelta(
        days=time_units['days'],
        seconds=time_units['seconds'],
        microseconds=time_units['microseconds'],
        milliseconds=time_units['milliseconds'],
        minutes=time_units['minutes'],
        hours=time_units['hours'],
        weeks=time_units['weeks']
    )

def datetime_from_relative_units(number: float|int, units: str) -> datetime:
    
    try:
        (unit, ref_datetime) = decode_time_units(units)
    except TypeError:
        raise TypeError('''Second argument must be a string in the format:
                        [time unit] since [date or date-time in ISO format]''')

    delta = timedelta_from_units(number, unit)

    return ref_datetime + delta


def generate_coords(number: int, 
                    spacing: float|int, 
                    midpoint: bool = False) -> list:
    first = spacing/2 if midpoint else 0
    return [(x * spacing + first) for x in range(number)]


def stem_str(*args):
        if not args:
            raise ValueError('At least one string must be provided.')
        filenames = [os.path.splitext(os.path.basename(path))[0] for path in args]  # If each path contains only a filename (no extension) this will still work.
        stem = filenames[0]
        if len(args) > 1:
            for f in filenames[1:]:
                match = SequenceMatcher(a=stem, b=f).find_longest_match()
                stem = stem[match.a: match.a + match.size]
        return stem
    

if __name__ == '__main__':
    value = 315960.
    units = 'seconds since 2020-01-25'
    decoded_units = decode_time_units(units)
    print('decode_time_units:', decoded_units)
    print('timedelta_from_units:', timedelta_from_units(value, decoded_units[0]))
    print('datetime_from_relative_units:', datetime_from_relative_units(value, units))
