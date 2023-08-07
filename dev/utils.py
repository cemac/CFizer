import re
from datetime import timedelta, datetime

def init_vocab():
    pass

def options_db_to_global_attrs():
    pass


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


if __name__ == '__main__':
    value = 315960.
    units = 'seconds since 2020-01-25'
    decoded_units = decode_time_units(units)
    print('decode_time_units:', decoded_units)
    print('timedelta_from_units:', timedelta_from_units(value, decoded_units[0]))
    print('datetime_from_relative_units:', datetime_from_relative_units(value, units))
