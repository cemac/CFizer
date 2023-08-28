from cfunits import Units


class TimeUnits(Units):
    '''
    Uses proleptic_gregorian as default calendar, as per ISO 8601:2004,
    even though CF standard calendar is a mixed Gregorian/Julian.
    '''
    def __init__(self, units:str=None, calendar:str='proleptic_gregorian', formatted=False, names=False, definition=False, _ut_unit=None):
        super().__init__(units=units, calendar=calendar, formatted=formatted, names=names, definition=definition, _ut_unit=_ut_unit)
        if not self.isvalid:
            raise ValueError('Units not valid according to cfunits module.')
    
    def time_unit(self) -> str:
        return self.formatted().split(' since ')[0]
    
    def since(self) -> str:
        return self.reftime.isoformat()
    
    def cf(self) -> str:
        return self.formatted()
    
    def base_units_match(self, other: Units) -> bool:
        return self.formatted().split(' since ')[0] == other.formatted().split(' since ')[0]
    
    def base_units_equivalent(self, other) -> bool:
        return Units(units=self.formatted().split(' since ')[0]).equivalent(
            Units(units=other.formatted().split(' since ')[0])
            )
    
    def ref_dates_match(self, other: Units) -> bool:
        return self.reftime == other.reftime
    
    def has_calendar(self) -> bool:
        try:
            return self.calendar is not None
        except AttributeError:
            return False

    def calendars_match(self, other: Units) -> bool:
        try:
            return self.calendar == other.calendar
        except AttributeError:
            return False
        except Exception as e:
            print('calendars_match function:', e)
            return False
