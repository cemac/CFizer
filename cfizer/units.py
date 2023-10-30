from cfunits import Units
import re
from typing import Union


def format_units(units: Union[Units, str]) -> str:
    """
    cfunits.Units.formatted() reorders component units weirdly.
    This formatter preserves the original ordering, but makes units look more
    CF-like.
    It also preserves ratios in their original form, e.g. 'kg kg-1' instead of
    '1'.
    """

    if isinstance(units, Units):
        units = units.units

    frac = units.split("/")
    unit_list = [u.replace("^", "") for u in re.split("[ .]", frac[0])]

    if len(frac) > 1:
        for denom in frac[1:]:
            for u in re.split("[ .]", denom):
                if "^" in u:
                    unit_list.append("-".join(u.split("^")))
                else:
                    unit_list.append(f"{u}-1")

    return " ".join(unit_list)


class TimeUnits(Units):
    """
    Uses proleptic_gregorian as default calendar, as per ISO 8601:2004,
    even though CF standard calendar is a mixed Gregorian/Julian.
    """

    def __init__(
        self,
        units: str = None,
        calendar: str = "proleptic_gregorian",
        formatted=False,
        names=False,
        definition=False,
        _ut_unit=None,
    ):
        super().__init__(
            units=units,
            calendar=calendar,
            formatted=formatted,
            names=names,
            definition=definition,
            _ut_unit=_ut_unit,
        )
        if not self.isvalid:
            raise ValueError("Units not valid according to cfunits module.")

    def time_unit(self) -> str:
        return self.formatted().split(" since ")[0]

    def since(self) -> str:
        return self.reftime.isoformat()

    def cf(self) -> str:
        return self.formatted()

    def base_units_match(self, other: Units) -> bool:
        return (
            self.formatted().split(" since ")[0]
            == other.formatted().split(" since ")[0]
        )

    def base_units_equivalent(self, other) -> bool:
        return Units(units=self.formatted().split(" since ")[0]).equivalent(
            Units(units=other.formatted().split(" since ")[0])
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
            print("calendars_match function:", e)
            return False
