from spacesim import constants as const
import datetime as dt


def utc_to_julian_date(epoch: dt.datetime) -> float:
    """Converts an epoch in UTC time to a Julian Date

    Args:
        epoch (dt.datetime): The epoch to convert

    Returns:
        float: A julian date.
    """
    y = epoch.year
    m = epoch.month
    d = epoch.day
    H = epoch.hour
    M = epoch.minute
    S = epoch.second
    
    J_0 = 367 * y - int((7 * (y + int((m + 9) / 12))) / 4) + int((275 * m) / 9) + d + 1721013.5
    hour_frac = H + (M / 60.0) + (S / 3600.0)
    day_frac = hour_frac / 24.0
    
    JD = J_0 + day_frac
    return JD

def gmst(epoch: dt.datetime) -> float:
    """Calculates the Greenwich Mean Sidereal Time (GMST) at the given epoch.
    
    Source: https://astronomy.stackexchange.com/questions/21002/how-to-find-greenwich-mean-sideral-time,
            The Explanatory Supplement to the Astronomical Almanac, 3rd ed. (1992), p. 50 
    
    Args:
        epoch (datetime): The epoch to calculate GMST for in UTC time.
    
    Returns:
        float: The GMST in seconds.
    """
    midnight = dt.datetime(epoch.year, epoch.month, epoch.day, 0, 0, 0)
    jd_delta: float = utc_to_julian_date(midnight) - utc_to_julian_date(const.J2000)
    
    T_u: float = jd_delta / 36525
    H_0: float = 24110.54841 + 8640184.812866 * T_u + 0.093104 * T_u**2 - 6.2e-6 * T_u**3      # GMST at midnight
    omega: float = 1.00273790935 + 5.9e-11 * T_u                                 # Earth's sidereal rotation rate

    day_seconds: int = epoch.hour * 3600 + epoch.minute * 60 + epoch.second
    H: float = H_0 + omega * day_seconds           # GMST in seconds
    
    return H % const.solar_day_length

def tle_epoch_to_datetime(julian: str) -> dt.datetime:
    """Converts a Julian date in the form `YYDDD.DDDDDDDD` into datetime format,
    where `YY` are last two digits of the year and `DD...` is the day fraction.
    
    Args:
        julian (str): A Julian date.

    Returns:
        datetime: The corresponding date time object.
    """
    julian_year = int(julian[0:2])
    year: int = 2000 + julian_year if julian_year < 57 else 1900 + julian_year
    
    day_fraction = float(julian[2:])
    seconds_in_day: int = 24 * 3600
    # TLE epoch starts at one, not zero.
    julian_seconds: float = day_fraction * seconds_in_day - seconds_in_day
    day_of_year = dt.timedelta(seconds=julian_seconds)
    
    year_start = dt.datetime(year, 1, 1)
    return year_start + day_of_year

def datetime_nearest_sec(time: dt.datetime) -> dt.datetime:
    """Rounds a datetime object to the nearest second.
    
    Args:
        time (datetime): The datetime object to round.

    Returns:
        datetime: The rounded datetime object.
    """
    rnd_sec = 1 if time.microsecond >= 500000 else 0
    seconds = time.second + rnd_sec
    
    return dt.datetime(time.year, time.month, time.day, time.hour, time.minute, seconds)

if __name__ == "__main__":
    date = dt.datetime(2023, 9, 14, 0, 0, 0)
    print(gmst(date))
    print(utc_to_julian_date(const.J2000))