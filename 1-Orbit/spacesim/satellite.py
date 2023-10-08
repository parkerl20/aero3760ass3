from spacesim import orbit as orb
from spacesim import celestial_body as cb
from spacesim import time_util as tu

from typing import Callable
import datetime as dt
import numpy as np
import math


class Satellite():
    """A class to represent a planet satellite.
    """
    
    def __init__(
        self,
        name: str,
        tle_1: str,
        tle_2: str,
        primary_body: cb.CelestialBody,
        orbit_dynamics: Callable[[float, np.ndarray, orb.Orbit], np.ndarray] = None
    ) -> None:
        """Initializes a planet satellite object from a two line element.

        Args:
            name (str): The name of the satellite.
            tle_1 (str): The first line of the two line element.
            tle_2 (str): The second line of the two line element.
            planet (cb.Planet, optional): The planet the satellite is orbiting.
            orbit_dynamics (Callable, optional): A function that simulates the orbit dynamics
            of the satellite. Must be in the form `f(t: float, v: np.ndarray, orbit: Orbit) -> np.ndarray`.
        """
        self.name = name
        self.satellite_number: int = int(tle_1[2:7])
        self.intl_designator: str = tle_1[9:17]
        self.launch_year: int = int(tle_1[9:11]) + 2000
        epoch: dt.datetime = tu.tle_epoch_to_datetime(tle_1[18:32])         # in UTC time
        
        i: float = float(tle_2[8:16])
        rt_asc: float = float(tle_2[17:25])
        e: float = float("0." + tle_2[26:33])
        arg_p: float = float(tle_2[34:42])
        
        mean_anomaly: float = float(tle_2[43:51])
        mean_motion: float = float(tle_2[52:63])
        
        seconds_in_day: int = 24 * 3600
        n = (2 * math.pi * mean_motion) / seconds_in_day        # mean motion in rad/s
        mu = primary_body.gravitational_parameter
        
        a = (mu / n**2) ** (1/3)
        theta = orb.Orbit.mean_to_true_anomaly(mean_anomaly, e)
        
        self.orbit = orb.Orbit(a, e, i, rt_asc, arg_p, theta, primary_body, epoch, self.name, orbit_dynamics)
        
        return

    def __str__(self) -> str:
        right_justify = 21
        
        sat_num_str = f"Satellite Number".rjust(right_justify)
        intl_des_str = f"Intl. Designator".rjust(right_justify)
        
        return (
            f"{sat_num_str}: {self.satellite_number}\n"
            f"{intl_des_str}: {self.intl_designator}\n"
            f"{self.orbit}"
        )