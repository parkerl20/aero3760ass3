from spacesim import satellite as sat
from typing import Any, Callable
import numpy as np


class SatelliteSensor():
    """A class that can simulate a sensor on a satellite.
    
    The measurement field is initialized to None and must 
    be set by the sensor's simulator.
    """
    def __init__(
        self,
        name: str,
        simulator: Callable[
            ['SatelliteSensor', 'sat.RealTimeSatellite', np.ndarray, np.ndarray],
            None
        ],
        frequency: float = None
    ) -> None:
        """A sensor that works on a satellite

        Args:
            name (str): The name of the sensor.
            simulator (Callable[[SatelliteSensor, orb.Orbit, np.ndarray, np.ndarray], bool]): A function that
                simulates the sensor making measurements. Returns True if the sensor
                makes a measurement, and False if it does not. Takes the sensor, satellite,
                position, and velocity as arguments.
            frequency (float): The operating frequency of the sensor. 
                Defaults to None. 
        """
        self.name = name
        self.simulator = simulator
        self.mesurement: any = None
        self.frequency = frequency
        
        return

    def make_measurement(
        self,
        satellite: sat.RealTimeSatellite,
        r: np.ndarray,
        v: np.ndarray
    ) -> bool:        
        return self.simulator(self, satellite, r, v)
    
    def get_measurement(self) -> any:
        return self.mesurement
    
    def sampling_period(self) -> float:
        if self.frequency is None:
            return 0
        else:
            return 1 / self.frequency