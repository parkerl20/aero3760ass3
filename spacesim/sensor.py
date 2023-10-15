from spacesim import orbit as orb
from typing import Callable


class SatelliteSensor():
    """A class that can simulate a sensor on a satellite.
    """
    def __init__(
        self,
        simulator: Callable[['SatelliteSensor', orb.Orbit], None],
        frequency: float = None
    ) -> None:
        """A sensor that works on a satellite

        Args:
            simulator (Callable[[SatelliteSensor, orb.Orbit], any]): A function that
                simulates the sensor making measurements.
            frequency (float): The operating frequency of the sensor. 
                Defaults to None. 
        """
        self.simulator = simulator
        self.mesurement: any = None
        self.frequency = frequency
        
        return

    def make_measurement(self, measurement: any) -> None:        
        self.mesurement = measurement
        return
    
    def get_measurement(self) -> any:
        return self.mesurement
    
    def sampling_period(self) -> float:
        return 1 / self.frequency