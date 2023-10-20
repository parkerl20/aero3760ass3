from spacesim import orbit as orb
from spacesim import celestial_body as cb
from spacesim import time_util as tu
from spacesim import sensor
from spacesim import util

from typing import Any, Callable
from scipy import integrate
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

class SatelliteAlgorithm():
    """An algorithm that can run on a satellite.
    
    An algorithm can be any generic object that acts 
    on a satellite's state and/or sensor data.
    
    The function that simulates the algorithm must take
    the algorithm, the satellite, and a dictionary of
    sensors and their measurements as arguments. If a
    sensor did not make a measurement since the last
    time the algorithm was run, then it will not be
    included in the dictionary.
    """
    def __init__(
        self,
        name: str,
        algorithm: any,
        algo_func: Callable[
            [
                float,
                'SatelliteAlgorithm',
                'RealTimeSatellite',
                dict
            ], 
            None
        ],
        logger: util.DictLogger = None
    ) -> None:
        self.name = name
        self.algorithm = algorithm
        self.algo_func = algo_func
        self.t_last = 0         # the last time the algorithm was run
        self.logger = logger
        return
    
    def __call__(
        self,
        t: float,
        satellite: 'RealTimeSatellite',
        measurements: dict[str, sensor.SatelliteSensor]
    ) -> None:
        """Runs the algorithm on the satellite.

        Args:
            t (float): The current time.
            satellite (RealTimeSatellite): The satellite to run the algorithm on.
            measurements (dict[str, sensor.SatelliteSensor]): A dictionary of
                recent measurements from the satellite's sensors. If a sensor
                has not made a measurement since the last time the algorithm
                was run, then it will not be included in the dictionary.
        """
        self.algo_func(t, self, satellite, measurements, self.logger)
        self.t_last = t
        return
      
    
class RealTimeSatellite(orb.Orbit):
    """An orbit that can be iterated through in real time.
    
    A real time satellite can carry sensors.
    A real time satellite can perform algorithms
    that act on the sensor data.
    """
    def __init__(
        self,
        a: float,
        e: float,
        i: float,
        rt_asc: float,
        arg_p: float,
        true_anomaly: float,
        body: cb.CelestialBody,
        epoch: dt.datetime,
        *,
        name: str = None,
        orb_dyn: Callable[[float, np.ndarray, 'RealTimeSatellite'], np.ndarray] = None,
        colour: str = None,
        propagation_length: float = 20.0,
        propagation_step: float = 1.0
    ) -> None:
        """Initializes a RealTimeOrbit object.
        
        The orbit related arguments define the initial position of
        the orbit.
        
        Args:
            a (float): The semi-major axis of the orbit in meters.
            e (float): The eccentricity of the orbit.
            i (float): The inclination of the orbit in degrees.
            rt_asc (float): The right ascension of the orbit in degrees.
            arg_p (float): The argument of perigee of the orbit in degrees.
            theta (float): The true anomaly of the orbit in degrees.
            body (cb.Planet): The body being orbitted.
            epoch (dt.datetime): The epoch of the start of the orbit.
            name (str, optional): The name of the orbit. Defaults to None.
            orb_dyn (Callable, optional): A function that simulates the 
                dynamics of the orbit.
            If it is None, then a default function is used. Defaults to None.
            colour (str, optional): The colour of the orbit for plotting. 
                Defaults to None.
        """        
        super().__init__(a, e, i, rt_asc, arg_p, true_anomaly, body, epoch, name, orb_dyn, colour)
        
        self.current_r_eci = self.init_r_eci
        self.current_v_eci = self.init_v_eci
        self.current_time = 0
        
        self.propagation_length = propagation_length
        self.propagation_step = propagation_step
        
        # Dictionary of sensor and the last time they made a measurement
        self.sensors: dict[str, list[sensor.SatelliteSensor, float]] = dict()
        self.algorithms: dict[str, SatelliteAlgorithm] = dict()
        
        self._start_flag = True
        
        return
    
    def attach_sensor(self, sensor: sensor.SatelliteSensor) -> None:
        # -1 causes the sensor to make a measurement on the first iteration
        sensor_info = [sensor, -1]
        self.sensors[sensor.name] = sensor_info
        return
    
    def add_algorithm(self, algorithm: SatelliteAlgorithm) -> None:
        self.algorithms[algorithm.name] = algorithm
        return
    
    def get_sensor(self, name: str) -> sensor.SatelliteSensor:
        return self.sensors[name][0]
    
    def inc_propagate(
        self
    ) -> np.ndarray:
        """Propagates the orbit of a satellite.

        Returns:
            np.ndarray: The position of the satellite in ECI coordinates.
        """
        y0 = np.concatenate((self.current_r_eci, self.current_v_eci)).flatten()
        t_span = [0, self.propagation_length]
        
        # Calculate evaluation points      
        t_steps = np.arange(0, self.propagation_length + self.propagation_step, self.propagation_step)
        for sensor, t_measure in self.sensors.values():
            t_diff = self.current_time - t_measure
            init_sample = max(sensor.sampling_period() - t_diff, 0)
            # Time at which the sensor will make a measurement
            t_sample = np.arange(
                init_sample,
                self.propagation_length,
                sensor.sampling_period()
            )
            
            t_steps = np.union1d(t_steps, t_sample)
        
        solution = integrate.solve_ivp(
            self.orbit_dynamics,
            t_span,
            y0,
            t_eval=t_steps,
            args=(self,)
        )
        
        r = solution.y[0:3]
        v = solution.y[3:6]       
        t = solution.t + self.current_time
        
        # Simulate sensor measurements
        for i in range(len(t) - 1):
            position = r[:, i].flatten()
            velocity = v[:, i].flatten()
            time = t[i]
            
            self.current_r_eci = position
            self.current_v_eci = velocity
            
            sensor_measurements = dict()
            
            for sensor_name in self.sensors.keys():
                sensor, t_measure = self.sensors[sensor_name]
                t_diff = time - t_measure
                
                # Make measurement
                if t_diff >= sensor.sampling_period():
                    if sensor.make_measurement(self, position, velocity):
                        # print(f"IMU measurement at:\t{time}")
                        self.sensors[sensor.name][1] = time
                        sensor_measurements[sensor.name] = sensor
                
            # Run algorithms
            if len(sensor_measurements) > 0:
                for algorithm in self.algorithms.values():
                    algorithm(
                        time,
                        self,
                        sensor_measurements
                    )
                    np.set_printoptions(linewidth=1000, threshold=np.inf, precision=3)
                    # Test
                    # print(f"Algorithm: {algorithm.name}")
                    # print(f"Time: {time:.3f} s")
                    # print(f"State: {algorithm.algorithm.get_state().flatten()}")
                    
                    # pos_residual = r[:,i].flatten() - algorithm.algorithm.get_state()[:3].flatten()
                    # vel_residual = v[:,i].flatten() - algorithm.algorithm.get_state()[3:].flatten()
                    
                    # print(f"Position residual: {pos_residual.flatten()} m")
                    # print(f"Velocity residual: {vel_residual.flatten()} m/s")
                    
                    
        self.current_r_eci = r[:, -1].flatten()
        self.current_v_eci = v[:, -1].flatten()
        self.current_time = t[-1]
        
        return r[:,-1], v[:,-1], t[-1]
    
    def reset(self) -> None:
        """Resets the satellite to its initial state.
        """
        self.current_r_eci = self.init_r_eci
        self.current_v_eci = self.init_v_eci
        self.current_time = 0
        
        for sensor in self.sensors.values():
            sensor[1] = -1
        
        return
    
    def __iter__(self):
        return self
    
    def __next__(self) -> tuple[np.ndarray, np.ndarray, float]:
        if self._start_flag:
            self._start_flag = False
            return (
                self.current_r_eci.flatten(),
                self.current_v_eci.flatten(),
                self.current_time
            )
        
        return self.inc_propagate()