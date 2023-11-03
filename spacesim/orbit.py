from spacesim import celestial_body as cb
from spacesim import orbital_transforms as ot
from spacesim import time_util as tu

from typing import Callable
import scipy.integrate as integrate
import scipy.optimize as optimize
import numpy as np
import datetime as dt
import math


class Orbit():
    """An orbit is defined by the six classical orbital elements.
    
    Attributes:
        a (float): The semi-major axis of the orbit in meters.
        e (float): The eccentricity of the orbit.
        i (float): The inclination of the orbit in degrees.
        rt_asc (float): The right ascension of the orbit in degrees.
        arg_p (float): The argument of perigee of the orbit in degrees.
        true_anomaly (float): The true anomaly of the orbit in degrees
        period (float): The period of the orbit in seconds.
        body (cb.CelestialBody): The body being orbitted.
        epoch (dt.datetime): The epoch of the start of the orbit.
        name (str): The name of the orbit.
        orbit_dynamics (Callable): A function that simulates the dynamics of the orbit.
    """
    
    def __init__(self,
        a: float,
        e: float,
        i: float,
        rt_asc: float,
        arg_p: float,
        theta: float,
        body: cb.CelestialBody,
        epoch: dt.datetime,
        name: str = None,
        orb_dyn: Callable[[float, np.ndarray, "Orbit"], np.ndarray] = None,
        colour: str = None
    ) -> None:
        """Initializes an orbit object.
        
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
            orb_dyn (Callable, optional): A function that simulates the dynamics of the orbit.
            If it is None, then a default function is used. Defaults to None.
            colour (str, optional): The colour of the orbit for plotting. Defaults to None.
        """
        self.body: cb.CelestialBody = body
        self.semi_major_axis: float = a
        self.eccentricity: float = e
        self.inclination: float = i
        self.right_ascension: float = rt_asc
        self.argument_periapsis: float = arg_p
        self.true_anomaly: float = theta
        self.epoch = epoch
        self.name = name
        self.colour = colour
        
        if orb_dyn is None:
            self.orbit_dynamics = Orbit.__default_dynamics
        else:
            self.orbit_dynamics = orb_dyn
        
        mu = self.body.gravitational_parameter
        
        self.init_r_eci, self.init_v_eci = ot.elements_to_ECI(
            a,
            e,
            i,
            rt_asc,
            arg_p,
            theta
        )
        
        self.init_attitude =np.array([[0, 0 ,0 ,0]]).T
        
        # Calculated instance variables
        self.period: float = 2 * math.pi * (mu / a**3) ** (-1/2)
        
        # Propergation cache
        self._propergate_period = [0, -1]
        self._propergate_step = -1
        self._prop_in_km = False
        self._frame = ""
        self._orbit_r: np.ndarray = None
        self._orbit_v: np.ndarray = None
        self._orbit_t: np.ndarray = None
        self._new_dynamics = False
        self._analytical = False
        
        return
    
    def __str__(self) -> str:
        right_justify = 21
        name_str = f"Name".rjust(right_justify)
        epoch_str = f"Epoch".rjust(right_justify)
        a_str = f"Semi-major axis".rjust(right_justify)
        e_str = f"Eccentricity".rjust(right_justify)
        i_str = f"Inclination".rjust(right_justify)
        rt_asc_str = f"Right ascension".rjust(right_justify)
        arg_p_str = f"Argument of periapsis".rjust(right_justify)
        theta_str = f"True anomaly".rjust(right_justify)
        period_str = f"Period".rjust(right_justify)
        
        
        return (f"{name_str}: {self.name}\n"
                f"{epoch_str}: {tu.datetime_nearest_sec(self.epoch)} UTC\n"
                f"{a_str}: {self.semi_major_axis:.3f} metres\n"
                f"{e_str}: {self.eccentricity:.4e}\n"
                f"{i_str}: {self.inclination:.3f} deg\n"
                f"{rt_asc_str}: {self.right_ascension:.3f} deg\n"
                f"{arg_p_str}: {self.argument_periapsis:.3f} deg\n"
                f"{theta_str}: {self.true_anomaly:.3f} deg\n"
                f"{period_str}: {self.period:.2f} secs\n")
    
    def propagate(
        self,
        t: float,
        *,
        t_start: float = 0,
        use_km: bool = False,
        coord_frame: str = "ECI",
        max_step: int = 20,
        analytical: bool = False,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Propagates the orbit of a satellite. 
        
        Checks the cache to see if the orbit has already been propagated 
        to the time specified. If not, then the orbit is propagated and 
        the results are cached.
        
        All units used and returned are SI.

        Args:
            t (float): The time to propagate to in seconds.
            t_start (float, optional): The time to start the propagation
                from in seconds.
            use_km (bool, optional): Whether to return the results in km.
                Defaults to False.
            coord_frame (str, optional): The coordinate frame to return the
                position in. Defaults to "ECI".
            max_step (int, optional): The maximum number of steps to take.
                Defaults to 20.

        Returns:
            tuple[np.ndarray, np.ndarray, np.ndarray]: The position, velocity
                and time steps arrays.
        """
        # TODO: add caching for analytical solution
        if analytical:  
            if (self._analytical 
                and self._propergate_period[0] == t_start 
                and self._propergate_period[1] == t
                and self._propergate_step == max_step
                and self._frame == coord_frame
            ):    
                r = self._orbit_r
                v = self._orbit_v
                t = self._orbit_t
                
                if use_km and not self._prop_in_km:
                    r /= 1000
                    v /= 1000
                elif not use_km and self._prop_in_km:
                    r *= 1000
                    v *= 1000
                    
                return r, v, t
            else:
                return self.__propagate_analytical(
                    t,
                    t_start,
                    use_km,
                    coord_frame,
                    max_step
                )
            
        
        # If no changes needed i.e. fully cached
        if self.__is_cached(t, t_start, use_km, coord_frame) and not analytical:
            return self._orbit_r, self._orbit_v, 0, self._orbit_t
        
        return self.__propagate(t, t_start, use_km, coord_frame, max_step)
    
    def __is_cached(
        self,
        t: float,
        t_start: float = 0,
        use_km: bool = False,
        coord_frame: str = "ECI"
    ) -> bool:
        """Private method to check if the orbit has already been propergated 
        to the time specified.

        Args:
            t (float): The time to propagate to in seconds.
            t_start (float, optional): The time to start the propagation 
                from in seconds.
            use_km (bool, optional): Whether to return the results in km. 
                Defaults to False.
            coord_frame (str, optional): The coordinate frame to return the 
                position in. Defaults to "ECI".

        Returns:
            bool: True if the orbit has already been propergated to the time 
                specified, else False.
        """
        new_prop_period = [t_start, t]
        return (
            all([a == b] for a, b in zip(self._propergate_period, new_prop_period))
            and self._prop_in_km == use_km
            and self._frame == coord_frame
            and not self._new_dynamics
        )
    
    @staticmethod
    def __default_dynamics(
        t: float,
        state_vector: np.ndarray, 
        orbit: "Orbit"
    ) -> np.ndarray:
        """Default function to simulate the dynamics of an orbit in a 
            restricted two body
        system. Works in conjuction with `scipy.integrate.solve_ivp`.

        Args:
            t (float): The current time in seconds.
            state_vector (np.ndarray): The current state of the orbit.
            orbit (Orbit): The orbit object being simulated.

        Returns:
            np.ndarray: The next state of the orbit.
        """
        mu = orbit.body.gravitational_parameter
        
        r = state_vector[0:3]
        v = state_vector[3:6]
        
        r_mag = np.sqrt(r.dot(r))
        
        # Dynamics
        r_dot = v
        v_dot = -(mu / (r_mag ** 3)) * r      # Two body equation
        
        # TODO: Attidude dynamics - flatten at end
        att = state_vector[6:10]       # 1D array length 4
        att_dot = np.array([[0, 0 ,0 ,0]]).T
        att_dot = att_dot.flatten()

        return [*r_dot, *v_dot, *att_dot]
    
    def __propagate(
        self,
        t: float,
        t_start: float = 0,
        use_km: bool = False,
        coord_frame: str = "ECI",
        max_step: int = 20
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Private method to propergate the orbit of a satellite.

        Args:
            t (float): The time to propagate to in seconds.
            coord_frame (str, optional): The coordinate frame to return 
                the position in. Defaults to "ECI".
            max_step (int, optional): The maximum number of steps to take. 
                Defaults to 20.

        Returns:
            tuple[np.ndarray, np.ndarray, np.ndarray]: The position, velocity
                and time steps arrays.
        """
        if coord_frame != "ECI" and coord_frame != "perifocal":
            raise ValueError(f'Invalid coordinate frame "{coord_frame}".')
        
        # if coord_frame == "ECI":
        #     init_r, init_v = ot.elements_to_ECI(
        #         self.semi_major_axis,
        #         self.eccentricity,
        #         self.inclination,
        #         self.right_ascension,
        #         self.argument_periapsis,
        #         self.true_anomaly
        #     )
        # elif coord_frame == "perifocal":
        #     init_r, init_v = ot.elements_to_perifocal(
        #         self.semi_major_axis,
        #         self.eccentricity,
        #         self.true_anomaly
        #     )
        
        init_r = self.init_r_eci
        init_v = self.init_v_eci
        init_att = self.init_attitude
        
        y0 = np.concatenate((init_r, init_v, init_att)).flatten()
        t_span = [t_start, t]
        solution = integrate.solve_ivp(
            self.orbit_dynamics,
            t_span,
            y0,
            max_step=max_step,
            args=(self,)
        )
        
        r = solution.y[0:3]
        v = solution.y[3:6]
        att = solution.y[6:10]
        t = solution.t
        
        if use_km:
            r = r / 1000
            v = v / 1000
        
        # Cache results
        self._orbit_r = r
        self._orbit_v = v
        self._orbit_t = t
        self._prop_in_km = use_km
        self._propergate_period = [t_start, t]
        self._frame = coord_frame

        return r, v, att, t

    def set_dynamics(
        self, 
        dynamics: Callable[[float, np.ndarray, "Orbit"], np.ndarray]
    ) -> None:
        """Sets the dynamics function to use when propergating the orbit.

        Args:
            dynamics (Callable[[float, np.ndarray, "Orbit"], np.ndarray]): 
                The dynamics function.
        """
        self.orbit_dynamics = dynamics
        self._new_dynamics = True
        return
    
    def __propagate_analytical(
        self,
        t: float,
        t_start: float = 0,
        use_km: bool = False,
        coord_frame: str = "ECI",
        max_step: int = 20,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Private method to propergate the orbit of a satellite 
        using an analytical solution. Suitable for pure keplerian orbits.

        Args:
            t (float): The time to propagate to in seconds.
            t_start (float, optional): The time to start the propagation 
                from in seconds.
            use_km (bool, optional): Whether to return the results in km. 
                Defaults to False.
            coord_frame (str, optional): The coordinate frame to return the 
                position in. Defaults to "ECI".
            max_step (int, optional): The step size to use when propergating

        Returns:
            tuple[np.ndarray, np.ndarray, np.ndarray]: The position, velocity
                and time steps arrays.
        """
        e = self.eccentricity
        theta = math.radians(self.true_anomaly)
        mu = self.body.gravitational_parameter
        a = self.semi_major_axis
        
        # Get initial mean anomaly
        E: float = math.acos( (e + math.cos(theta)) / (1 + e * math.cos(theta)) )
        E = E if theta < math.pi else 2 * math.pi - E
        
        M_e = E - e * math.sin(E)
        init_M_e = M_e if M_e > 0 else 2 * math.pi + M_e
        n = math.sqrt(mu / a**3)
        
        samples = int((t - t_start) / max_step) + 1
        
        r = np.zeros((samples, 3))      # take transpose later
        v = np.zeros((samples, 3))      # easier to iterate over rows
        t_steps = np.arange(t_start, t, max_step)
        
        # init_M_e = (M_e + n * t_start) % (2 * math.pi)
        h = math.sqrt(mu * a * (1 - e**2))
        
        # Propagate in perifocal frame
        for i, t_i in enumerate(t_steps):
            M_e = (init_M_e + n * t_i) % (2 * math.pi)
            theta = Orbit.mean_to_true_anomaly(math.degrees(M_e), e, degrees=False)
            
            r_i_mag = h**2 / (mu * (1 + e * math.cos(theta)))
            v_i_coeff = mu / h
            
            r_i = r_i_mag * np.array([math.cos(theta), math.sin(theta), 0])
            v_i = v_i_coeff * np.array([-math.sin(theta), e + math.cos(theta), 0])
            
            r[i] = r_i
            v[i] = v_i
        
        r = r.T if not use_km else r.T / 1000
        v = v.T if not use_km else v.T / 1000
        
        if coord_frame == "ECI":
            i = self.inclination
            rt_asc = self.right_ascension
            arg_p = self.argument_periapsis
            
            r = ot.perifocal_to_ECI_matrix(i, rt_asc, arg_p) @ r
            v = ot.perifocal_to_ECI_matrix(i, rt_asc, arg_p) @ v
        
        # Cache results
        self._analytical = True
        self._frame = coord_frame
        self._prop_in_km = use_km
        self._propergate_period = [t_start, t]
        self._propergate_step = max_step
        
        self._orbit_r = r
        self._orbit_v = v
        self._orbit_t = t_steps
        
        return r, v, 0, t_steps
               
    @staticmethod
    def mean_to_true_anomaly(
        mean_anomaly: float,
        e: float,
        degrees: bool = True
    ) -> float:
        """Finds the true anomaly of a satellite in its orbit
        using a taylor series expansion.

        Args:
            mean_anomaly (float): The mean anomaly in degrees.
            e (float): The eccentricity of the orbit.

        Returns:
            float: The true anomaly
        """
        M_e: float = math.radians(mean_anomaly)
        E: float = (
            M_e 
            + (e - e**3 / 8 + e**5 / 192) * math.sin(M_e)
            + (e**2 / 2 - e**4 / 6 + e**6 / 48) * math.sin(2 * M_e)
            + ((3 / 8) * e**3 - (27 / 128) * e**5) * math.sin(3 * M_e)
            + (e**4 / 3 - (4 / 15) * e**6) * math.sin(4 * M_e)
            + (125 / 384) * e**5 * math.sin(5 * M_e)
            + (27 / 80) * e**6 * math.sin(6 * M_e)
        )
        
        theta: float = math.acos( (e - math.cos(E)) / (e * math.cos(E) - 1) )
        theta = theta if M_e < math.pi else 2 * math.pi - theta
        
        return math.degrees(theta) if degrees else theta