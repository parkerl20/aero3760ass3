from spacesim import orbital_transforms as ot
from spacesim import orbit as ob
from spacesim import constants as const

from collections import namedtuple
from scipy import constants as spconst
import datetime as dt
import numpy as np
import math


class OrbitObservatory():
    """A ground station that can track orbits
    """
    
    def __init__(self, name: str, location: tuple[float,float,float]) -> None:
        """Initialise the ground station

        Args:
            name (str): The name of the ground station
            location (tuple[float,float,float]): The latitude (deg), longitude (deg)
            and elevation (m) of the ground station.

        """
        latitude, longitude, elevation = location
        
        self.name: str = name
        self.location = location
        
        llh = np.array([[latitude, longitude, elevation]]).T        
        self.ecef: np.ndarray = ot.LLH_to_ECEF(llh, use_km=False)
        
        #---- Laser parameters
        self.laser_wavelength = 0     # metres
        self.laser_hertz = 0               # Hz
        
        #---- Atmospheric parameters
        self.pressure = 0           # hPa
        self.temperature = 0        # degrees Celsius
        self.humidity = 0           # percent
        
        #---- Measurement noise
        self.range_std = 0.0042              # metres
        self.elevation_std = 0.000277778     # degrees
        self.azimuth_std = 0.000277778     # degrees
        return
    
    def observe_orbit(
        self,
        orbit: ob.Orbit,
        t: float,
        *,
        t_start: float = 0,
        use_km: bool = False,
        max_step: int = 10,
        analytical: bool = False,
        noise: bool = False,
    ) -> namedtuple:
        """Observes an orbit by propagating it from epoch

        Args:
            orbit (ob.Orbit): The orbit to observe
            t (float): The time to observe the orbit to
            use_km (bool, optional): Whether to use km or m. Defaults to False.
            t_start (float, optional): The time to start observing the orbit from. 
                Defaults to 0.
            max_step (int, optional): The maximum step size to use when propagating 
                the orbit. Defaults to 10. A lower max_step will result in a more accurate 
                propagation but will take longer.
        
        Returns:
            Each return field (excluding visiblity period) is a list of numpy arrays
            with each index corresponding to a visibility period.
        
            Returns a `namedtuple` with the following fields:
                elevation (`list[np.ndarray]`): The elevation of the satellite in degrees
                azimuth (`list[np.ndarray]`): The azimuth of the satellite in degrees
                range (`list[np.ndarray]`): The range of the satellite.
                time (`list[np.ndarray]`): The time since propagration epoch in seconds
                visibilty_period (`list[tuple[datetime]]`): A list of tuples of the visibility periods
        """
        # TODO: In future can observe from a starting datetime object
        r, _, t = orbit.propagate(
            t,
            use_km=False,
            t_start=t_start,
            max_step=max_step,
            analytical=analytical,
        )
        
        # Convert to numpy array later. Python lists more efficient at growing
        elevations, curr_elevations = [], []
        azimuths, curr_azimuths = [], []
        ranges, curr_ranges = [], []
        times, curr_times = [], []
        
        t = t.tolist()
        
        # Track visibility periods
        visibility: list[tuple] = []                # Visibility time stamps
        visible_flag = False
        visible_start: dt.datetime = None
        
        for i in range(len(t)):
            r_i = r[:,i].reshape(3,1)      # Convert to (3 x 1), indexing column returns 1D
            sat_ecef = ot.ECI_to_ECEF(r_i, t[i], orbit.epoch)
            sat_enu = ot.ECEF_to_ENU(self.location, sat_ecef, use_km=False)
            enu_polar = ot.cartesian_to_polar(sat_enu).flatten()
            
            elev, azim, rng = enu_polar
            
            if noise:
                turbulance_var = atmospheric_turbulence_noise(elev)
                
                rng_error = self.range_std + math.sqrt(turbulance_var)
                
                elev += np.random.normal(0, self.elevation_std)
                azim += np.random.normal(0, self.azimuth_std)
                rng += np.random.normal(0, rng_error)
                
            # check visbility flag, set accordingly
            if not visible_flag and elev > 0:
                # visibility period started
                visible_flag = True
                visible_start = orbit.epoch + dt.timedelta(seconds=t[i])
            elif visible_flag and elev < 0:
                # visibility period ended
                visible_flag = False
                visible_end = orbit.epoch + dt.timedelta(seconds=t[i])
                visibility.append((visible_start, visible_end))
                
                elevations.append(np.array(curr_elevations))
                azimuths.append(np.array(curr_azimuths))
                ranges.append(np.array(curr_ranges))
                times.append(np.array(curr_times))
                
                curr_elevations, curr_azimuths, curr_ranges, curr_times = [], [], [], []
            
            if visible_flag:
                rng /= 1000 if use_km else 1
                
                curr_elevations.append(elev)
                curr_azimuths.append(azim)
                curr_ranges.append(rng)
                curr_times.append(t[i])
        
        # If ended while visible, add last visibility period
        if visible_flag:
            visible_end = orbit.epoch + dt.timedelta(seconds=t[-1])
            visibility.append((visible_start, visible_end))
            
            elevations.append(np.array(curr_elevations))
            azimuths.append(np.array(curr_azimuths))
            ranges.append(np.array(curr_ranges))
            times.append(np.array(curr_times))
        
        # namedtuple for return
        ObeserverReturn = namedtuple(
            "ObserverReturn",
            [
                "elevation",
                "azimuth", 
                "range", 
                "time",
                "visibility_period"
            ]
        )
         
        return ObeserverReturn(elevations, azimuths, ranges, times, visibility)
    
    def apply_SLR_measurement_noise(
        self,
        r_eci: np.ndarray,
        epoch: dt.datetime,
        *,
        use_km: bool = False,
        atmospheric_noise: bool = True
    ) -> tuple[np.ndarray, tuple]:
        """Applies SLR measurement from a ground station to a preexisting
        ECI vector
        
        Source: pg 104 of lecture 5 thesis paper

        Args:
            r (np.ndarray): _description_
            epoch (dt.datetime): _description_
            range_std (int, optional): _description_. Defaults to 30.
            elevation_std (float, optional): _description_. Defaults to 0.015.
            azimuth_std (float, optional): _description_. Defaults to 0.015.
            use_km (bool, optional): _description_. Defaults to False.

        Returns:
            tuple[np.ndarray, tuple]: The noisy ECI vector and the noise applied
                in the form (elevation_std, azimuth_std, range_std)
        """
        r_eci = r_eci.reshape((3, 1))       # ensure correct vector shape
        r_enu = ot.ECI_to_ENU(self.location, r_eci, epoch, use_km=use_km)
        
        r_ground_measure = ot.cartesian_to_polar(r_enu)
        elev = r_ground_measure[0]
        
        #------- Apply noise
        if atmospheric_noise:
            range_std = self.range_std + atmospheric_turbulence_noise(elev)
        else:
            range_std = self.range_std

        range_std /= 1000 if use_km else 1
        
        r_ground_measure[0] += np.random.normal(0, self.elevation_std)
        r_ground_measure[1] += np.random.normal(0, self.azimuth_std)
        r_ground_measure[2] += np.random.normal(0, range_std)
        
        # Convert back to ECI
        r_enu_noise = ot.polar_to_cartesian(r_ground_measure)
        r_eci_noise = ot.ENU_to_ECI(self.location, r_enu_noise, epoch, use_km=use_km)
        
        return r_eci_noise, (self.elevation_std, self.azimuth_std, range_std)


def atmospheric_turbulence_noise(elevation: float) -> float:
    """Calculates the atmospheric turbulence noise for a given elevation

    Args:
        elevation (float): The elevation in degrees.

    Returns:
        float: The atmospheric turbulence noise in meters.
    """
    return (2.601e-6 * (90 - elevation)**3  + 1.38) * 1e-3
  
  
def normal_point_info(filename: str) -> namedtuple:
    """Reads and returns information from the provided
    normal point file.

    Args:
        filename (str): The file to read from

    Returns:
        namedtuple: A named tuple containing the following fields:
            station_name (str): The name of the ground station
            laser_hertz (float): The laser frequency in Hz
            laser_wavelength (float): The laser wavelength in micro-metres
            range_std (float): The range standard deviation in metres
            humidity (float): The relative humidity in percent
            temperature (float): The temperature in degrees Celsius
            pressure (float): The atmospheric pressure in hPa
            location (tuple[float,float,float]): The latitude, longitude and elevation
                of the ground station in degrees, degrees and metres respectively.
    """
    name = ""
    laser_hertz = 0
    laser_wavelength = 0
    np_range_std = 0
    humidity = 0
    temperature = 0
    pressure = 0
    location = (0, 0, 0)    # lat, long, elev
    
    with open(filename, 'r') as f:
        for line in f:
            if len(line) < 2:
                continue
            
            line = line.strip().split()
            prefix = line[0].lower()
            
            if prefix == "other":
                lng = float(line[1])
                lat = float(line[2])
                elev = float(line[3])
                location = (lat, lng, elev)
                name = line[4]
            elif prefix == "c0":
                laser_wavelength = float(line[2])
            elif prefix == "c1":
                laser_hertz = float(line[5])
            elif prefix == "11":
                np_time_std = float(line[7])
                np_range_std = np_time_std * spconst.c * 1e-12
            elif prefix == "20":
                pressure = float(line[2])
                temperature = float(line[3]) - const.C_ZERO_K
                humidity = float(line[4])
    
    NormalPoint = namedtuple(
        "NormalPoint",
        [
            "station_name",
            "laser_hertz",
            "laser_wavelength",
            "range_std",
            "humidity",
            "temperature",
            "pressure",
            "location"
        ]
    )
    
    return NormalPoint(
        name,
        laser_hertz,
        laser_wavelength,
        np_range_std,
        humidity,
        temperature,
        pressure,
        location
    )

# --------------- Atmospheric range error functions
def atmospheric_range_error(
    latitude: float,
    height: float,
    wavelength: float,
    elevation: float,
    pressure: float,
    temperature: float,
    humidity: float,
) -> float:
    """Calculates the atmospheric range error
    for an SLR system.

    Args:
        latitude (float): The latitude of the ground station in degrees
        height (float): The height of the ground station in metres
        wavelength (float): The wavelength of the laser in micro-metres
        elevation (float): The elevation of the SLR laser in degrees
        pressure (float): The atmospheric pressure in hPa
        temperature (float): The atmospheric temperature in degrees Celsius
        humidity (float): The relative humidity in percent
        use_km (bool, optional): Whether to use km or m. Defaults to False.

    Returns:
        float: The atmospheric range error in metres.
    """
    A = 0.002416579
    f_h = hydrostatic_dispersion(wavelength)
    f_s = f_hydro_denom(latitude, height)
    hydrostatic_error = A * f_h * pressure / f_s
    
    d_nh = non_hydrostatic_error(
        latitude,
        height,
        wavelength,
        temperature,
        humidity
    )
    
    d_zenith = hydrostatic_error + d_nh
    
    e_r = np.radians(elevation)
    # Black and Eisner mapping function
    MR = 1 / math.sqrt(1 - (math.cos(e_r) / (1.001))**2)
    
    return d_zenith * MR

def hydrostatic_dispersion(wavelength: float) -> float:
    """Calculates the hydrostatic dispersion.

    Args:
        wavelength (float): The wavelength of the laser in micro-metres

    Returns:
        float: The hydrostatic dispersion.
    """
    # Dispersion constants
    k_0 = 238.0185          # μ/m^2
    k_1 = 19990.975
    k_2 = 57.362
    k_3 = 579.551744
    CO2 = 0.99995995        # CO2 constant
    
    s = 1 / wavelength      # wave number
    
    # Implement equation
    A = k_1 * (k_0 + s**2) / (k_0 - s**2)**2
    B = k_3 * (k_2 + s**2) / (k_2 - s**2)**2
    
    return 1e-2 * (A + B) * CO2

def f_hydro_denom(latitude: float, height: float) -> float:
    H = ot.geodetic_radius(latitude) + height
    lat_r = np.radians(latitude)
    return 1 - 0.00266 * np.cos(2 * lat_r) - 0.00000028 * H

def non_hydrostatic_error(
    latitude: float,
    height: float,
    wavelength: float,
    temperature: float,
    humidity: float
) -> float:
    """_summary_

    Args:
        latitude (float): _description_
        height (float): _description_
        wavelength (float): _description_
        temperature (float): _description_
        humidity (float): _description_

    Returns:
        float: _description_
    """
    f_nh = non_hydrostatic_dispersion(wavelength)
    f_h = hydrostatic_dispersion(wavelength)
    e_s = water_vapour_pressure(temperature, humidity)
    f_s = f_hydro_denom(latitude, height)
    
    return 1e-4 * (5.316 * f_nh - 3.759 * f_h) * e_s / f_s

def water_vapour_pressure(
    temp: float,
    humidity: float
) -> float:
    """Calculates the water vapour pressure.

    Args:
        temp (float): The atmospheric temperature in degrees Celsius
        humidity (float): The relative humidity in percent

    Returns:
        float: The water vapour pressure.
    """
    return (6.1078 * ((7.5 * temp) / (temp + 237.3)) ** 10) * humidity / 100

def non_hydrostatic_dispersion(wavelength: float) -> float:
    """Calculates the non-hydrostatic dispersion.

    Args:
        wavelength (float): The wavelength of the laser in micro-metres

    Returns:
        float: The non-hydrostatic dispersion.
    """
    # Dispersion constants
    w_0 = 295.235        
    w_1 = 2.6422        # μm^2
    w_2 = -0.032380     # μm^4
    w_3 = 0.004028      # μm^6
    
    s = 1 / wavelength      # wave number
    return 0.003101 * (w_0 + (3 * w_1 * s**2) + (5 * w_2 * s**4) + (7 * w_3 * s**6))