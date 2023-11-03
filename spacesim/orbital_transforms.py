"""Coordinate transformations between different frames of reference.
"""
from spacesim import rotations as R
from spacesim import constants as const
from spacesim import time_util as tutil

import numpy as np
import datetime as dt
import math


def elements_to_perifocal(a: float, e: float, theta: float) -> tuple[np.ndarray, np.ndarray]:
    """Converts the classical orbital elements of an orbit into the perifocal frame.

    Args:
        a (float): The semimajor axis of the orbit.
        e (float): The eccentricity of the orbit.
        theta (float): The true anomaly in degrees.

    Returns:
        tuple[np.ndarray, np.ndarray]: The position and velocity vectors in 
        the perifocal frame as (3 x 1) vectors.
    """
    theta_rad = np.radians(theta)
    h = np.sqrt(const.MU_EARTH * a * (1 - e**2))                        # Specific angular momentum
    r = (h**2 / const.MU_EARTH) * (1 / (1 + e * np.cos(theta_rad)))     # Radius
    v_p = const.MU_EARTH / h                                            # Perifocal velocity

    r_perifocal = r * np.array([[np.cos(theta_rad), np.sin(theta_rad), 0]], dtype=np.float64)
    v_perifocal = v_p * np.array([[-np.sin(theta_rad), e + np.cos(theta_rad), 0]], dtype=np.float64)
    
    return r_perifocal.T, v_perifocal.T

def perifocal_to_ECI_matrix(i: float, rt_asc: float, arg_p: float) -> np.ndarray:
    """Creates a matrix to convert a vector in the perifocal frame to the ECI frame.

    Args:
        i (float): Inclination of the orbit in degrees.
        rt_asc (float): Right ascension of the ascending node in degrees.
        arg_p (float): Argument of perigee in degrees.

    Returns:
        np.ndarray: A (3 x 3) matrix to convert a vector in the 
            perifocal frame to the ECI frame.
    """
    # *_r denotes radians
    i_r = np.radians(i)
    rt_asc_r = np.radians(rt_asc)
    arg_p_r = np.radians(arg_p)
    
    M11 = -np.sin(rt_asc_r) * np.cos(i_r) * np.sin(arg_p_r) + np.cos(rt_asc_r) * np.cos(arg_p_r)
    M12 = -np.sin(rt_asc_r) * np.cos(i_r) * np.cos(arg_p_r) - np.cos(rt_asc_r) * np.sin(arg_p_r)
    M13 = np.sin(rt_asc_r) * np.sin(i_r)
    
    M21 = np.cos(rt_asc_r) * np.cos(i_r) * np.sin(arg_p_r) + np.sin(rt_asc_r) * np.cos(arg_p_r)
    M22 = np.cos(rt_asc_r) * np.cos(i_r) * np.cos(arg_p_r) - np.sin(rt_asc_r) * np.sin(arg_p_r)
    M23 = -np.cos(rt_asc_r) * np.sin(i_r)
    
    M31 = np.sin(i_r) * np.sin(arg_p_r)
    M32 = np.sin(i_r) * np.cos(arg_p_r)
    M33 = np.cos(i_r)
    
    return np.array(
        [
            [M11, M12, M13],
            [M21, M22, M23],
            [M31, M32, M33]
        ], 
        dtype=np.float64
    )
    
def perifocal_to_ECI(x_perifocal: np.ndarray, i: float, rt_asc: float, arg_p: float) -> np.ndarray:
    """Converts a vector in the perifocal frame to the ECI frame.

    Args:
        x_perifocal (np.ndarray): A (3 x 1) vector in the perifocal frame.
        i (float): Inclination of the orbit in degrees.
        rt_asc (float): Right ascension of the ascending node in degrees.
        arg_p (float): Argument of perigee in degrees.

    Returns:
        np.ndarray: A (3 x 1) vector in the ECI frame.
    """    
    return perifocal_to_ECI_matrix(i, rt_asc, arg_p) @ x_perifocal

def elements_to_ECI( a: float, e: float, i: float, rt_asc: float,
                    arg_p: float, theta: float) -> tuple[np.ndarray, np.ndarray]:
    """Converts the classical orbital elements of an orbit into the ECI frame.

    Args:
        theta (float): The true anomaly in degrees.
        a (float): The semimajor axis of the orbit.
        e (float): The eccentricity of the orbit.
        i (float): Inclination of the orbit in degrees.
        rt_asc (float): Right ascension of the ascending node in degrees.
        arg_p (float): Argument of perigee in degrees.

    Returns:
        tuple[np.ndarray, np.ndarray]: The position and velocity vectors in
        the ECI frame as (3 x 1) vectors.
    """
    r_p, v_p = elements_to_perifocal(a, e, theta)
    r_ECI = perifocal_to_ECI(r_p, i, rt_asc, arg_p)
    v_ECI = perifocal_to_ECI(v_p, i, rt_asc, arg_p)
    
    return r_ECI, v_ECI

def ECI_to_elements(r: np.ndarray, v: np.ndarray, mu: float) -> np.ndarray:
    """Finds the orbital elements of an orbit given the position and velocity
    vectors in the ECI frame.

    Args:
        r (np.ndarray): The position vector in the ECI frame.
        v (np.ndarray): The velocity vector in the ECI frame.
        mu (float): The gravitational parameter of the central body.

    Returns:
        np.ndarray: The orbital elements as a (6 x 1) vector.
        In the following order:
            semimajor axis,\n
            eccentricity,\n
            inclination,\n
            right ascension of the ascending node,\n
            argument of perigee,\n
            true anomaly.
    """
    
    r = r.flatten()
    v = v.flatten()
    
    r_mag = math.sqrt(r.dot(r))
    r_hat = r / r_mag
    
    h = np.cross(r, v)                          # Specific angular momentum
    h_mag = math.sqrt(h.dot(h))
    
    e_vec = np.cross(v, h) / mu - r_hat
    e = np.linalg.norm(e_vec)                   # eccentricity
    
    k_hat = np.array([0,0,1], dtype=np.float64)
    N = np.cross(k_hat, h)                      # Node line
    N_mag = np.linalg.norm(N)
    
    a = h_mag**2 / (mu * (1 - e**2))            # semimajor axis
    
    RAAN = math.acos(N[0] / N_mag)              # Right accension
    RAAN = RAAN if N[1] >= 0 else 2 * math.pi - RAAN
    
    omega = math.acos(N.dot(e_vec) / (N_mag * e))         # Argument of perigee
    omega = omega if e_vec[2] >= 0 else 2 * math.pi - omega
    
    inclination = math.acos(h[2] / h_mag)                 # Inclination
    
    e_hat = e_vec / e
    theta = math.acos(e_hat.dot(r_hat))                   # True anomaly
    theta = theta if r.dot(v) >= 0 else 2 * math.pi - theta
    
    return np.array([[
        a,
        e,
        np.degrees(inclination),
        np.degrees(RAAN),
        np.degrees(omega),
        np.degrees(theta)
    ]], dtype=np.float64).T
    
def ECEF_to_ECI(ecef: np.ndarray, t: float, epoch: dt.datetime = None) -> np.ndarray:
    
    if epoch:
        ERA = tutil.gmst(epoch) * const.ROT_V_EARTH
    else:
        ERA = 0
        
    theta = const.ROT_V_EARTH * t + ERA
    rot_mat = R.rot3_z(theta)
    return rot_mat @ ecef        
    
def ECI_to_ECEF(eci: np.ndarray, t: float, epoch: dt.datetime = None) -> np.ndarray:
    """Converts a vector in the ECI frame to the ECEF frame.

    Args:
        eci (np.ndarray): A vector in the ECI frame.
        t (float): The time in seconds since the epoch.
        epoch (datetime): The epoch of the start of the orbit.

    Returns:
        np.ndarray: A vector in the ECEF frame.
    """
    eci = eci.reshape((3,1))    # Ensure column vector
    
    if epoch:
        ERA = tutil.gmst(epoch) * const.ROT_V_EARTH
    else:
        ERA = 0
        
    theta = const.ROT_V_EARTH * t + ERA
    rot_mat = R.rot3_z(theta).T
    return rot_mat @ eci

def polar_to_cartesian(polar: np.ndarray) -> np.ndarray:
    """Converts a vector in the polar coordinate system to the cartesian coordinate system.
    
    The unit of the calculated cartesian vector is the same as the unit of the range.

    Args:
        polar (np.ndarray): A (3 x 1) vector containing elevation (deg), 
        azimuth (deg), and range.

    Returns:
        np.ndarray: A (3 x 1) vector in the cartesian coordinate system.
    """
    elev, azim, rng = polar.flatten()
    az_rad = np.radians(azim)
    ele_rad = np.radians(elev)

    x = rng * np.cos(az_rad) * np.cos(ele_rad)
    y = rng * np.sin(az_rad) * np.cos(ele_rad)
    z = rng * np.sin(ele_rad)
    
    cart = np.array([[x, y, z]])
    return cart.T

def cartesian_to_polar(cart: np.ndarray) -> np.ndarray:
    """Converts a vector in the cartesian coordinate system to the polar coordinate system.
    
    The unit of range is the same as the unit of the cartesian vector.

    Args:
        cart (np.ndarray): A (3 x 1) vector in the cartesian coordinate system.

    Returns:
        np.ndarray: A (3 x 1) vector containing elevation (deg), azimuth (deg), and range.
    """
    x,y,z = cart.flatten()
    rng = np.sqrt(x**2 + y**2 + z**2)
    azimuth = np.arctan2(y, x)
    elevation = np.arcsin(z / rng)
    
    return np.array([[np.degrees(elevation), np.degrees(azimuth), rng]], dtype=np.float64).T

def geocentric_LLH_to_ECEF(llh: np.ndarray, *, use_km: bool = False) -> np.ndarray:
    """Converts the geocentric LLH of a point to the ECEF frame.

    Args:
        llh (np.ndarray): A (3 x 1) vector containing latitude (deg), longitude (deg), and height.
        use_km (bool, optional): Whether to use km as the unit of height. Defaults to False.

    Returns:
        np.ndarray: A (3 x 1) position vector in the ECEF frame.
    """
    lat, lng, height = llh.flatten()    
    lat_rad = np.radians(lat)
    lng_rad = np.radians(lng)
    
    R = const.R_EARTH / 1000 if use_km else const.R_EARTH
    
    x = (R + height) * np.cos(lng_rad) * np.cos(lat_rad)
    y = (R + height) * np.sin(lng_rad) * np.cos(lat_rad)
    z = (R + height) * np.sin(lat_rad)
    
    return np.array([[x, y, z]], dtype=np.float64).T
    
def ECEF_to_geocentric_LLH(ecef: np.ndarray, *, use_km: bool = False) -> np.ndarray:
    """Converts a position vector in the ECEF frame to geocentric LLH coordinates.

    Args:
        ecef (np.ndarray): A (3 x 1) position vector in the ECEF frame.
        use_km (bool, optional): Whether to use km as the unit of height. Defaults to False.

    Returns:
        np.ndarray: A (3 x 1) vector containing latitude (deg), longitude (deg), and height.
    """
    x, y, z = ecef.flatten()
    R = const.R_EARTH / 1000 if use_km else const.R_EARTH
    
    radial = np.sqrt(x**2 + y**2 + z**2)
    lat_rad = np.arcsin(z / radial)
    lng_rad = np.arctan2(y, x)
    h = radial - R
    
    return np.array([[np.degrees(lat_rad), np.degrees(lng_rad), h]], dtype=np.float64).T
    
def LLH_to_ECEF(llh: np.ndarray, *, use_km: bool = False) -> np.ndarray:
    """Converts a vector in the geodedic LLH frame to the ECEF frame.
    
    Args:
        llh (np.ndarray): A (3 x 1) vector containing latitude (deg), longitude (deg) 
        and height (m).
        use_km (bool, optional): Whether to use km as the unit of altitude. 
        Defaults to False.

    Returns:
        np.ndarray: A (3 x 1) vector in the ECEF frame.
    """
    lat, lng, height = llh.flatten()
    lat_rad = np.radians(lat)
    lng_rad = np.radians(lng)
    
    R_N = geodetic_radius(lat)
    
    if use_km:
        height /= 1000
        R_N /= 1000
        
    e_2 = const.e2_EARTH
    
    x = (R_N + height) * np.cos(lat_rad) * np.cos(lng_rad)
    y = (R_N + height) * np.cos(lat_rad) * np.sin(lng_rad)
    z = (R_N * (1 - e_2) + height) * np.sin(lat_rad)
    
    return np.array([[x, y, z]], dtype=np.float64).T

def ECEF_to_LLH(ecef: np.ndarray, *, use_km: bool = False) -> np.ndarray:
    """Applies Ferrari's solution to convert a coordinate in the ECEF frame to the
    geodetic LLH frame.
    
    Source: https://en.wikipedia.org/wiki/Geographic_coordinate_conversion#The_application_of_Ferrari's_solution

    Args:
        eci (np.ndarray): A (3 x 1) vector in the ECI frame.
        use_km (bool, optional): Whether to use km as the unit of altitude. Defaults to False.

    Returns:
        np.ndarray: A (3 x 1) vector containing latitude (deg), longitude (deg) and height.
    """
    ecef *= 1000 if use_km else 1   # Convert to metres temporarily if required
    
    x,y,z = ecef.flatten()
    a = const.R_EQAT_EARTH
    b = const.R_POLE_EARTH
    e_2 = (a**2 - b**2) / a**2
    e_2_prime = (a**2 - b**2) / b**2
    p = np.sqrt(x**2 + y**2)
    F = 54 * b**2 * z**2
    G = p**2 + (1 - e_2) * z**2 - e_2 * (a**2 - b**2)
    c = (e_2**2 * F * p**2) / (G**3)
    s = (1 + c + np.sqrt(c**2 + 2 * c))**(1/3)
    k = s + 1 / s + 1
    P = F / (3 * k**2 * G**2)
    Q = np.sqrt(1 + 2 * e_2**2 * P)
    r_0 = -(P * e_2 * p) / (1 + Q) + np.sqrt(0.5 * a**2 * (1 + 1 / Q) - P * (1 - e_2) * z**2 / (Q * (1 + Q)) - 0.5 * P * p**2)
    U = np.sqrt((p - e_2 * r_0)**2 + z**2)
    V = np.sqrt((p - e_2 * r_0)**2 + (1 - e_2) * z**2)
    z_0 = (b**2 * z) / (a * V)
    
    # Final values
    h = U * (1 - b**2 / (a * V))
    latitude = np.arctan((z + e_2_prime * z_0) / p)
    longitude = np.arctan2(y, x)
    
    h /= 1000 if use_km else 1
    
    return np.array([[np.degrees(latitude), np.degrees(longitude), h]], dtype=np.float64).T

def iter_ECEF_to_LLH(ecef: np.ndarray, *, use_km: bool = False, 
                     iterations: int = 5) -> np.ndarray:
    """Applies an iterative method to convert a coordinate in the ECEF frame to the
    geodetic LLH frame.

    Args:
        ecef (np.ndarray): A (3 x 1) vector in the ECEF frame.
        use_km (bool, optional): Whether distance is measured in km.

    Returns:
        np.ndarray: A (3 x 1) vector containing latitude (deg), longitude (deg) and height.
    """
    x,y,z = ecef.flatten()
    R_xy = np.sqrt(x**2 + y**2)
    e_2 = const.e2_EARTH
    
    # initial values
    phi_deg,_,_ = ECEF_to_geocentric_LLH(ecef, use_km=use_km).flatten()
    phi = np.radians(phi_deg)
    h = (R_xy / np.cos(phi)) - geodetic_radius(phi_deg)
    
    # Perform newton's method
    for _ in range(iterations):
        N = geodetic_radius(phi_deg)
        phi = np.arctan((z / R_xy) / (1 - e_2 * N / (N + h)))
        h = (R_xy / np.cos(phi)) - N
        phi_deg = np.degrees(phi)
    
    lng = np.arctan2(y, x)

    return np.array([[phi_deg, np.degrees(lng), h]], dtype=np.float64).T
    
def geodetic_radius(lat_geodetic: float) -> float:
    """Finds the radius of earth at a given latitude, accounting
    for Earth's oblateness.

    Args:
        lat_geodetic (float): The geodetic latitude in degrees.

    Returns:
        float: The radius of earth in metres.
    """
    lat_rad = np.radians(lat_geodetic)
    a = const.R_EQAT_EARTH
    e_2 = const.e2_EARTH
    return a / np.sqrt(1 - e_2 * (np.sin(lat_rad) ** 2))

def ECEF_to_ENU_matrix(lat: float, lng: float) -> np.ndarray:
    """Provides a rotational matrix to transform a vector in the
    ECEF frame to the ENU frame.
    
    The transpose of the matrix can be used to convert from the ENU
    frame to the ECEF frame.
    
    Source: https://en.wikipedia.org/wiki/Geographic_coordinate_conversion#Geodetic_to/from_ENU_coordinates

    Args:
        lat (float): The latitude in degrees.
        lng (float): The longitude in degrees

    Returns:
        np.ndarray: A matrix to convert a vector in the ECEF frame
        to the ENU frame.   
    """
    lat_rad = np.radians(lat)
    lng_rad = np.radians(lng)
    
    M11 = -np.sin(lng_rad)
    M12 = np.cos(lng_rad)
    M13 = 0
    
    M21 = -np.sin(lat_rad) * np.cos(lng_rad)
    M22 = -np.sin(lat_rad) * np.sin(lng_rad)
    M23 = np.cos(lat_rad)
    
    M31 = np.cos(lat_rad) * np.cos(lng_rad)
    M32 = np.cos(lat_rad) * np.sin(lng_rad)
    M33 = np.sin(lat_rad)
    
    return np.array([
        [M11, M12, M13],
        [M21, M22, M23],
        [M31, M32, M33]
    ])

def ECEF_to_ENU(
    ground_pos: tuple[float, float, float],
    satellite_pos: np.ndarray,
    *,
    use_km: bool = False,
    geocentric: bool = False
) -> np.ndarray:
    """Converts a satellite's position in the ECEF frame to the East-North-Up (ENU) frame.
    
    Source: https://en.wikipedia.org/wiki/Geographic_coordinate_conversion#Geodetic_to/from_ENU_coordinates

    Args:
        ground_pos (tuple[float, float, float]): The ground station's latitude (deg), longitude (deg)
        and height.
        satellite_pos (np.ndarray): The satellite's position in the ECEF frame as a (3 x 1) vector.
        use_km (bool, optional): Whether to use km as the unit of altitude and range. Defaults to False.
        geocentric (bool, optional): Whether a geocentric or geodetic model of the earth should be used.
        Defaults to False.
        
    Returns:
        np.ndarray: The satellite's position in the ENU frame as a (3 x 1) vector.
    """
    ground_lat, ground_lng, _ = ground_pos
    ground_llh = np.array([[*ground_pos]], dtype=np.float64).T
    ground_ecef = LLH_to_ECEF(ground_llh, use_km=use_km)
    
    if geocentric:
        ground_ecef = geocentric_LLH_to_ECEF(ground_llh, use_km=use_km)
    
    satellite_rel_ecef = satellite_pos - ground_ecef
    return ECEF_to_ENU_matrix(ground_lat, ground_lng) @ satellite_rel_ecef
    
def ENU_to_ECEF(
    ground_pos: tuple[float, float, float],
    enu: np.ndarray,
    *,
    use_km: bool = False,
    geocentric: bool = False
) -> np.ndarray:
    """Converts a satellite's position in the East-North-Up (ENU) frame to the ECEF frame.
    
    Source: https://en.wikipedia.org/wiki/Geographic_coordinate_conversion#Geodetic_to/from_ENU_coordinates

    Args:
        ground_pos (tuple[float, float]): The ground station's latitude (deg), longitude (deg) and height.
        enu (np.ndarray): A (3 x 1) vector in the ENU frame.
        use_km (bool, optional): Whether to use km as the unit of altitude and range. Defaults to False.
        geocentric (bool, optional): Whether to use the geocentric or geodetic frame. Defaults to False.
        
    Returns:
        np.ndarray: The ECEF position of the satellite.
    """
    ground_lat, ground_lng, _ = ground_pos
    ground_llh = np.array([[*ground_pos]], dtype=np.float64).T
    ground_ecef = LLH_to_ECEF(ground_llh, use_km=use_km)
    
    if geocentric:
        ground_ecef = geocentric_LLH_to_ECEF(ground_llh, use_km=use_km)
    
    return ground_ecef + ECEF_to_ENU_matrix(ground_lat, ground_lng).T @ enu

def ENU_to_ECI(
    ground_pos: tuple[float, float, float],
    enu: np.ndarray,
    epoch: dt.datetime,
    *,
    use_km: bool = False,
    geocentric: bool = False
) -> np.ndarray:
    """Converts a satellite's position in the East-North-Up (ENU) frame to the ECI frame.

    Args:
        ground_pos (tuple[float, float, float]): The ground station's latitude (deg), 
        longitude (deg) and height.
        enu (np.ndarray): A (3 x 1) vector in the ENU frame.
        epoch (dt.datetime): The time of the observation.
        use_km (bool, optional): Whether to use km as the unit of altitude and range. 
        Defaults to False.
        geocentric (bool, optional): Whether to use the geocentric or geodetic model.

    Returns:
        np.ndarray: A (3 x 1) vector in the ECI frame.
    """
    ecef = ENU_to_ECEF(ground_pos, enu, use_km=use_km, geocentric=geocentric)
    return ECEF_to_ECI(ecef, 0, epoch)

def ECI_to_ENU(
    ground_pos: tuple[float, float, float],
    eci: np.ndarray,
    epoch: dt.datetime,
    *,
    use_km: bool = False,
    geocentric: bool = False
) -> np.ndarray:
    ecef = ECI_to_ECEF(eci, 0, epoch)
    return ECEF_to_ENU(ground_pos, ecef, use_km=use_km, geocentric=geocentric)

def ECI_to_mapping_error(
    r_eci: np.ndarray,
    r_true: np.ndarray,
    v_true: np.ndarray
) -> tuple[float, float, float]:
    """Finds the In-track, Cross-track and radial position errors
    of a satellite. Used for determining mapping error.

    Args:
        r_eci (np.ndarray): Estimated position in the ECI frame.
        r_true (np.ndarray): True position in the ECI frame.
        v_true (np.ndarray): True velocity in the ECI frame.

    Returns:
        tuple[float, float, float]: In-track error, cross-track error, radial error.
    """
    r_eci =  r_eci.flatten()
    r_true = r_true.flatten()
    v_true = v_true.flatten()
    
    r_error = r_eci - r_true
    
    h = np.cross(r_true, v_true)    # h is normal to the orbital plane
    
    v_unit = v_true / np.linalg.norm(v_true)
    h_unit = h / np.linalg.norm(h)
    r_unit = r_true / np.linalg.norm(r_true)
    
    # vector projection - In track error
    I_vec = (np.dot(r_error, v_true) / np.linalg.norm(v_true) ** 2) * v_true
    delta_I = np.linalg.norm(I_vec)
    
    # determine direction of in-track error
    if np.dot(I_vec, v_unit) < 0:
        delta_I *= -1
    
    # vector projection - Radial error
    proj_r = (np.dot(r_eci, r_true) / np.linalg.norm(r_true) ** 2) * r_true
    R_vec = proj_r - r_eci   
    delta_R = np.linalg.norm(R_vec)
    
    if np.dot(R_vec, r_unit) < 0:
        delta_R *= -1
    

    # vector projection - Cross track error
    C_vec = (np.dot(r_eci, h) / np.linalg.norm(h) ** 2) * h   # cross-track error
    delta_C = np.linalg.norm(C_vec)
    
    if np.dot(C_vec, h_unit) < 0:
        delta_C *= -1
    
    # print(delta_C)
    
    return delta_I, delta_C, delta_R

def ECI_to_azimuth_error(
    euler_truths: np.ndarray,
    euler_estimate: np.ndarray
) -> tuple[float, float]:
    """Finds the Azimuth and Elevation errors of a 
    satellite. Used for determining mapping error.

    Args:
        euler_truths (np.ndarray): Propagated euler truths
        euler_estimate (np.ndarray): Estimated euler from NLLS
        

    Returns:
        tuple[float, float, float]: Azimuth error, Elevation error.
    """

    # # Errors from roll/pitch/yaw over time with Weightings
    # errors = euler_truths - euler_estimate 

    # # Pitch = elevation, yaw = azimuth

    # delta_azimuth = errors[2]
    # delta_elevation = errors[1]

    # return delta_azimuth, delta_elevation

    return 0, 0

