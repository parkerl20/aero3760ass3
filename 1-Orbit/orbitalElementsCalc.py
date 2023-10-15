# -------- Add spacesim to script path
import sys
import os

current_file_path = os.path.abspath(__file__)
parent_directory = os.path.dirname(current_file_path)
grandparent_directory = os.path.dirname(parent_directory)
sys.path.append(grandparent_directory)

from spacesim import constants as const
import numpy as np

def orbitalParamaters(lat, lon, swathe_width):
    # For one satellite
    # Coordinates of the centre of NSW
    lon = np.deg2rad(lon)
    lat = np.deg2rad(lat)
    swathe_width = np.deg2rad(swathe_width)
    i = lat - swathe_width/2 # Inclination of the satellite

    # A target wll pass directly over a target of ground station on the earth if:
    L_node = lon - np.arcsin(np.tan(lat)/np.tan(i)) 
    
    return 0


def minimumAltitude(swathe_width, alpha):
    '''
    minimumAltitude - Calculates the minimium altitude a satellite can be positioned at above ground while still 
                      covering the required swathe_width, based on the spatial resolution of the sensor
    Inputs: swathe_width (float) = Degrees of latitude that the satellite must be able to observe
            alpha (float) = Spatial resoltuion of camera

    Outputs: altitude (float) = Minimum altitude the satellite can be placed at
    '''
    altitude = (swathe_width * 10000/90)/alpha
    return altitude + const.R_EARTH


def repeatingGroundTrackAltitude(j, k, i):
    """Calculates the altitude required to have a satellite with a groundtrack that repeats

    Args:
        j (float): Number of orbits per day (should be around 14 or 15 for LEO)
        k (float): Number of days after which the orbit ground track repeats
        i (float): Inclination of the orbit

    Returns:
        altitude (float): Altitude of satellite
    """
    # Start with an initial guess of altitude H0, which ignores Earth oblateness
    H0 = const.MU_EARTH**(1/3) * ((2 * np.pi * j)/(const.D_EARTH * k))**(-2/3) - const.R_EARTH
    alt_0 = H0 + const.R_EARTH

    # Constant used for calculations 
    k_2 = 0.75 * (360/2 * np.pi) * const.J2_EARTH * (const.MU_EARTH**(1/2)) * (const.R_EARTH**2)
    
    # Rotation rate of Earth
    L_dot = 360     # Degrees per sidereal day

    # Rate of change of the ascending node
    omega_dot = -2 * k_2 * (alt_0**(-7/2)) * np.cos(np.deg2rad(i))

    # Rate of change of perigee
    w_dot = k_2 * (alt_0**(-7/2)) * (5 * (np.cos(np.deg2rad(i))**2) - 1)

    # Rate of change of mean anomaly
    M_dot = k_2 * (alt_0**(-7/2)) * (3 * (np.cos(np.deg2rad(i))**2) - 1)

    # Mean angular motion
    n = (j/k) * (L_dot - omega_dot) - (w_dot + M_dot)

    # Final altitude calculation for a repeating ground track
    H = (const.MU_EARTH**(1/3)) * ((n * np.pi)/(180 * const.D_EARTH))**(-2/3) - const.R_EARTH
    alt = H + const.R_EARTH

    return alt


def flyOverRightAscension(lat, lon, i):
    lat = np.deg2rad(lat)
    i = np.deg2rad(i)

    rt_asc1 = lon - np.rad2deg(np.arcsin((np.arctan2(lat, 1))/(np.arctan2(i, 1))))
    rt_asc2 = lon - (np.rad2deg(np.arcsin((np.arctan2(lat, 1))/(np.arctan2(i, 1)))) - 180.0)
    return rt_asc1