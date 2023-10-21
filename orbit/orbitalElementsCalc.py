# -------- Add spacesim to script path
import sys
import os

current_file_path = os.path.abspath(__file__)
parent_directory = os.path.dirname(current_file_path)
grandparent_directory = os.path.dirname(parent_directory)
sys.path.append(grandparent_directory)

from spacesim import constants as const
import numpy as np
import matplotlib.pyplot as plt

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


def minimumAltitudeTest(swathe_width, alpha):
    '''
    minimumAltitude - Calculates the minimium altitude a satellite can be positioned at above ground while still 
                      covering the required swathe_width, based on the spatial resolution of the sensor
    Inputs: swathe_width (float) = Degrees of latitude that the satellite must be able to observe
            alpha (float) = Spatial resoltuion of camera

    Outputs: altitude (float) = Minimum altitude the satellite can be placed at
    '''
    altitude = (swathe_width * 10000/90)/alpha
    return altitude + const.R_EARTH


def minimumAltitude(swathe_width):
    '''
    minimumAltitude - Calculates the minimium altitude a satellite can be positioned at above ground while still 
                      covering the required swathe_width, based on the FOV of the camera
    Inputs: swathe_width (float) = Degrees of latitude that the satellite must be able to observe

    Outputs: min_altitude (float) = Minimum altitude the satellite can be placed at
    '''
    # Constants regarding camera sensor
    n = 4                        # Number of cameras
    R_E = const.R_EARTH                # Earth radius
    lambda_max = np.deg2rad(swathe_width)   # Swath Width
    beta = np.deg2rad(21.7)      # Camera FOV in degrees
    FOV = beta * n               # Total FOV

    # Minimum altitude required to meet coverage requirements
    min_altitude = (R_E * (np.sin(0.5*(lambda_max+FOV)) - np.sin(FOV/2)) / np.sin(FOV/2))
    print(min_altitude)
    return min_altitude + R_E


def repeatingGroundTrackAltitude(j, k, i, e):
    """Calculates the altitude required to have a satellite with a groundtrack that repeats

    Args:
        j (float): Number of orbits per day (should be around 14 or 15 for LEO)
        k (float): Number of days after which the orbit ground track repeats
        i (float): Inclination of the orbit
        e (float): Eccentricity of the orbit

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
    omega_dot = -2 * k_2 * (alt_0**(-7/2)) * np.cos(np.deg2rad(i)) * (1 - e**2)**(-2)

    # Rate of change of perigee
    w_dot = k_2 * (alt_0**(-7/2)) * (5 * (np.cos(np.deg2rad(i))**2) - 1) * (1 - e**2)**(-2)

    # Rate of change of mean anomaly
    M_dot = k_2 * (alt_0**(-7/2)) * (3 * (np.cos(np.deg2rad(i))**2) - 1) * (1 - e**2)**(-3/2)

    # Mean angular motion
    n = (j/k) * (L_dot - omega_dot) - (w_dot + M_dot)

    # Final altitude calculation for a repeating ground track
    H = (const.MU_EARTH**(1/3)) * ((n * np.pi)/(180 * const.D_EARTH))**(-2/3) - const.R_EARTH
    alt = H + const.R_EARTH

    return alt


def choosingOrbitsPerDay(a_min, i, j, k, show_results):
    # Array to record each altitude 
    alt = []

    # Call the repeatingGroundTrackAltitude function for each number of j orbits in a day
    for orbit_num in range(len(j)):
        alt.append((repeatingGroundTrackAltitude(j[orbit_num], k, i, e=0) - const.R_EARTH)/1000)

    # Select the orbit per day number that is nearest above coverage
    counter = 0
    while(1):
        if(alt[counter] < a_min - const.R_EARTH):
            break
        else:
            counter += 1
    j_final = j[counter-1]

    # Make lables readable 
    Set_Label_Sizes()

    # Create a figure and axis
    fig, ax = plt.subplots(figsize=(8, 6))

    # Add labels and title
    ax.set_xlabel('Number of Orbits in a Day (j)')
    ax.set_ylabel('Altitude (km)')
    ax.set_title('Altitude vs. Number of Orbits in a Day for repeatable ground track')
    
    # Add a grid for better readability
    ax.grid(True, linestyle='--', alpha=0.7)

    # Add horizontal lines to each point
    for altitude in alt:
        ax.hlines(altitude, xmin=min(j), xmax=max(j), colors='r', linewidth=2.5)

    # Minimum altitude for coverage
    a_min = (a_min - const.R_EARTH)
    ax.hlines(a_min, xmin=min(j), xmax=max(j), colors='g', linestyles='dotted', label='Altitude for coverage requirements', linewidth=2.5)

    # Plot the results of different altitudes
    ax.scatter(j, alt, label='Altitude (km)', color='b', marker='o', s=100)

    # Show legend
    ax.legend()

    # Display the plot
    if(show_results):
        plt.show()

    return j_final


def choosingRepeatEccentricity(a_min, i, j, k, show_results):
    # Array to keep track of altitudes
    perigee_radius = []

    # Array to keep track of eccentricities
    eccentricity = []

    # Converting minimum altitude for coverage for easier manipulation
    a_min = a_min - const.R_EARTH

    # Tolerance for eccentricity change
    TOL = 0.001

    # Select the eccentricity that gives the radius at perigee closest to minimum altitude for coverage
    counter, e = 0, 0
    while(1):
        # Semi-major 
        a = repeatingGroundTrackAltitude(j, k, i, e)

        # Semi-minor
        b = ((a**2 * (1 - e**2))**(1/2) - const.R_EARTH)/1000
        rp = (a * (1 - e) - const.R_EARTH)/1000

        # Append arrays for plotting
        perigee_radius.append(rp)
        eccentricity.append(e)

        # Coverage requirement
        if(rp < a_min):
            e -= TOL    # Gets the eccentricity value just before minimum coverage altitude
            break
        else:
            counter += 1
            e += TOL

    # # Make lables readable 
    Set_Label_Sizes()

    # # Create a figure and axis
    fig, ax = plt.subplots(figsize=(8, 6))

    # # Add labels and title
    ax.set_xlabel('Eccentricity of orbit (degrees)')
    ax.set_ylabel('Altitude (km)')
    ax.set_title('Altitude vs. Eccentricity of the orbit')
    
    # # Add a grid for better readability
    ax.grid(True, linestyle='--', alpha=0.7)

    # # Add horizontal lines to each point
    for altitude in perigee_radius:
        ax.hlines(altitude, xmin=min(eccentricity), xmax=max(eccentricity), colors='r', linewidth=2.5)

    # Minimum altitude for coverage
    ax.hlines(a_min, xmin=min(eccentricity), xmax=max(eccentricity), colors='g', label='Altitude for coverage requirements', linewidth=2.5)

    # Finalised eccentricity 
    ax.hlines(perigee_radius[-1], xmin=min(eccentricity), xmax=max(eccentricity), colors='y', label='Final perigee altitude', linewidth=2.5)

    # Plot the results of different altitudes
    ax.scatter(eccentricity, perigee_radius, label='Altitude (km)', color='b', marker='o', s=100)

    # # Show legend
    ax.legend()

    # # Display the plot
    if(show_results):
        plt.show()

    return e

def flyOverRightAscension(lat, lon, i):
    """Calculates the right ascension necessary to travel directly over a ground station point

    Args:
        lat (float): Latitude of fly over point
        lon (float): Longitude of fly over point
        i (degree): Inclination of orbit

    Returns:
        rt_asc: Right ascension of orbit
    """
    # Convert to radians
    lat = np.deg2rad(lat)
    i = np.deg2rad(i)

    # Calculate right ascension. There are 2 values for highest and lowest point inclination point
    rt_asc1 = lon - np.rad2deg(np.arcsin((np.arctan2(lat, 1))/(np.arctan2(i, 1))))
    rt_asc2 = lon - (np.rad2deg(np.arcsin((np.arctan2(lat, 1))/(np.arctan2(i, 1)))) - 180.0)

    # Selects rt_asc based on if the orbit is prograde or retrograde
    rt_asc = rt_asc1 if (i > 0) else rt_asc2

    return rt_asc


def Set_Label_Sizes():
    """Generate global paramaters for plotting
    """
    # Set the title, axis labels, axis ticks, and label sizes
    plt.rcParams['axes.titlesize'] = 28     # Increase title size
    plt.rcParams['axes.labelsize'] = 24     # Increase axis label size
    plt.rcParams['xtick.labelsize'] = 20    # Increase x-tick label size
    plt.rcParams['ytick.labelsize'] = 20    # Increase y-tick label size
    plt.rcParams['xtick.major.size'] = 15   # Increase x-tick size
    plt.rcParams['ytick.major.size'] = 15   # Increase y-tick size
    plt.rcParams['legend.fontsize'] = 18    # Increase the legend font size
