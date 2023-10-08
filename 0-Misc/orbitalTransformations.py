#--Function of various coordiante transformations-------------------------------------------------

import numpy as np 
import math
import datetime

# Cartesian to Polar coordinates
def car_to_pol(X_CAR):
    x, y, z = X_CAR

    # Transformation
    r = np.sqrt(x**2 + y**2 + z**2)
    theta = np.arctan2(y, x)  
    phi = np.arcsin(z/r)  

    # Convert to degrees
    theta = np.rad2deg(theta)
    phi = np.rad2deg(phi)

    # Output
    X_POL = np.array([r, theta, phi])
    return X_POL

# Cartesion to Polar Coordiantes with ENU frame
def car_to_pol_enu(X_CAR):
    y, x, z = X_CAR

    # Transformation
    r = np.sqrt(x**2 + y**2 + z**2)
    theta = np.arctan2(x, y)  
    #phi = np.arcsin(z/r) 
    phi = np.arctan2(z, r) 

    # Convert to degrees
    theta = np.rad2deg(theta)
    phi = np.rad2deg(phi)

    # Output
    X_POL = np.array([r, theta, phi])
    return X_POL


# Polar to Cartesian coordinates
def pol_to_car(X_POL):
    r, theta, phi = X_POL
    # if(theta < 0):
    #     theta = theta + 360
    theta = np.deg2rad(theta)
    phi = np.deg2rad(phi)

    # Transformation
    x = r*np.cos(theta)*np.cos(phi)
    y = r*np.sin(theta)*np.cos(phi)
    z = r*np.sin(phi)

    # Output
    X_CAR = np.array([x, y, z])
    return X_CAR

# Polar to Cartesian coordinates
def pol_to_car_enu(X_POL):
    r, theta, phi = X_POL
    # if(theta < 0):
    #     theta = theta + 360
    theta = np.deg2rad(theta)
    phi = np.deg2rad(phi)

    # Transformation
    x = r*np.cos(theta)*np.cos(phi)
    y = r*np.sin(theta)*np.cos(phi)
    z = r*np.sin(phi)

    # Output
    X_CAR = np.array([y, x, z])
    return X_CAR


# ECEF frame to ECI frame
def ecef_to_eci(X_ECEF, t):
    # Angular velocity of the Earth (rad/s)
    w_ie = 7.292116*(10**(-5))  

    # Transformation matrix
    l1 = np.array([np.cos(w_ie * t), -np.sin(w_ie * t), 0])
    l2 = np.array([np.sin(w_ie * t), np.cos(w_ie * t), 0])
    l3 = np.array([0, 0, 1])

    ECEF_to_ECI = np.array([l1, l2, l3])

    # Output
    X_ECI = np.dot(ECEF_to_ECI, X_ECEF)
    return X_ECI

# ECEF to LLH coordinates
# Constants
a = 6378137  # Equartorial radius (m)
def ecef_to_llh(X_ECEF):
    x, y, z = X_ECEF
    
    # Latitude
    lat = np.arctan2(z, np.sqrt(x**2 + y**2))
    lat = np.rad2deg(lat)

    # Longitude
    long = np.arctan2(y, x)
    long = np.rad2deg(long)

    # Height
    R = np.sqrt((x**2) + (y**2) + (z**2))
    height = R - a

    llh = np.array([lat, long, height])
    
    return llh



# ECI to ECEF Coordinates
# Initalise variables
t = 180 # Time since vernal equinox (s)
def eci_to_ecef(X_ECI, t):
    # t = -t
    # Angular velocity of the earth (rad/s)
    w_ie = 7.292116*(10**(-5)) 

    # Transformation matrix
    l1 = np.array([np.cos(w_ie * t), np.sin(w_ie * t), 0])
    l2 = np.array([-np.sin(w_ie * t), np.cos(w_ie * t), 0])
    l3 = np.array([0, 0, 1])

    ECI_to_ECEF = np.array([l1, l2, l3])

    # Output
    X_ECEF = np.dot(ECI_to_ECEF, X_ECI)
    return X_ECEF



# LGDV to ECEF coordinates
def lgdv_to_ecef(X_LGDV):
    lat, lon, h = X_LGDV
    h = h * 1000    # Accounts for m
    lat = np.deg2rad(lat)
    lon = np.deg2rad(lon)

    # Variables
    a = 6378137 # radius of the Earth
    e_2 = 6.69437999014*(10**(-3))
    r_n = a/(1-e_2*(np.sin(lat)**2))**(1/2)

    # Transformations
    x = (r_n + h) * np.cos(lon) * np.cos(lat)
    y = (r_n + h) * np.sin(lon) * np.cos(lat)
    z = ((r_n * (1 - e_2)) + h) * (np.sin(lat))

    # Output
    X_ECEF = np.array([x, y, z])

    return X_ECEF/1000


# LGCV to ECEF coordinates, usually for groundstations
def lgcv_to_ecef(LLH):
    """
    Converts a position in the Geocentric Latitude, Longitude, and Height coordinate system
    to a position vector in the Earth Center Earth Fixed (ECEF) coordinate system.
    
    Args:
    - LLH (list or numpy array): [latitude, longitude, height]
    
    Returns:
    - numpy array: ECEF coordinates [rx, ry, rz]
    """
    
    # Constants
    Re = 6378.137  # Radius of the earth in m
    
    # Extract latitude, longitude, and height
    latitude, longitude, height = LLH
    
    # Convert degrees to radians
    lam = np.deg2rad(latitude)
    phi = np.deg2rad(longitude)
    
    # Compute rx, ry, rz
    R = height + Re
    rx = R * np.cos(lam) * np.cos(phi)
    ry = R * np.cos(lam) * np.sin(phi)
    rz = R * np.sin(lam)
    
    # Final output
    X_ECEF = np.array([rx, ry, rz])
    
    return X_ECEF


# ECEF to ENU Coordinates
# Takes in a relative vector rather than ecef
def ecef_to_enu_v2(X_REL, X_GROUND):
    lat, lon, h = X_GROUND
    lat = np.deg2rad(lat)
    lon = np.deg2rad(lon)

    # Transformation
    l1 = np.array([-np.sin(lon), np.cos(lon), 0])
    l2 = np.array([-np.sin(lat) * np.cos(lon), -np.sin(lat) * np.sin(lon), np.cos(lat)])
    l3 = np.array([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)])

    # Transformation matrix
    X_REL_to_X_ENU = np.array([l1, l2, l3])

    # Output
    X_ENU = np.dot(X_REL_to_X_ENU, X_REL)
    return X_ENU


# ECEF to ENU Coordinates
def ecef_to_enu(X_ECEF, X_GROUND_ECEF, X_REL_to_X_ENU):
    # Relative vector
    X_REL = np.subtract(X_ECEF, X_GROUND_ECEF)

    # Dot product with transformation matrix
    X_ENU_ECEF = np.dot(X_REL_to_X_ENU, X_REL)
    return X_ENU_ECEF

def t_matrix_X_REL_to_X_ENU(X_GROUND):
    lat, lon, _ = X_GROUND
    lat = np.deg2rad(lat)
    lon = np.deg2rad(lon)

    l1 = np.array([-np.sin(lon), np.cos(lon), 0])
    l2 = np.array([-np.sin(lat) * np.cos(lon), -np.sin(lat) * np.sin(lon), np.cos(lat)])
    l3 = np.array([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)])
    
    return np.array([l1, l2, l3])

def t_matrix_X_ENU_to_X_REL(X_GROUND):
    return t_matrix_X_REL_to_X_ENU(X_GROUND).T  # Transposition of other matrix


def enu_to_ecef(X_ENU, X_GROUND_ECEF, X_REL_to_X_ECEF):
    # Dot product with transformation matrix
    X_REL = np.dot(X_REL_to_X_ECEF, X_ENU)

    # Add ground station coordinates
    X_ECEF = X_REL + X_GROUND_ECEF

    return X_ECEF

def enu_to_ecefv3(X_ENU, X_GROUND):
    # Variables
    x, y, z = X_ENU
    lat, lon, h = X_GROUND
    lat = np.deg2rad(lat)
    lon = np.deg2rad(lon)
    X_GROUND_ECEF = lgdv_to_ecef(X_GROUND)

    l1 = np.array([-np.sin(lon), -np.sin(lat) * np.cos(lon), np.cos(lat) * np.cos(lon)])
    l2 = np.array([np.cos(lon), -np.sin(lat) * np.sin(lon), np.cos(lat) * np.sin(lon)])
    l3 = np.array([0, np.cos(lat), np.sin(lat)])

    # Arrays needed
    Q = np.array([l1, l2, l3])
    Enu_vec = np.array([x, y, z])   # Accounts for east-north reversal
    X_ECEF = np.add(np.dot(Q, Enu_vec), X_GROUND_ECEF)
    return X_ECEF


# ECEF to Geodetic LLA
def ecef_to_geodetic_lla(X_ECEF):
    x, y, z = X_ECEF

    # WGS84 constants for semi-major and semi-minor axes (in Km)
    a = 6378.1370
    b = 6356.7523142
    
    # Calculate the first eccentricity
    e2 = 1 - (b**2 / a**2)
    
    # Calculate intermediate values
    r = math.sqrt(x**2 + y**2)
    E2 = a**2 - b**2
    F = 54 * b**2 * z**2
    G = r**2 + (1 - e2) * z**2 - e2 * E2
    c = (e2**2 * F * r**2) / (G**3)
    
    # Check for invalid values before proceeding
    sqrt_term = c**2 + 2 * c
    
    s = (1 + c + math.sqrt(sqrt_term))**(1/3)
    P = F / (3 * (s + 1/s + 1)**2 * G**2)
    Q = math.sqrt(1 + 2 * e2**2 * P)
    # r0 = (-P * e2 * r) / Q + math.sqrt(0.5 * a**2 * (1 + 1/Q) - P * (1 - e2) * z**2 / (Q * (1 + Q)) - 0.5 * P * r**2)
    r0 = (-P * e2 * r) / (1 + Q) + math.sqrt(0.5 * a**2 * (1 + 1/Q) - P * (1 - e2) * z**2 / (Q * (1 + Q)) - 0.5 * P * r**2)
    U = math.sqrt((r - e2 * r0)**2 + z**2)
    V = math.sqrt((r - e2 * r0)**2 + (1 - e2) * z**2)
    z0 = (b**2 * z) / (a * V)
    
    # Calculate the geodetic coordinates
    e_dash_squared = (a**2 - b**2) / b**2
    latitude = math.atan2((z + e_dash_squared * z0), r)
    longitude = math.atan2(y, x)
    altitude = U * (1 - b**2 / (a * V))
    
    # Return results in degrees and meters
    X_LGDV = np.array([math.degrees(latitude), math.degrees(longitude), altitude])
    return X_LGDV


def eci_pos_to_ecef_pos(XX_ECI, SAT, t):
    """ECI to ECEF position coordinates, accounts for time since vernal equinox

    Args:
        XX_ECI (Matrix): Coordinates of satellite in ECI
        SAT (Class): Orbital parameters
        t (float): Time since vernal equinox + epoch

    Returns:
        XX_ECEF: ECEF matrix of satellite coordinates
    """
    # Where XX_ECI is a matrix of ECI coordinates given by solve_ivp
    XX_ECI_rows = len(XX_ECI)
    shape = (XX_ECI_rows, 3)   # 3 columns for x, y and z axis
    XX_ECEF = np.empty(shape)   # Initialise matrix for ECEF position coordinates
    t_vernal = t_diff(SAT)
    for i in range(XX_ECI_rows):
        XX_ECEF[i] = eci_to_ecef(XX_ECI[i, 0:3], t_vernal + t[i])

    return XX_ECEF

def t_diff(SAT):
    JD2000 = 2451545    # J2000 Julian day
    J0 = JulianD_Midnight(SAT)
    days = J0 - JD2000
    seconds = SAT.day_frac * 24 * 60 * 60
    frac = days/36525
    w = 1.00273790935 + frac * (5.9 * 10**11)
    t = 24110.54841 + (8640184.812866 * frac) + (0.093104 * (frac**2)) - ((6.2 * 10**(-6)) * frac**3) + (w * seconds)
    t = t % 86400
    return t

def JulianD_Midnight(SAT):
    year = SAT.year + 2000
    date = datetime.date(year, 1, 1) + datetime.timedelta(days=SAT.day-1)
    # date = datetime.datetime(year, date_calc.month, date_calc.day, 0, 0, 0) # GMST date
    J0 = 367*year - int((7*(year + int((date.month + 9)/12)))/4) + int(275*date.month/9) + date.day + 1721013.5
    return J0


