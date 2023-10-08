#--TLE parsing and orbital parameters----------------------------------------------

import math
import numpy as np
from scipy.integrate import solve_ivp


class Satellite:  
    """Satellite class that takes in all necessary variables from the given TLE
    """
    def __init__(self, _name, _year, _day, _day_frac, _i, _omega, _e, _w, _m_e, _n):
        # Name of satellite
        self.name = _name

        # Epoch 
        self.year = _year 
        self.day = _day
        self.day_frac = _day_frac

        # Orbital Elements
        self.i = _i  # Inclination
        self.omega = _omega  # Right ascension of the ascending node
        self.e = _e   # Eccentricity
        self.w = _w    # Argument of perigee
        self.m_e = _m_e   # Mean anomaly
        self.n = _n    # Mean motion
        self.mu = 398600    # Gravitational constant


def read_TLE(filename):
    """Parses TLE variables to the satellite class

    Args:
        filename (str): path to a text file

    Returns:
        class: satelite class with parsed variables
    """
    with open(filename, "r") as file:
        lines = file.readlines()
        name = lines[0].strip()
        
        # Line 1 data
        year = int(lines[1][18:20])
        day = int(lines[1][20:23])
        day_frac = float(lines[1][23:32])
        
        # Line 2 data
        i = float(lines[2][8:16])
        omega = float(lines[2][17:25])
        e = float("0." + lines[2][26:33].replace(" ", ""))  # Get rid of space
        w = float(lines[2][34:42])
        m_e = float(lines[2][43:51])
        n = float(lines[2][52:63])

    # Creates SAT Satellite Object
    SAT = Satellite(name, year, day, day_frac, i, omega, e, w, m_e, n)

    return SAT

    

def element_to_perifocal(SAT):
    """Transforms TLE data into the perifocal frame

    Args:
        SAT (class): Satellite class

    Returns:
        X_PERI: Array of the Perifocal coordinates at epoch
        SAT: Satellite class
    """
    SAT.E = np.deg2rad(newt(0, SAT)) # Eccentricity anomaly with newtons method
    SAT.theta = np.arccos((SAT.e - np.cos(SAT.E))/(SAT.e*np.cos(SAT.E) - 1))    # True anomaly
    SAT.period = 2*np.pi/SAT.n  # Period of the orbit
    SAT.a = (SAT.mu)**(1/3)/(SAT.n*2*np.pi/(24*3600))**(2/3)   # Semi-major axis [Km]
    SAT.h = (SAT.a * SAT.mu * (1-SAT.e **2))**(1/2)    # Angular momentum [Km^2/s]

    # Calculate the magnitude of the radius vector
    r_mag = (SAT.h**2/SAT.mu) * (1/(1 + SAT.e * np.cos(SAT.theta)))

    # Calculate the components of the radius vector
    r_peri = r_mag * np.array([np.cos(SAT.theta), np.sin(SAT.theta), 0])

    # Calculate the components of the velocity vector
    v_peri = (SAT.mu/SAT.h) * np.array([-np.sin(SAT.theta), (SAT.e + np.cos(SAT.theta)), 0])

    X_PERI = np.array([r_peri, v_peri])
    return X_PERI, SAT

def perifocal_to_ECI(X_PERI, SAT):
    """Converts from the perifocal frame into the ECI frame for plotting

    Args:
        X_PERI (array): Perifocal frame
        SAT (class): Satellite class

    Returns:
        X_ECI: Array of ECI coordinates at epoch
    """
    # Transpose the input vectors for easier calculation
    r = np.transpose(X_PERI[0])
    v = np.transpose(X_PERI[1])

    # Euler angles
    om = np.deg2rad(SAT.omega)  # Argument of perigee
    i = np.deg2rad(SAT.i)  # Inclination
    w = np.deg2rad(SAT.w)  # Argument of perigee

    # Define transformation matrix Q using the Euler angles
    l1 = [-np.sin(om)*np.cos(i)*np.sin(w) + np.cos(om)*np.cos(w), -np.sin(om)*np.cos(i)*np.cos(w) - np.cos(om)*np.sin(w), np.sin(om)*np.sin(i)]
    l2 = [np.cos(om)*np.cos(i)*np.sin(w) + np.sin(om)*np.cos(w), np.cos(om)*np.cos(i)*np.cos(w) - np.sin(om)*np.sin(w), -np.cos(om)*np.sin(i)]
    l3 = [np.sin(i)*np.sin(w), np.sin(i)*np.cos(w), np.cos(i)]
    Q = np.array([l1, l2, l3])  

    # Transform the position and velocity vectors to the ECI frame
    r_ECI = np.dot(Q, r)
    v_ECI = np.dot(Q, v)

    # Output
    X_ECI = [r_ECI, v_ECI]
    return X_ECI
    

# Helper functions
def solve_ivp_fun(t, X_ECI):
    """Helper function for solve_ivp, integrates the function according to the two-body equation

    Args:
        t (int): Time between 0 and 48 hours in seconds
        X_ECI (array): ECI coordinates of satellite

    Returns:
        dXdt: Array of function to integrate
    """
    mu = 398600 # Gravitational constant
    x = X_ECI[0:3]  # Position data
    xdot = X_ECI[3:6]   # Velocity data
    xddot = -mu/(np.linalg.norm(x)**3) * x  # Two-body equation
    dXdt = np.concatenate([xdot, xddot])    # Function to integrate
    return dXdt

def solve_ivp_fun_j2(t, X_ECI):
    """Helper function for solve_ivp, integrates the function according to the two-body equation
    accounting for the J2 pertubation

    Args:
        t (int): Time between 0 and 48 hours in seconds
        X_ECI (array): ECI coordinates of satellite
        SAT (class): Satellite class

    Returns:
        dXdt: Array of function to integrate
    """
    mu = 398600 # Gravitational constant
    x = X_ECI[0:3]  # Position data
    xdot = X_ECI[3:6]   # Velocity data

    # J2 additions
    j2 = 0.00108263 # J2 constant
    r_e = 6378.137 # Earth radius
    r = np.linalg.norm(x)

    # Two body equation with J2
    ddx_j2 = -mu/(r**3) * (1 - (3/2)*j2*(r_e/r)**2 * (5*(x[2] /r)**2 - 1)) * x[0]
    ddy_j2 = -mu/(r**3) * (1 - (3/2)*j2*(r_e/r)**2 * (5*(x[2] /r)**2 - 1)) * x[1]
    ddz_j2 = -mu/(r**3) * (1 - (3/2)*j2*(r_e/r)**2 * (5*(3*(x[2] /r)**2 - 1))) * x[2]

    # Output
    xddot = [ddx_j2, ddy_j2, ddz_j2]
    dXdt = np.concatenate([xdot, xddot])
    return dXdt

def newt(guess, SAT):
    """Newton's approximation method for calculating the eccentricity anomaly

    Args:
        guess (_type_): _description_
        SAT (_type_): _description_

    Returns:
        _type_: _description_
    """
    f0 = guess - (SAT.e * math.sin(guess)) - SAT.m_e
    f0_dash = 1 - (SAT.e * math.cos(guess))
    y = guess - f0/f0_dash
    err = abs(y - guess)
    guess = y
    if err > 1e-5:  # Tolerance for the function
        return newt(guess, SAT)
    return y

def print_satellite_elements(SAT):
    """Prints all satellite elements to the terminal

    Args:
        SAT (class): Satellite information
    """
    for name, value in vars(SAT).items():
        print(f"{name} = {value}")
