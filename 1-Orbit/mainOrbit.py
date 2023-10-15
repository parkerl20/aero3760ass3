import numpy as np
from orbitalElementsCalc import orbitalParamaters
from orbitalElementsCalc import minimumAltitude
from simulateOrbit import simulateOrbit
from orbitalElementsCalc import repeatingGroundTrackAltitude
# from spacesim import constants as =

def mainOrbit():
   # For example of multiple satellite plotting, run demo.py

 # General pseudocode for completeing the orbital design processes
   """"
   1. Start with a circular orbit of altitude 160km, inclination 0 degrees
   2. Analytically calculate the field of view (FOV) as well as frequency and
      duration around the region of interest: 146.9211, -31.2532 (lat/lon of NSW centre)
   3. Record measurements
   4. Repeat steps 2 and 3 for every altitude and inclination combination, altitude 160-1000km
      and inclination 0-90 and 0-(-90) degrees
   5. Repeat stes 2-4 for both a 1 satellite and 2 satellite orbit configuration (possibly also
      changing eccentricity for the single satellite)
   6. Compare results to satellite requirements and liase with Engineers in other sections
   7. Determine finalised orbital paramaters
   8. Propogate satellite orbit in the ECI frame
   9. Determine thrusters needed to counteract the J2 effect
   10. Send final ECI matrix to other sections as a "Truth"
   """
   # For one satellite
   # Coordinates of the centre of NSW
   lon = 148.5
   lat = -33

   # Other paramaters
   swathe_width = 10.0  # Swathe width for one satellite to cover NSW
   alpha = 2.4        # Spacial resolution of the chosen camera

   # Minimum altitude of satelltie to meet swathe_width requirement
   a_min = minimumAltitude(swathe_width, alpha)

   # Altitude for a repeating ground track orbit based on j: orbital periods and k:number of days
   a = repeatingGroundTrackAltitude(j=15, k=1, i=lat)

   # Calculate orbital paramaters based on a latitude, longitude and swathe_width
   orbitalParamaters(lat, lon, swathe_width)
   print(f"altitude: {a}")
   # Define the orbital elements
   # a = a*1000 + 6378130  # Semi-major axis in meters
   # e = 0.0  # Eccentricity
   # i = 38.0  # Inclination in degrees
   # rt_asc = 90.0  # Right ascension of the ascending node in degrees
   # arg_p = 60.0  # Argument of periapsis in degrees
   # theta = 0.0  # True anomaly in degrees

   a = 647236 + 6378130  # Semi-major axis in meters
   e = 0.0  # Eccentricity
   i = -33.0  # Inclination in degrees
   rt_asc = 114.0  # Right ascension of the ascending node in degrees
   arg_p = 0.0  # Argument of periapsis in degrees
   theta = 0.0  # True anomaly in degrees

   # Simulate the orbit
   r, v, t = simulateOrbit(a, e, i, rt_asc, arg_p, theta)

   return r, v, t


if __name__ == "__main__":
    mainOrbit()