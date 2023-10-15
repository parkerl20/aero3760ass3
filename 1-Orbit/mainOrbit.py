import numpy as np
from orbitalElementsCalc import flyOverRightAscension
from orbitalElementsCalc import minimumAltitude
from simulateOrbit import simulateOrbit
from orbitalElementsCalc import repeatingGroundTrackAltitude
# from spacesim import constants as =

def mainOrbit():
   # Coordinates of the centre of NSW
   lon = 148.5
   lat = -32.0

   # Other paramaters
   swathe_width = 10.0  # Swathe width for one satellite to cover NSW
   alpha = 2.4        # Spacial resolution of the chosen camera

   # Minimum altitude of satelltie to meet swathe_width requirement
   a_min = minimumAltitude(swathe_width, alpha)

   # Altitude for a repeating ground track orbit based on j: orbital periods and k:number of days
   a = repeatingGroundTrackAltitude(j=15, k=1, i=-32.0)

   # Right ascension for a satellite with a required fly over point
   rt_asc = flyOverRightAscension(lat, lon, i=lat)
   # rt_asc = 0

   i = lat  # Inclination in degrees
   e, arg_p, theta = 0.0, 0.0, 0.0  # Other orbital paramaters
   
   # Simulate the orbit
   simulateOrbit(a, e, i, rt_asc, arg_p, theta)

   return 0


if __name__ == "__main__":
   mainOrbit()