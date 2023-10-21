import numpy as np
from .orbitalElementsCalc import flyOverRightAscension
from .orbitalElementsCalc import minimumAltitude
from .orbitalElementsCalc import choosingOrbitsPerDay
from .simulateOrbit import simulateOrbit
from .orbitalElementsCalc import repeatingGroundTrackAltitude
from .orbitalElementsCalc import choosingRepeatEccentricity
from .orbitalElementsCalc import minimumAltitudeTest

def mainOrbit(show_results):
   # Boundarys of NSW 
   nsw_bounds = {
    "lat_min": -28.6,
    "lat_max": -37.4,
    "lon_min": 141,
    "lon_max": 156.64
   }

   # Coordinates of middle point to fly over
   fly_over = {
    "lat": (nsw_bounds["lat_min"] + nsw_bounds["lat_max"])/2,
    "lon": (nsw_bounds["lon_min"] + nsw_bounds["lon_max"])/2
   }

   # Other paramaters
   SWATHE_WIDTH = 6.2  # Swathe width for one satellite to cover NSW right to left (calculated through GEE)
   i = fly_over["lat"] - 1.5      # Inclination of satellite in degrees
   alpha = 2.4               # Spacial resolution of the chosen camera
   SINGLE_DAY = 1            # Orbit repeats every 1 day

   # Minimum altitude of satelltie to meet swathe_width requirement
   # a_min = minimumAltitude(SWATHE_WIDTH)
   # print(a_min)
   a_min = minimumAltitudeTest(SWATHE_WIDTH, alpha=2.4)
   # print(a_min)

   # Altitude for a repeating ground track orbit based on j: orbital periods and k:number of days
   orbits_per_day = [13, 14, 15, 16]
   j = choosingOrbitsPerDay(a_min, i, orbits_per_day, SINGLE_DAY, show_results)

   # Choose ideal eccentricity for minimum altitude that still meets coverage requirements
   e = choosingRepeatEccentricity(a_min, i, j, SINGLE_DAY, show_results)   

   # Chooses 15 orbits a day as the ideal number
   a = repeatingGroundTrackAltitude(j, SINGLE_DAY, i, e)

   # Right ascension for a satellite with a required fly over point
   rt_asc = flyOverRightAscension(fly_over["lat"], fly_over["lon"], i)

   # Other orbital paramaters
   arg_p, theta = 0.0, 0.0  
   
   # Simulate the orbit
   results = simulateOrbit(a, e, i, rt_asc, arg_p, theta, show_results)

   return results


if __name__ == "__main__":
   mainOrbit(1)