import numpy as np
from .orbitalElementsCalc import flyOverRightAscension
from .orbitalElementsCalc import minimumAltitude
from .orbitalElementsCalc import choosingOrbitsPerDay
from .simulateOrbit import simulateOrbit
from .orbitalElementsCalc import repeatingGroundTrackAltitude
from .orbitalElementsCalc import choosingRepeatEccentricity

def mainOrbit():
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
   swathe_width = abs(nsw_bounds["lat_min"] - nsw_bounds["lat_max"])  # Swathe width for one satellite to cover NSW left to right
   i = fly_over["lat"]       # Inclination of satellite in degrees
   alpha = 2.4               # Spacial resolution of the chosen camera
   SINGLE_DAY = 1            # Orbit repeats every 1 day

   # Minimum altitude of satelltie to meet swathe_width requirement
   a_min = minimumAltitude(swathe_width, alpha)

   # Altitude for a repeating ground track orbit based on j: orbital periods and k:number of days
   orbits_per_day = [13, 14, 15, 16]
   j = choosingOrbitsPerDay(a_min, i, orbits_per_day, k=SINGLE_DAY)
   
   # Choose ideal eccentricity for minimum altitude that still meets coverage requirements
   e = choosingRepeatEccentricity(a_min, i, j, SINGLE_DAY)   

   # Chooses 15 orbits a day as the ideal number
   a = repeatingGroundTrackAltitude(j, SINGLE_DAY, i, e)

   # Right ascension for a satellite with a required fly over point
   rt_asc = flyOverRightAscension(fly_over["lat"], fly_over["lon"], i)

   # Other orbital paramaters
   arg_p, theta = 0.0, 0.0  
   
   # Simulate the orbit
   results = simulateOrbit(a, e, i, rt_asc, arg_p, theta)

   return results


if __name__ == "__main__":
   mainOrbit()