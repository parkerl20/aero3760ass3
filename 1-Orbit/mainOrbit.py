import numpy as np
from orbitalElementsCalc import orbitalParamaters

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
    lon = 146.9211
    lat = -31.2535
    swathe_width = 4.609

    # Calculate orbital paramaters based on a latitude, longitude and swathe_width
    orbitalParamaters(lat, lon, swathe_width)

    # Final return skeleton
    eci = np.array([0, 0, 0])
    t = np.linspace(0, 10, 1)

    return eci, t


if __name__ == "__main__":
    mainOrbit()