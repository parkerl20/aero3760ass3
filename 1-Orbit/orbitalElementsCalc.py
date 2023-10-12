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
    print(L_node)