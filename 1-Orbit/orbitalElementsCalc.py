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


def minimumAltitude(swathe_width, alpha):
    '''
    minimumAltitude - Calculates the minimium altitude a satellite can be positioned at above ground while still 
                      covering the required swathe_width, based on the spatial resolution of the sensor
    Inputs: swathe_width (float) = Degrees of latitude that the satellite must be able to observe
            alpha (float) = Spatial resoltuion of camera

    Outputs: alt (float) = Minimum altitude the satellite can be placed at
    '''
    altitude = swathe_width/alpha
    return altitude