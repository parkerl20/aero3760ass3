import numpy as np
import matplotlib.pyplot as plt
import math


def calculateCameraFOV():

    # Data from the microHSI 425 Sensor

    H_size = 15 * 0.001 # micrometers into mm
    V_size = 15 * 0.001 # micrometers into mm

    focal = 25 # mm

    H_pixels = 640
    V_pixels = 512

    h = H_size * H_pixels
    w = V_size * V_pixels

    FOV_H = 2 * np.arctan(w/(2*focal))
    FOV_V = 2 * np.arctan(h/(2*focal))

    print("Field of View, Vertical: ", np.rad2deg(FOV_V))

    return 0


def calculateSpatialRes():

    n = 4                        # Number of cameras
    R_E = 6371000                # Earth radius
    lambda_max = np.deg2rad(6)   # Swath Width
    beta = np.deg2rad(21.7)      # Camera FOV in degrees
    FOV = beta * n               # Total FOV

    h = R_E * (np.sin(0.5*(lambda_max+FOV)) - np.sin(FOV/2)) / np.sin(FOV/2)

    print("Altitude: ", h/1000, "km")

    swath_width = 40075/360 * np.rad2deg(lambda_max)        # swath width in kilometres

    alpha = swath_width / h

    return alpha


def calculateMappingError():

    R_E = 6371  # km
    H = 500    # km

    elevation_deg = 60
    elevation_rad = np.deg2rad(elevation_deg)

    delta_phi_deg = 0.0005
    delta_phi_rad = np.deg2rad(delta_phi_deg)

    delta_eta_deg = 0.0005
    delta_eta_rad = np.deg2rad(delta_eta_deg)

    p_rad = np.arcsin(R_E/(R_E+H))
    p_deg = np.rad2deg(p_rad)

    eta_rad = np.arcsin(np.cos(elevation_rad) * np.sin(p_rad))
    eta_deg = np.rad2deg(eta_rad)

    lam_deg = 90 - elevation_deg - p_deg
    lam_rad = np.deg2rad(lam_deg)

    print("Lambda:", lam_deg, "eta:", eta_deg)
    print("FOR MANUAL CALCULATION: sin(lambda):", np.sin(lam_rad), "sin(eta):", np.sin(eta_rad))

    # division = 0.4335692286 for elevation 10 deg -- manually calculated because Python sucks.
    division = 0.0451649 # for elevation 60 deg (worst case)

    D = R_E * division

    # Calculate Mapping Error

    print("\n\n Mapping Error \n\n")

    azimuth_error = delta_phi_rad * D * np.sin(eta_rad)
    nadir_error = delta_eta_rad * D / np.sin(elevation_rad)

    print("Azimuth error: ", azimuth_error)
    print("Nadir angle error: ", nadir_error)

    R_T = R_E
    R_S = R_E + H

    delta_I = 0.36 * 0.01 * 0.001 #0.00134
    delta_C = delta_I #0.00134
    delta_R = delta_I #0.00267

    delta_Rt = 0.0001

    phi_deg = 30
    phi_rad = np.deg2rad(phi_deg)

    H_mapping = np.arcsin(np.sin(lam_rad)*np.sin(phi_rad))
    G_mapping = np.arcsin(np.sin(lam_rad)*np.cos(phi_rad))
 
    intrack_error = delta_I * R_T / R_S * np.cos(H_mapping)
    crosstrack_error = delta_C * R_T / R_S * np.cos(G_mapping)
    radial_error = delta_R * np.sin(eta_rad) / np.sin(elevation_rad)

    target_altitude = delta_Rt / np.tan(elevation_rad)

    print("In-track error: ", intrack_error)
    print("Cross-track error: ", crosstrack_error)
    print("Radial error: ", radial_error)
    print("Target altitude: ", target_altitude)

    # Sample data
    data = [azimuth_error, nadir_error, intrack_error, crosstrack_error, radial_error, target_altitude]

    # Calculate RMS
    rms_result = calculate_rms(data)
    print(f"The root mean square (FINAL MAPPING ACCURACY) is {rms_result:.6f}.")

    # Calculate Pointing Error

    '''
        SMAD states that,
            say our esnsor has a 1 deg FOV, if our pointing error is
            at least four times less than that (i.e. 0.25 deg), we 
            can be assured the target is within the FOV on essentially 
            every observation (within 6SD, which is virtually always).
        Page 212 of PDF, 195 of SMAD
    '''

    print("\n\n Pointing Error \n\n")

    azimuth_pointing_error = delta_phi_deg * np.sin(eta_rad)
    nadir_pointing_error = delta_eta_deg

    print("Azimuth pointing error: ", azimuth_pointing_error)
    print("Nadir angle pointing error: ", nadir_pointing_error)

    Y_I = np.arccos(np.cos(phi_rad)*np.sin(eta_rad))
    Y_C = np.arccos(np.sin(phi_rad)*np.sin(eta_rad))
    
    intrack_pointing = delta_I / D * np.sin(Y_I)
    crosstrack_pointing = delta_C / D * np.sin(Y_C)
    radiai_pointing = delta_R / D * np.sin(eta_rad)

    print("In-track error: ", intrack_pointing)
    print("Cross-track error: ", crosstrack_pointing)
    print("Radial error: ", radiai_pointing)
    print("Target altitude: N/A ")

    # Sample data
    data_pointing = [azimuth_pointing_error, nadir_pointing_error, 
                     intrack_pointing, crosstrack_pointing, radiai_pointing]

    # Calculate RMS
    rms_result_pointing = calculate_rms(data_pointing)
    print(f"The root mean square (FINAL POINTING ACCURACY) is {rms_result_pointing:.6f}.")

    return delta_I, delta_C, delta_R



def nadirMappingError(std_dev_x, std_dev_y, num_points):

    std_dev_x *= 1000 # Convert from km to m
    std_dev_y *= 1000 # Convert from km to m

    # Set the mean and standard deviations
    mean = 0

    # Generate random data points
    east_errors = np.random.normal(mean, std_dev_x, num_points)
    north_errors = np.random.normal(mean, std_dev_y, num_points)

    # Create the plot
    plt.figure(figsize=(8, 6))
    plt.scatter(east_errors, north_errors, facecolors='none', edgecolors='b')
    plt.axhline(y=0, color='k', linestyle='--', linewidth=0.5)
    plt.axvline(x=0, color='k', linestyle='--', linewidth=0.5)
    plt.title('Nadir Mapping Errors', fontsize=16)
    plt.xlabel('East Error (m)', fontsize=12)
    plt.ylabel('North Error (m)', fontsize=12)
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.show()


def swathEdgeMappingError(std_dev_x, std_dev_y, num_points):

    np.random.seed(0)  # Setting a seed for reproducibility

    # Set the mean and standard deviations
    mean = 0
    std_dev_x *= 2500 # Convert from km to m
    std_dev_y *= 2500 # Convert from km to m

    # Generate random data points
    east_errors = np.random.normal(mean, std_dev_x, num_points)
    north_errors = np.random.normal(mean, std_dev_y, num_points)

    # Create the plot
    plt.figure(figsize=(8, 6))
    plt.scatter(east_errors, north_errors, facecolors='none', edgecolors='b')
    plt.axhline(y=0, color='k', linestyle='--', linewidth=0.5)
    plt.axvline(x=0, color='k', linestyle='--', linewidth=0.5)
    plt.title('Edge of FOV Mapping Errors', fontsize=16)
    plt.xlabel('East Error (m)', fontsize=12)
    plt.ylabel('North Error (m)', fontsize=12)
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.show()


def singlePointError(std_dev, size):

    # Generating data with normal distribution and some noise
    mean = 0
    data = np.random.normal(mean, std_dev, size)
    noise = np.random.normal(0, 1, size)
    data_with_noise = data + noise

    # Plotting the histogram
    plt.figure(figsize=(8, 6))
    plt.hist(data_with_noise, bins=30, density=True, alpha=0.7, color='g')
    plt.title('Monte Carlo on Mapping Accuracy Confidence on One Point', fontsize=16)
    plt.xlabel('Mapping Accuracy Error (m)', fontsize=12)
    plt.ylabel('Frequency', fontsize=12)
    plt.grid(axis='y', alpha=0.75)
    plt.show()


def calculate_rms(data):
    """
    Calculate the root mean square of a list of values.
    """
    n = len(data)
    if n < 1:
        raise ValueError("Input list must contain at least one element.")

    sum_of_squares = sum(x ** 2 for x in data)
    rms = math.sqrt(sum_of_squares / n)
    return rms



def main():

    observations = 5000
    
    calculateCameraFOV()

    spatialRes = calculateSpatialRes()
    print("Spatial Resolution: ", spatialRes)

    std_dev_x, std_dev_y, std_dev_radial = calculateMappingError()
    # singlePointError(std_dev_radial, observations)
    # nadirMappingError(std_dev_x, std_dev_y, observations)
    # swathEdgeMappingError(std_dev_x, std_dev_y, observations)

    


if __name__ == "__main__":
    main()