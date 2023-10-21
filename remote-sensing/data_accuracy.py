import numpy as np
import matplotlib.pyplot as plt


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
    H = 1000    # km

    elevation_deg = 10
    elevation_rad = np.deg2rad(elevation_deg)

    delta_phi_deg = 0.06
    delta_phi_rad = np.deg2rad(delta_phi_deg)

    delta_eta_deg = 0.03
    delta_eta_rad = np.deg2rad(delta_eta_deg)

    p_rad = np.arcsin(R_E/(R_E+H))
    p_deg = np.rad2deg(p_rad)

    eta_rad = np.arcsin(np.cos(elevation_rad) * np.sin(p_rad))
    eta_deg = np.rad2deg(eta_rad)

    lam_deg = 90 - elevation_deg - p_deg
    lam_rad = np.deg2rad(lam_deg)

    division = 0.4335692286  # manually calculated because Python sucks.

    D = R_E * division

    azimuth_error = delta_phi_rad * D * np.sin(eta_rad)
    nadir_error = delta_eta_rad * D / np.sin(elevation_rad)

    print("Azimuth error: ", azimuth_error)
    print("Nadir angle error: ", nadir_error)


    return 0


def nadirMappingError():

    # Set the mean and standard deviations
    mean = 0
    std_dev_x = 1.5     # North Error (m)
    std_dev_y = 2.0     # East Error (m)
    num_points = 10000  # Observations taken

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


def swathEdgeMappingError():
    return 0


def singlePointError():

    # Generating data with normal distribution and some noise
    mean = 0
    std_dev = 0.02
    size = 10000
    data = np.random.normal(mean, std_dev, size)
    noise = np.random.normal(0, 1, size)
    data_with_noise = data + noise

    # Plotting the histogram
    plt.figure(figsize=(8, 6))
    plt.hist(data_with_noise, bins=30, density=True, alpha=0.7, color='g')
    plt.title('Histogram of Data with Normal Distribution and Noise', fontsize=16)
    plt.xlabel('Mapping Accuracy Error (m)', fontsize=12)
    plt.ylabel('Frequency', fontsize=12)
    plt.grid(axis='y', alpha=0.75)
    plt.show()



def main():
    
    spatialRes = calculateSpatialRes()
    print("Spatial Resolution: ", spatialRes)

    calculateMappingError()
    # singlePointError()
    # nadirMappingError()
    swathEdgeMappingError()


if __name__ == "__main__":
    main()