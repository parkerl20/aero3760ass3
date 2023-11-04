import numpy as np
import matplotlib.pyplot as plt
from . import mapping_error


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




def nadirMappingError(std_dev_x, std_dev_y, num_points):

    std_dev_x *= 1000  # Convert from km to m
    std_dev_y *= 1000  # Convert from km to m

    # Set the mean and standard deviations
    mean = 0

    # Generate random data points
    east_errors = np.random.normal(mean, std_dev_x, num_points)
    north_errors = np.random.normal(mean, std_dev_y, num_points)

    # Calculate standard deviation
    std_dev = np.sqrt(east_errors ** 2 + north_errors ** 2)
    
    print("Nadir std dev: ", np.mean(std_dev))

    # Create the plot    
    fig = plt.figure(figsize=(8, 6))
    plt.scatter(east_errors, north_errors, c=std_dev, cmap='viridis', marker='o')
    plt.colorbar(label='Standard Deviation')
    plt.axhline(y=0, color='k', linestyle='--', linewidth=0.5)
    plt.axvline(x=0, color='k', linestyle='--', linewidth=0.5)
    plt.title('Nadir Mapping Errors', fontsize=16)
    plt.xlabel('East Error (m)', fontsize=12)
    plt.ylabel('North Error (m)', fontsize=12)
    plt.grid(True, linestyle='--', alpha=0.5)
    fig.savefig('./figures/rs_nadir_error.png')
    
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

    # Calculate standard deviation
    std_dev = np.sqrt(east_errors ** 2 + north_errors ** 2)
    
    print("Swath edge std dev: ", np.mean(std_dev))

    # Create the plot
    plt.figure(figsize=(8, 6))
    plt.scatter(east_errors, north_errors, c=std_dev, cmap='viridis', marker='o')
    plt.colorbar(label='Standard Deviation')
    plt.axhline(y=0, color='k', linestyle='--', linewidth=0.5)
    plt.axvline(x=0, color='k', linestyle='--', linewidth=0.5)
    plt.title('Edge of FOV Mapping Errors', fontsize=16)
    plt.xlabel('East Error (m)', fontsize=12)
    plt.ylabel('North Error (m)', fontsize=12)
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.savefig('./figures/rs_swath_error.png')

    plt.show()



def main():

    observations = 5000
    
    calculateCameraFOV()

    spatialRes = calculateSpatialRes()
    print("Spatial Resolution: ", spatialRes)

    roll_residual_std = 0.01616645
    pitch_residual_std = 0.0134385
    yaw_residual_std = 0.020882695
    cross_track_error = 0.264


    nadirMappingError(pitch_residual_std, yaw_residual_std, observations)
    swathEdgeMappingError(cross_track_error, cross_track_error, observations)

    


if __name__ == "__main__":
    main()