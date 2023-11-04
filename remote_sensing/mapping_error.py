import math
import numpy as np


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
            say our sensor has a 1 deg FOV, if our pointing error is
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

