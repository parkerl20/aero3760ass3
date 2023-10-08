import position_truth as position_truth

def main():
    
    print("Starting Position determination...")
    
    position_truth.generate_truth_orbit()
    
    '''
    
    1. Find ECEF State vector from GNSS for the postion measurement
    2. Find ECEF State vector from SLR for the noise model
    3. Using Extended Kalman filter class choose process, state and observation covariance
    4. Pass each respective state vector through the EKF
    5. Reach State Vector Estimate in ECEF frame 
    
    return estimatedPosition
    
    '''
    
    return 0


if __name__ == "__main__":
    main()
    