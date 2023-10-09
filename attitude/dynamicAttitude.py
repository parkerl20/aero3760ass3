import numpy as np


"""
dynamicAttitude.py

Computes dynamic attitude determination using Extended Kalman Filter
To be integrated with static attitude determination using sensor fusion to account for acceration drift
"""


# Get initial guess from static attitude determination algorithm



# -------------- Initialise vectors -----------------------
# 3 axis stabilised system
# state vector: roll pitch yaw + rates from gyro
# m dimensional vector which includes all variables necessary for accuracte attitude determination 
m = 6
X = np.arange(m)

# observation vector from sensors
# generate measurements from observations (depending on sensor, could be multiple obs for one measurement)
y = None # length n
n = 10
# observational model vector - predicted values of obs based on estimated values of state vector elements 
z = np.zeros(n) 
 
# convert to state to quaternions - which are time dependent
# state vector x also has errors: may be constant (eg biases) or time varying (quaternions) 


## --------------  EKF algorithm - iterate until converge -----------------------

# while not converged (check tolerance)
    # get measurement
    # Predict new state
    # jacobians
    # predict new covariance
    # simulate measurement
    # Update observation model vector z
    # Update residual (difference between truth and estimate)
    # Find Kalman gain
    # Update State
    # Update covariance



def main():
    
   
    return 0


if __name__ == "__main__":
    main()