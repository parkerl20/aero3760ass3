import numpy as np

'''

NLLS for Static Attitude Determination

'''

def nlls_gps(rho_obs,satellites):
    '''
    NLLS_GPS - calculates the position of a ground station given pseudorange measurements
    and satellite positions using Non-Linear Least Sqaures and the Pseudorange Equation.

    Inputs: rho_obs (numpy.ndarray) = contains pseudorange measurements from different satellites (n x 1 numpy array)
            satellites (numpy.ndarray) = contains the satellites that took the respectiv pseudorange measurements (n x 3 numpy array)

    Outputs: pos_est (numpy.ndarray) = contains the final estimation of the ground station (1x3 numpy array)

    '''

    pos_est = np.array([0.0,0.0,0.0,0.0])    # Initial guess for the satellite position and clock bias
    
    rel_tol = 1e-10     # Initialising acceptable error
    del_x   = 1         # Setting current residual to be necessarily high

    iter     = 0        # Initialising iteration counter
    max_iter = 100      # Setting max number of iterations

    # While error less than tolerance continue to iterate
    while np.linalg.norm(del_x) > rel_tol:
        # Calculating pseudoranges analytically using the estimated position         
        rho_anal = np.sqrt((satellites[:,0] - pos_est[0])**2 +                   
                           (satellites[:,1] - pos_est[1])**2 +                   
                           (satellites[:,2] - pos_est[2])**2) + pos_est[3]    
        
        # Initialising Jacobian of NLLS
        H = np.ones((len(rho_obs),4))                                               

        for i in range(0,len(rho_obs)):

            # Calculating Jacobian and filling H
            H[i,0:3] = np.array([                                                   
                -(satellites[i,0] - pos_est[0]),                                                         
                -(satellites[i,1] - pos_est[1]),                                                         
                -(satellites[i,2] - pos_est[2])                                                          
            ]) / np.sqrt(                                                           
                    (satellites[i,0] - pos_est[0])**2 +                                                  
                    (satellites[i,1] - pos_est[1])**2 +                          
                    (satellites[i,2] - pos_est[2])**2                            
                )                                                                   

        # Calculating difference between observed pseudoranges and calculated pseudoranges
        del_rho = rho_obs - rho_anal    

        # Calculating residual between estimated position and true position using NLLS       
        try:                                    
            del_x = np.linalg.inv(H.T @ H) @ H.T @ del_rho        
        except:
            return np.array([np.nan, np.nan, np.nan, np.nan]) , np.nan                      

        # Adjust initial guess
        pos_est += del_x                                                         
        iter += 1       


        # If iteration count exceed max iterations break loop
        if iter > max_iter:
            print("iteration could not meet tolerance")
            break

        # Geometric Dilution of Precision
        V = np.linalg.inv(H.T @ H)

        pdop = np.sqrt(V[0,0] + V[1,1] + V[2,2])
        gdop = np.sqrt(V[0,0] + V[1,1] + V[2,2] + V[3,3])

    return pos_est, gdop
