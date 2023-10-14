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




def find_nlls_gs_estimate(pseudo_sat, pos_sat, pos_init):
    """
    Finds Non Linear Least Squares estimate of ground station position, and PDOP
    Inputs:
        pseudo_sat: pseudorange data for each satellite
        pos_sat: ECI positions of satellite
        pos_init: initial guess of ground station in ECI frame

    Outputs:
        pos_opt: final position estimate of ground station in LLH
        pdop: PDOP for NLLSE
        pos_store: All the computed position estimates over the iterations
        i_nlls: Number of iterations for convergence
    """
    
    # Initialise arrays
    pos_store = np.array(pos_init)
    fx = np.zeros_like(pseudo_sat)
    dpdx = np.zeros_like(pseudo_sat)
    dpdy = np.zeros_like(pseudo_sat)
    dpdz = np.zeros_like(pseudo_sat)

    # Initial iteration params
    i_nlls = 0
    max_iter = 200
    tol =  1e-8
    dpos = 100
    pos_opt = pos_init

    # until converges
    while (i_nlls < max_iter) and (np.sum(np.abs(dpos) > tol)):
        
        # for each satellite find calculated range fx based on previous estimate
        for i in range(0, pseudo_sat.size):
            fx[i] = np.sqrt((pos_sat[i,0] - pos_store[0,i_nlls])**2 + (pos_sat[i, 1] - pos_store[1,i_nlls])**2\
                            + (pos_sat[i,2] - pos_store[2,i_nlls])**2)
            
        # difference between measured and calculated
        dp = pseudo_sat.T - fx.T

        # Partial differentiation of p - building H matrix (without clock bias)
        for i in range(0, pseudo_sat.size):
            dpdx[i] = -(pos_sat[i,0] - pos_store[0,i_nlls]) / fx[i]
            dpdy[i] = -(pos_sat[i,1] - pos_store[1,i_nlls]) / fx[i]
            dpdz[i] = -(pos_sat[i,2] - pos_store[2,i_nlls]) / fx[i]


        H = np.array([dpdx, dpdy, dpdz]).T
        # Equal weights matrix
        W = np.eye(pseudo_sat.size)
        
        # PDOP        
        pdop = np.sqrt(np.trace( np.linalg.inv( H.T @ H )))

        # NLLS calculation for change in position for next iterations
        dpos = np.linalg.inv(H.T @ W @ H) @ H.T @ W @ dp
        dpos = dpos.reshape((pos_opt.shape))
        

        if np.sum(np.abs(dpos) > tol):
            pos_opt = pos_opt + dpos
            i_nlls +=1
            pos_store = np.append(pos_store, pos_opt, axis = 1)
        else:
            break
    
    print()
    print(f"Number of iterations to converge: {i_nlls}")
    print("PDOP:",pdop)
    return pos_opt, pdop, pos_store, i_nlls
