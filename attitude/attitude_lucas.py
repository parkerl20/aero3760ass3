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

    pos_est = np.array([0.0,0.0,0.0])    # Initial guess for the satellite position and clock bias
    
    rel_tol = 1e-10     # Initialising acceptable error
    del_x   = 1         # Setting current residual to be necessarily high

    iter     = 0        # Initialising iteration counter
    max_iter = 100      # Setting max number of iterations

    # While error less than tolerance continue to iterate
    while np.linalg.norm(del_x) > rel_tol:
        # Calculating pseudoranges analytically using the estimated position         
        rho_anal = np.sqrt((satellites[:,0] - pos_est[0])**2 +                   
                           (satellites[:,1] - pos_est[1])**2 +                   
                           (satellites[:,2] - pos_est[2])**2)    
        
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
            return np.array([np.nan, np.nan, np.nan]) , np.nan                      

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

    return pos_est, pdop



def calculate_attitude(noise_bool):
    
    # Obtain Euler angles (an absolute truth)
    # Unit vector given by sun sensor or star tracker
    
    # Add noise based on a standard deviation inaccuracy given by the sensor
    # Add noise based on thermal, space, etc.
    
    # Iterate through NLLS
    
    start_index = 0
    end_index = 1000

    # extract range values to use as psuedo-ranges
    psuedo_range = np.zeros((4,int(end_index-start_index)))
    for i in range(0,4):
        psuedo_range[i,:] = ltpp[i,2,start_index:end_index]


    # extract corresponding ecef satellite positions
    satellite_ecef_pos = np.zeros((4,3,int(end_index-start_index)))
    for i in range(0,4):
        satellite_ecef_pos[i,:,:] = ecef_pos[i,:,start_index:end_index]


    pos_est = np.zeros((3,int(end_index-start_index)))
    pos_acc = np.zeros(int(end_index-start_index))
    pdop = np.zeros(int(end_index-start_index))
    t_over = t_step[0,start_index:end_index]

    # calculate ground station position estimate and gdop
    for i in range(0,len(psuedo_range[0,:])):

        if noise_bool:
            rho_obs = psuedo_range[:,i] + np.array([np.random.normal(0,0.1),np.random.normal(0,0.1),np.random.normal(0,0.1),np.random.normal(0,0.1)])     
        else:
            rho_obs = psuedo_range[:,i]

        satellites = np.array([
            satellite_ecef_pos[0, :, i],
            satellite_ecef_pos[1, :, i],
            satellite_ecef_pos[2, :, i],
            satellite_ecef_pos[3, :, i]
        ])

        pos_est[:,i], pdop[i] = hp.nlls_gps_3(rho_obs,satellites)
        pos_acc[i] = np.linalg.norm(pos_est[:,i] - ground_station_ecef)

    # fig = plt.figure()
    # ax = fig.add_subplot()
    # ax.plot(t_over, pdop,  color = 'r',linewidth=0.5,label = 'gibbs')
    # ax.legend()
    # plt.grid()

    # plt.show()

    return pos_acc, t_over, pdop
