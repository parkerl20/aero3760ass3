
import numpy as np

'''

NLLS for Static Attitude Determination

'''
"""
Planning for function
"""
# Euler angles
# y = np.array([phi, theta, psi])


# m = lgcv reference vectors


# # reference vectors in LGCV
# y = C_lgcv @ m

# C_lgcv = np.array([[np.cos(psi)*np.cos(theta), np.sin(psi)*np.cos(theta), -np.sin(theta)],
#               [np.cos(psi)*np.sin(theta)*np.sin(phi) - np.sin(psi)*np.cos(phi), np.sin(psi)*np.sin(theta)*np.sin(phi) + np.cos(psi)*np.cos(phi), np.cos(theta)*np.sin(phi)],
#             [np.cos(psi)*np.sin(theta)*np.cos(phi) + np.sin(psi)*np.sin(phi), np.sin(psi)*np.sin(theta)*np.cos(phi) - np.cos(psi)*np.sin(phi), np.cos(theta)*np.cos(phi)]])

# # for each vector measurement
# # m is vector measurements
# # number of measurements
# n = 2
# m = np.zeros(n, 3)





def nlls_gps_weights(att_data, att_init):
    """
    Finds Non Linear Least Squares estimate of attitude in Euluer angles, and PDOP
    Inputs:
        
        att_data: attitude data for each satellite in Euluer (we can change this to be lgcv or something using matrix above)
        
        att_init: initial guess of attitude of sat in Euler angles

    Outputs:
        pos_opt: final position estimate of ground station in LLH
        pdop: PDOP for NLLSE
        pos_store: All the computed position estimates over the iterations
        i_nlls: Number of iterations for convergence
    """
    
    # Initialise arrays
    att_store = np.array(att_init)


    # fx = np.zeros_like(pseudo_sat)
    # dpdx = np.zeros_like(pseudo_sat)
    # dpdy = np.zeros_like(pseudo_sat)
    # dpdz = np.zeros_like(pseudo_sat)

    # Initial iteration params
    i = 0
    max_iter = 200
    tol =  1e-8
    datt = 100
    att_opt = np.array([att_init])

    # until converges
    while (i < max_iter) and (np.sum(np.abs(datt) > tol)):
        
        # for each angle find calculated  based on previous estimate
        psi = np.deg2rad(att_opt[-1,0])
        theta = np.deg2rad(att_opt[-1,1])
        phi =  np.deg2rad(att_opt[-1,2])

        # difference between measured and calculated
        dy = att_data - att_opt

        print(dy)

        # Building H matrix

        H11 = 0
        H12 = -np.cos(psi)*np.sin(theta)*att_data[i,0] - np.sin(psi)*np.sin(theta)*att_data[i,1] - np.cos(theta)*att_data[i,2]
        H13 = -np.sin(psi)*np.cos(phi)*att_data[i,0] + np.cos(psi)*np.cos(phi)*att_data[i,1]

        H21 = (np.cos(psi)*np.sin(theta)*np.cos(phi) + np.sin(psi)*np.sin(phi))*att_data[i,0] + (np.sin(psi)*np.sin(theta)*np.cos(phi) - np.cos(psi)*np.sin(phi))*att_data[i,1] + np.cos(theta)*np.cos(phi)*att_data[i,2]
        H22 = np.cos(psi)*np.cos(theta)*np.sin(phi)*att_data[i,0] + np.sin(psi)*np.cos(theta)*np.sin(phi)*att_data[i,1] - np.sin(theta)*np.sin(phi)*att_data[i,2]
        H23 = (-np.sin(psi)*np.sin(theta)*np.sin(phi) - np.cos(psi)*np.cos(phi))*att_data[i,0] + (np.cos(psi)*np.sin(theta)*np.sin(phi) - np.sin(psi)*np.cos(phi))*att_data[i,1]

        H31 = (-np.cos(psi)*np.sin(theta)*np.sin(phi) + np.sin(psi)*np.cos(phi))*att_data[i,0] + (-np.sin(psi)*np.sin(theta)*np.sin(phi) - np.cos(psi)*np.cos(phi))*att_data[i,1] - np.cos(theta)*np.sin(phi)*att_data[i,2]
        H32 = np.cos(psi)*np.cos(theta)*np.cos(phi)*att_data[i,0] + np.sin(psi)*np.cos(theta)*np.cos(phi)*att_data[i,1] - np.sin(theta)*np.cos(phi)*att_data[i,2]
        H33 = (-np.sin(psi)*np.sin(theta)*np.cos(phi) + np.cos(psi)*np.sin(phi))*att_data[i,0] + (np.cos(psi)*np.sin(theta)*np.cos(phi) + np.sin(psi)*np.sin(phi))*att_data[i,1]

        H = np.array([
            [H11, H12, H13 ],
            [H21, H22, H23],
            [H31, H32, H33]])



        # Equal weights matrix
        W = np.eye(att_data.size)
        W = np.eye((3))
        
        print(H.shape)
        # PDOP        
        pdop = np.sqrt(np.trace( np.linalg.inv( H.T @ H )))

        # NLLS calculation for change in attitude for next iterations
        datt = np.linalg.inv(H.T @ W @ H) @ H.T @ W @ dy

        # confused about dimension sizing
        datt = datt.reshape((att_opt.shape))
        

        if np.sum(np.abs(datt) > tol):
            att_opt = att_opt + datt
            i +=1
            att_store = np.append(att_store, att_opt, axis = 1)
        else:
            break
    
    print()
    print(f"Number of iterations to converge: {i}")
    print("PDOP:",pdop)
    return att_opt, pdop, att_store, i



att_init = (np.array([11,32,-45]))
att_data = np.array([[10,31,-41],[9,30,-47],[14,32,-44]])
# needs to be in radians
nlls_gps_weights(att_data, att_init)
