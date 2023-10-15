import matplotlib.pyplot as plt
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

# # for each vector measurement
# # m is vector measurements
# # number of measurements
# n = 2
# m = np.zeros(n, 3)





def nlls_gps_weights(vector_obs, ref_vectors_lgcv, att_init):
    """
    known vectors ref_vectors_lgcv m
    y measured vectors in body frame vector_obs
    should be performing nlls on y which is in the body frame vectors

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
    att_store = np.array([att_init])


    # fx = np.zeros_like(pseudo_sat)
    # dpdx = np.zeros_like(pseudo_sat)
    # dpdy = np.zeros_like(pseudo_sat)
    # dpdz = np.zeros_like(pseudo_sat)

    # Initial iteration params
    i = 0
    max_iter = 100
    tol =  1e-8
    datt = 100
    att_opt = att_init

    # until converges
    while (i < max_iter) and (np.sum(np.abs(datt) > tol)):
        
        # for each sensor????
        # how does this correspond to satellites from before
        # like how do we get our input data???

        # for each angle find calculated  based on previous estimate
        psi = np.deg2rad(att_opt[0])
        theta = np.deg2rad(att_opt[1])
        phi =  np.deg2rad(att_opt[2])


        C_lgcv = np.array([[np.cos(psi)*np.cos(theta), np.sin(psi)*np.cos(theta), -np.sin(theta)],
                      [np.cos(psi)*np.sin(theta)*np.sin(phi) - np.sin(psi)*np.cos(phi), np.sin(psi)*np.sin(theta)*np.sin(phi) + np.cos(psi)*np.cos(phi), np.cos(theta)*np.sin(phi)],
                    [np.cos(psi)*np.sin(theta)*np.cos(phi) + np.sin(psi)*np.sin(phi), np.sin(psi)*np.sin(theta)*np.cos(phi) - np.cos(psi)*np.sin(phi), np.cos(theta)*np.cos(phi)]])
        
        y = C_lgcv @ ref_vectors_lgcv.T

        # difference between measured and calculated
        # this is wrong - basically i don't understand how do more iterations than i have data??
        # look over other nlls again 
        # should take all att data - 
        dy = (vector_obs.T - y)
        # i think this is wrong dimension

        print(dy)
        print("dy", dy.shape)
        dy = dy.reshape(dy.size)

        # Building H matrix for each vector measurement
        
        H = []
        for m in vector_obs:
            H11 = 0
            H12 = -np.cos(psi)*np.sin(theta)*m[0] - np.sin(psi)*np.sin(theta)*m[1] - np.cos(theta)*m[2]
            H13 = -np.sin(psi)*np.cos(phi)*m[0] + np.cos(psi)*np.cos(phi)*m[1]

            H21 = (np.cos(psi)*np.sin(theta)*np.cos(phi) + np.sin(psi)*np.sin(phi))*m[0] + (np.sin(psi)*np.sin(theta)*np.cos(phi) - np.cos(psi)*np.sin(phi))*m[1] + np.cos(theta)*np.cos(phi)*m[2]
            H22 = np.cos(psi)*np.cos(theta)*np.sin(phi)*m[0] + np.sin(psi)*np.cos(theta)*np.sin(phi)*m[1] - np.sin(theta)*np.sin(phi)*m[2]
            H23 = (-np.sin(psi)*np.sin(theta)*np.sin(phi) - np.cos(psi)*np.cos(phi))*m[0] + (np.cos(psi)*np.sin(theta)*np.sin(phi) - np.sin(psi)*np.cos(phi))*m[1]

            H31 = (-np.cos(psi)*np.sin(theta)*np.sin(phi) + np.sin(psi)*np.cos(phi))*m[0] + (-np.sin(psi)*np.sin(theta)*np.sin(phi) - np.cos(psi)*np.cos(phi))*m[1] - np.cos(theta)*np.sin(phi)*m[2]
            H32 = np.cos(psi)*np.cos(theta)*np.cos(phi)*m[0] + np.sin(psi)*np.cos(theta)*np.cos(phi)*m[1] - np.sin(theta)*np.cos(phi)*m[2]
            H33 = (-np.sin(psi)*np.sin(theta)*np.cos(phi) + np.cos(psi)*np.sin(phi))*m[0] + (np.cos(psi)*np.sin(theta)*np.cos(phi) + np.sin(psi)*np.sin(phi))*m[1]

            # H = np.append(H, [
            # H = np.array([
            H.append(np.array([
                [H11, H12, H13 ],
                [H21, H22, H23],
                [H31, H32, H33]]))
        
        H = np.array(H)
        print(H.shape)
        H = H.reshape(vector_obs.size,3)
        # make H for each sensor
        # H doesn't work because needs to be square
        # H = np.array([H,H,H])



        # Equal weights matrix for 
        W = np.eye(vector_obs.size)
        print(vector_obs.shape)

        # W = np.eye(3)
        
        print(H.shape)
        # PDOP        
        pdop = np.sqrt(np.trace( np.linalg.inv( H.T @ H )))

        # NLLS calculation for change in attitude for next iterations


        # H = np.array([H,H,H])

        print((np.linalg.inv(H.T @ W @ H) @ H.T @ W).shape)
        print(dy)
        datt = np.linalg.inv(H.T @ W @ H) @ H.T @ W @ dy

        # this will probs output radians
        print(datt)

        datt_deg = np.rad2deg(datt)
        print(datt_deg)

        


        # confused about dimension sizing
        # datt = datt.reshape((att_opt.shape))
        

        if np.sum(np.abs(datt) > tol):
            print(att_opt)
            att_opt = att_opt + datt.T
            print(att_opt)
            
            i +=1
            att_store = np.append(att_store,[att_opt], axis = 0)
        else:
            break
    
    print()
    print(f"Number of iterations to converge: {i}")
    print("PDOP:",pdop)
    return att_opt, pdop, att_store, i



att_init = (np.array([11,32,-45]))

# say we have 4 sensors
att_data = np.array([[10,31,-41],[9,30,-47],[14,32,-44]])

vector_obs = np.array([[-0.0879, 0.5242, -0.6383],[-0.3319, 0.3281, 0.2055], [0.8465, 0.0540, 0.6485]])
ref_vectors_lgcv = np.array([[0.2,0.7,-0.4],[0.1,0.3,0.4],[0.7,-0.8,0.1] ])
# needs to be in radians
att_opt, pdop, att_store, i = nlls_gps_weights(vector_obs, ref_vectors_lgcv, att_init)

print("Final attitude estimation", att_opt)
print(att_store.shape)
plt.figure()
plt.plot(range(len(att_store)), att_store[:,0])
plt.plot(range(len(att_store)), att_store[:,1])
plt.plot(range(len(att_store)), att_store[:,2])
plt.show()
# justifying nlls
