import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial.transform import Rotation


'''

NLLS for Static Attitude Determination

'''
"""
Planning for function
"""
# Euler angles

# m = lgcv reference vectors


# # reference vectors in LGCV

# # for each vector measurement
# # m is vector measurements
# # number of measurements
# n = 2
# m = np.zeros(n, 3)


def euler2quat(x):
    """
    z,y,x
    yaw, pitch, roll
    psi, theta, phi
    """

    yaw = np.deg2rad(x[0])
    pitch = np.deg2rad(x[1])
    roll = np.deg2rad(x[2])
    qx = np.sin(roll/2) * np.cos(pitch/2) * np.cos(yaw/2) - np.cos(roll/2) * np.sin(pitch/2) * np.sin(yaw/2)
    qy = np.cos(roll/2) * np.sin(pitch/2) * np.cos(yaw/2) + np.sin(roll/2) * np.cos(pitch/2) * np.sin(yaw/2)
    qz = np.cos(roll/2) * np.cos(pitch/2) * np.sin(yaw/2) - np.sin(roll/2) * np.sin(pitch/2) * np.cos(yaw/2)
    qw = np.cos(roll/2) * np.cos(pitch/2) * np.cos(yaw/2) + np.sin(roll/2) * np.sin(pitch/2) * np.sin(yaw/2)
    
    return np.array([ qx, qy, qz,qw])


def quat2euler(q):
    # roll (x-axis rotation)
    #x,y,z,w
    (x, y, z, w) = (q[0], q[1], q[2], q[3])

    t0 = 2 * (w * x + y * z)
    t1 = 1 - 2 * (x * x + y * y)
    roll = np.arctan2(t0, t1)
    t2 = 2 * (w * y - z * x)
    t2 = 1 if t2 > 1 else t2
    t2 = -1 if t2 < -1 else t2
    pitch = np.arcsin(t2)
    t3 = 2 * (w * z + x * y)
    t4 = 1 - 2 * (y * y + z * z)
    yaw = np.arctan2(t3, t4)


    yaw_deg = np.rad2deg(yaw)
    pitch_deg = np.rad2deg(pitch)
    roll_deg = np.rad2deg(roll)
    return np.array([yaw_deg, pitch_deg, roll_deg])


def nlls_quaternion_weights(vector_obs, ref_vectors_lgcv, att_init):
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

    # initial guess in quaternion
    att_opt = euler2quat(att_init)
    print("Initial guess:", att_init)

    # Initial iteration params
    i = 0
    max_iter = 100
    tol =  1e-1
    datt = 1
    datts = []

    while (i < max_iter) and (np.linalg.norm(datt) > tol):
        
        # for each sensor????
        # how does this correspond to satellites from before
        # like how do we get our input data???

        # for each angle find calculated  based on previous estimate
        qx = att_opt[0]
        qy = att_opt[1]
        qz = att_opt[2]
        qw = att_opt[3]

        C_lgcv2body_quat = np.array([[qw**2 + qx**2 - qy**2 - qz**2, 2*(qx*qy - qw*qz), 2*(qx*qz + qw*qy)],
                                [2*(qx*qy + qw*qz), qw**2 - qx**2 + qy**2 - qz**2, 2*(qy*qz - qw*qx)],
                                [2*(qx*qz - qw*qy), 2*(qy*qz + qw*qx), qw**2 - qx**2 - qy**2 + qz**2]]).T

        # print(ref_vectors_lgcv)
        y = C_lgcv2body_quat @ ref_vectors_lgcv.T

        # print("y", y)

        # y is in body frame
        # we make some observation in body frame
        # then we have some reference vectors in lgcv
        # we use the angles we have to find the equivalent body frame measurements

        #y is vectors y[:,i] is a vector 

        # print("y", y.T)

        


        # print(dy)

        
        # # Building H matrix for each vector measurement
        # H is derivative of C_lgcv2body matrix
        # 𝑄(𝑞∘𝑝)⋅𝐼∗ + 𝑄̂ (𝑝∘𝑞−1)
 

        H = []
        # for each vector measurement in body frame
        # or for each vector measurement in lgcv frame
        for m in y.T:
            H11 = 2*(qw * m[0] - qz * m[1] + qy * m[2])
            H12 = 2*(qx * m[0] + qy * m[1] + qz * m[2])
            H13 = 2*(- qy * m[0] + qx * m[1] + qw * m[2])
            H14 = 2*(- qz * m[0] - qw * m[1] + qx * m[2])

            H21 = 2*(qz*m[0] + qw*m[1] - qx*m[2])
            H22 = 2*(qy*m[0] - qx*m[1] - qw*m[2])
            H23 = 2*(qx*m[0] + qy*m[1] + qz*m[2])
            H24 = 2*(qw*m[0] - qz*m[1] + qy*m[2])

            H31 = 2*(-qy*m[0] + qx*m[1] + qw*m[2])
            H32 = 2*(qz*m[0] + qw*m[1] - qx*m[2])
            H33 = 2*(-qw*m[0] + qz*m[1] -qy*m[2])
            H34 = 2*(qx*m[0] + qy*m[1] + qz*m[2])
 
            H.append(np.array([
                [H11, H12, H13, H14],
                [H21, H22, H23, H24],
                [H31, H32, H33, H34]
            ]))

        H = np.array(H)
        # print(H.shape)
        # Reshaping to be 2 dimensional
        H = H.reshape(vector_obs.size, 4)
        # make H for each sensor
        # H = np.array([H,H,H])

        # Equal weights matrix for 
        W = np.eye(vector_obs.size)

        # difference between measured and calculated
        dy = (vector_obs - y)
        # print(y.T)
        # print("vector obs", vector_obs)
        # print("dy",dy)
        # print("dy", dy.shape)
        dy = dy.reshape(dy.size)
         
        pdop = np.sqrt(np.trace( np.linalg.inv( H.T @ H )))
        datt = np.linalg.inv(H.T @ W @ H) @ H.T @ W @ dy

        # this will  output radians
        # print(datt)
        # print(datt_deg)
        # confused about dimension sizing
        # datt = datt.reshape((att_opt.shape))

        print("Euler: ", quat2euler(att_opt))
        print("Shapes:", att_opt.shape, datt.shape)

        att_opt += datt
        i += 1

        datts.append(np.abs(np.linalg.norm(datt)))
        att_store = np.append(att_store, [quat2euler(att_opt)], axis = 0)
    
    print()
    print(f"Number of iterations to converge: {i}")
    print("PDOP: ", pdop)
    return att_opt, pdop, att_store, i, datts


# att_init = (np.array([11,32,-45])[)
att_init = (np.array([0, 0, 0]))
vector_obs = np.array([[-0.0879, 0.5242, -0.6383], [-0.3319, 0.3281, 0.2055], [0.8465, 0.0540, 0.6485]])
ref_vectors_lgcv = np.array([[0.2,0.7,-0.4],[0.1,0.3,0.4],[0.7,-0.8,0.1] ])
# ref_vectors_lgcv = np.array([np.array([0.2, 0.7, -0.4]), np.array([0.1, 0.3, 0.4]), np.array([0.7, -0.8, 0.1])])

# needs to be in radians
att_opt, pdop, att_store, i, datts = nlls_quaternion_weights(vector_obs, ref_vectors_lgcv, att_init)

print("Initial attitude estimation: ", att_init)
print("Final attitude estimation; ", quat2euler(att_opt))

plt.figure()
plt.plot(range(len(att_store)), att_store[:,0], label = "Yaw")
plt.plot(range(len(att_store)), att_store[:,1], label = "Pitch")
plt.plot(range(len(att_store)), att_store[:,2], label = "Roll")
plt.legend()
plt.xlabel("Iterations")
plt.ylabel("Value (degrees)")
plt.show()

plt.figure()
plt.plot(range(i), datts, label = "Change in x")
plt.legend()
plt.xlabel("Iterations")
plt.ylabel("Value")
plt.show()

