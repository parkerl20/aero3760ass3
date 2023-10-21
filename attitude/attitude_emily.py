import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial.transform import Rotation

'''
NLLS for Attitude Determination
'''

# Euler angles

# m = lgcv reference vectors


# # reference vectors in LGCV

# # for each vector measurement
# # m is vector measurements
# # number of measurements
# n = 2
# m = np.zeros(n, 3)

  


def nlls_euler_weights(vector_obs, ref_vectors_lgcv, att_init):
    """
    known vectors ref_vectors_lgcv m
    y measured vectors in body frame vector_obs
    should be performing nlls on y which is in the body frame vectors

    Finds Non Linear Least Squares estimate of attitude in Euler angles, and PDOP
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
    max_iter =200
    tol =  1e-5
    datt = 100
    att_opt = att_init
    datts = []

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

        #y is vectors y[:,i] is a vector 

        print(y.T)

        # difference between measured and calculated
        dy = (vector_obs - y.T)

        print(y.T)
        print("vector obs", vector_obs)
        print("dy",dy)



        print("dy", dy.shape)
        dy = dy.reshape(dy.size)


        print(dy)

        
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
        # print(H.shape)
        # Reshaping to be 2 dimensional
        H = H.reshape(vector_obs.size,3)
        # make H for each sensor
        # H = np.array([H,H,H])

       # Equal weights matrix for 
        W = np.eye(vector_obs.size)
     
        # W = np.array([
        #     [0.2,0,0,0,0,0,0,0,0,0,0,0],
        #     [0,0.2,0,0,0,0,0,0,0,0,0,0],
        #     [0,0,0.2,0,0,0,0,0,0,0,0,0],
        #     [0,0,0,1,0,0,0,0,0,0,0,0],
        #     [0,0,0,0,1,0,0,0,0,0,0,0],
        #     [0,0,0,0,0,1,0,0,0,0,0,0],
        #     [0,0,0,0,0,0,1,0,0,0,0,0],
        #     [0,0,0,0,0,0,0,1,0,0,0,0],
        #     [0,0,0,0,0,0,0,0,1,0,0,0],
        #     [0,0,0,0,0,0,0,0,0,1,0,0],
        #     [0,0,0,0,0,0,0,0,0,0,1,0],
        #     [0,0,0,0,0,0,0,0,0,0,0,1]

        # ])

        W[0,0] = 0.2
        W[1,1] = 0.2
        W[2,2] = 0.2
        


        # print(vector_obs.shape)

        # W = np.eye(3)
        
        # print(H.shape)
        # PDOP        
        pdop = np.sqrt(np.trace( np.linalg.inv( H.T @ H )))

        # NLLS calculation for change in attitude for next iterations
        # H = np.array([H,H,H])
        # print(dy)
        datt = np.linalg.inv(H.T @ W @ H) @ H.T @ W @ dy

        # this will  output radians
        # print(datt)

        datt_deg = np.rad2deg(datt)
        # print(datt_deg)

        # confused about dimension sizing
        # datt = datt.reshape((att_opt.shape))
        

        if np.sum(np.abs(datt) > tol):
            att_opt = att_opt + datt.T
            datts.append(np.linalg.norm(datt))
            
            i +=1
            att_store = np.append(att_store,[att_opt], axis = 0)
        else:
            break
    
    print()
    print(f"Number of iterations to converge: {i}")
    print("PDOP:",pdop)
    return att_opt, pdop, att_store, i, datts

def euler2quat(x):
    """
    Converts Euler angles to quaternion
    Inputs:
        x: Euler angles in degrees yaw, pitch roll
            in order (z,y,x: yaw, pitch, roll: psi, theta, phi)
    Outputs:
        q: quaternion in form w x y z
    """

    yaw = np.deg2rad(x[0])
    pitch = np.deg2rad(x[1])
    roll = np.deg2rad(x[2])
    qx = np.sin(roll/2) * np.cos(pitch/2) * np.cos(yaw/2) - np.cos(roll/2) * np.sin(pitch/2) * np.sin(yaw/2)
    qy = np.cos(roll/2) * np.sin(pitch/2) * np.cos(yaw/2) + np.sin(roll/2) * np.cos(pitch/2) * np.sin(yaw/2)
    qz = np.cos(roll/2) * np.cos(pitch/2) * np.sin(yaw/2) - np.sin(roll/2) * np.sin(pitch/2) * np.cos(yaw/2)
    qw = np.cos(roll/2) * np.cos(pitch/2) * np.cos(yaw/2) + np.sin(roll/2) * np.sin(pitch/2) * np.sin(yaw/2)
    
    return np.array([ qw,qx, qy, qz])


def quat2euler(q):
    """
    Converts quaternions to Euler angles
    Inputs:
        q: quaternion in form w x y z
        
    Outputs:
        x: Euler angles in degrees yaw, pitch roll
            in order (z,y,x: yaw, pitch, roll: psi, theta, phi)
        
    """
    # roll (x-axis rotation)
    #x,y,z,w
    (w, x, y, z) = (q[0], q[1], q[2], q[3])

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

    # Initial iteration params
    i = 0
    max_iter = 100
    tol =  1e-8
    datt = 100

    # initial guess in quaternion
    att_opt = euler2quat(att_init)
    # att_opt = att_init

    datts = []

    # until converges
    while (i < max_iter) and (np.sum(np.abs(datt) > tol)):
        
        # for each sensor????
        # how does this correspond to satellites from before
        # like how do we get our input data???

        # for each angle find calculated  based on previous estimate
        qw = att_opt[0]
        qx = att_opt[1]
        qy = att_opt[2]
        qz = att_opt[3]

        x_eul = quat2euler(att_opt)
        psi = np.deg2rad(x_eul[0])
        theta = np.deg2rad(x_eul[1])
        phi = np.deg2rad(x_eul[2])

        # To convert LGCV measurements to body frame
        # Take the transpose as this formula is based on Euler lgcv transpose
        # therefore my derivatives are wrong
        # fix quaternions first
        C_lgcv2body_quat = np.array([[qw**2 + qx**2 - qy**2 - qz**2, 2*(qx*qy - qw*qz), 2*(qx*qz + qw*qy)],
                                [2*(qx*qy + qw*qz), qw**2 - qx**2 + qy**2 - qz**2, 2*(qy*qz - qw*qx)],
                                [2*(qx*qz - qw*qy), 2*(qy*qz + qw*qx), qw**2 - qx**2 - qy**2 + qz**2]]).T



        C_lgcv2body = np.array([[np.cos(psi)*np.cos(theta), np.sin(psi)*np.cos(theta), -np.sin(theta)],
                      [np.cos(psi)*np.sin(theta)*np.sin(phi) - np.sin(psi)*np.cos(phi), np.sin(psi)*np.sin(theta)*np.sin(phi) + np.cos(psi)*np.cos(phi), np.cos(theta)*np.sin(phi)],
                    [np.cos(psi)*np.sin(theta)*np.cos(phi) + np.sin(psi)*np.sin(phi), np.sin(psi)*np.sin(theta)*np.cos(phi) - np.cos(psi)*np.sin(phi), np.cos(theta)*np.cos(phi)]])
      
        # print(ref_vectors_lgcv)

        y = C_lgcv2body_quat @ ref_vectors_lgcv.T
        # this is correct ^ vector
        # y is vectors y[:,i] is a vector 

        # difference between measured and calculated: actual - observed
        # y is in body frame
        # we make some observation in body frame
        # then we have some reference vectors in lgcv
        # we use the angles we have to find the equivalent body frame measurements

        dy = vector_obs - y.T
        # dy is x y z, x y z, x y z
        # for NLLS formula
        dy = dy.reshape(dy.size)

        
        # # Building H matrix for each vector measurement
        # H is Jacobian derivatives of C_lgcv2body matrix
        H = []

        # for each observed vector measurement body lgcv frame
        # make H for each sensor: 
        # TODO fix up this
        for m in vector_obs:
            H11 = 2*(qw*m[0] + qz*m[1] - qy*m[2])
            H12 = 2*(qx*m[0] + qy*m[1] + qz*m[2])
            H13 = 2*(- qy*m[0] + qx*m[1] - qw*m[2])
            H14 = 2*(- qz*m[0] + qw*m[1] + qx*m[2])

            H21 = 2*(-qz*m[0] + qw*m[1] + qx*m[2])
            H22 = 2*(qy*m[0] - qx*m[1] + qw*m[2])
            H23 = 2*(qx*m[0] + qy*m[1] + qz*m[2])
            H24 = 2*(-qw*m[0] - qz*m[1] + qy*m[2])

            H31 = 2*(qy*m[0] - qx*m[1] + qw*m[2])
            H32 = 2*(qz*m[0] - qw*m[1] - qx*m[2])
            H33 = 2*(qw*m[0] + qz*m[1] - qy*m[2])
            H34 = 2*(qx*m[0] + qy*m[1] + qz*m[2])
 
            H.append(np.array([
                [H11, H12, H13, H14],
                [H21, H22, H23, H24],
                [H31, H32, H33, H34]]))
    
        # wxyz = 1234

        H = np.array(H)

        # Reshaping to be 2 dimensional
        H = H.reshape(vector_obs.size,4)

       # Equal weights matrix for 
        W = np.eye(vector_obs.size)
        # W[0,0] = 0.9
        # W[1,1] = 0.9
        # W[2,2] = 0.9
     
        # W = np.array([
        #     [0.2,0,0,0,0,0,0,0,0,0,0,0],
        #     [0,0.2,0,0,0,0,0,0,0,0,0,0],
        #     [0,0,0.2,0,0,0,0,0,0,0,0,0],
        #     [0,0,0,1,0,0,0,0,0,0,0,0],
        #     [0,0,0,0,1,0,0,0,0,0,0,0],
        #     [0,0,0,0,0,1,0,0,0,0,0,0],
        #     [0,0,0,0,0,0,1,0,0,0,0,0],
        #     [0,0,0,0,0,0,0,1,0,0,0,0],
        #     [0,0,0,0,0,0,0,0,1,0,0,0],
        #     [0,0,0,0,0,0,0,0,0,1,0,0],
        #     [0,0,0,0,0,0,0,0,0,0,1,0],
        #     [0,0,0,0,0,0,0,0,0,0,0,1]
        # ])

       
    
        # Dilution of Precision     
        pdop = np.sqrt(np.trace( np.linalg.inv( H.T @ H )))

        # NLLS calculation for change in attitude for next iterations
        datt = np.linalg.inv(H.T @ W @ H) @ H.T @ W @ dy
        print(datt)
        # datt is quaternion wxyz

        # print("datt",quat2euler(datt))
        # print("datt quaternion",(datt))

        # this like needs to be differentiated or something
        # for some reason datt is kind of stabilising but it should be going to zero


        if np.sum(np.abs(datt) > tol):
           
            # print("should be same", (euler2quat(quat2euler(att_opt))))
            # print("att_opt", quat2euler(att_opt))

            # can't just add like this!! because i'm making my quaternion longer
            # scaling so datt is always norm of 1 no wonder that doesn't work
            # need to multiply
            print(att_opt)
            att_opt = att_opt + datt.T
            print(att_opt)
            print("START")
            # print(quat2euler(att_opt))
            # print(quat2euler(datt/np.linalg.norm(datt)))
            # # print(quat2euler(datt/np.linalg.norm(datt)))
            # att_opt = quaternion_multiply( datt/np.linalg.norm(datt), att_opt)
            # att_opt = att_opt/np.linalg.norm(att_opt)

            # print(quat2euler(att_opt))


            # print(quat2euler(att_opt))


            # print(quat2euler(datt))
            datts.append(np.abs(np.linalg.norm(datt)))
            # print(att_opt)
            i +=1
            att_store = np.append(att_store,[quat2euler(att_opt)], axis = 0)
        else:
            break
    
    print()
    print(f"Number of iterations to converge: {i}")
    print("PDOP:",pdop)
    return att_opt, pdop, att_store, i, datts

def quaternion_multiply(quaternion1, quaternion0):
    w0, x0, y0, z0 = quaternion0
    w1, x1, y1, z1 = quaternion1
    return np.array([-x1*x0 - y1*y0 - z1*z0 + w1*w0,
                        x1*w0 + y1*z0 - z1*y0 + w1*x0,
                        -x1*z0 + y1*w0 + z1*x0 + w1*y0,
                        x1*y0 - y1*x0 + z1*w0 + w1*z0])

def main():

    att_init = (np.array([10,32,-45]))

    # att_data = np.array([[10,31,-41],[9,30,-47],[14,32,-44]])

    # say we have 4 sensors
    # this is in body frame
    vector_obs = np.array([[ 0.3,  0.5,  0.3],[-0.08, -0.01,   0.5], [ 0.4, -0.9, -0.4],[ 0.4, -0.9, -0.3]])
    vector_obs = np.array([[-0.3, 0.3, -0.5], [-0.4, -0.2, -0.01], [-0.2, 0.5, 0.9], [-0.2, -0.5, 0.9]])
    # vector_obs = np.array([[-0.0879, 0.5242, -0.6383],[-0.3319, 0.3281, 0.2055], [0.8465, 0.0540, 0.6485]])
    vector_obs = np.array([[ 0.48,  0.58  ,0.34],[-0.08 ,-0.09  , 0.49],[ 0.41 ,-0.90, -0.37]])
    # vector_obs = np.array([[-0.2, 1, -0.3],[-0.4, 0.3, 0.4], [0.5, -0.07, 0.3]])
    # [[ 0.48208398  0.58341093  0.34238388]
    #  [-0.08427267 -0.0996557   0.49291668]
    #  [ 0.41381351 -0.90926659 -0.37681912]
    #  [ 0.41381351 -0.90926659 -0.37681912]]
    # this is in lgcv
    ref_vectors_lgcv = np.array([[0.2,0.7,-0.4],[0.1,0.3,0.4],[0.7,-0.8,0.1] ])
    # needs to be in radians
    att_opt, pdop, att_store, i, datts = nlls_quaternion_weights(vector_obs, ref_vectors_lgcv, att_init)

    print("Final attitude estimation", quat2euler(att_opt))
    print(att_store.shape)
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



# justifying nlls



"""
Tasks TODO

add quaternions
run with jacobian quaternions

get to converge
test with different weights

figure out how this interacts: what data do i need
what actually makes sense
make some sample data
integrate with position

make nice code (obj oriented)



"""

if __name__ == '__main__':
    main()
