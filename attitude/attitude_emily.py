import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial.transform import Rotation
from pos_sun import *
from attitude_propagation import *
'''
NLLS for Attitude Determination
'''

# from ECI 
# get nadir  to lat long of earth
# need to figure out rotation between body frame and vector to lat long of nadir
# ECI 2 ECEF 2 ENU
# TODO
# integrate with position

# point sensors in different directions - different stars?
# account for sun going into dark - time of day
# satellite needs always point nadir - how does this translate ??
# scannign system
# change attitude over time - to do with facing towards earth
# weightings
# for over time - use previous output as initial gues


# yaw is perturbation
# pitch change 
# roll perturbation


# read eci from csv
# pos_eci = np.loadtxt('satellite_position.csv', delimiter=',')    


# convert to ecef
# use ot object
# convert to enu


# 1/1/2023 is epoch







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
    
    # normalise
    quat = np.array([ qw,qx, qy, qz])
    quat = quat/np.linalg.norm(quat)
    return quat


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

    Finds Non Linear Least Squares estimate of attitude in quaternions , DOP
    Inputs:
        
        vector_obs: observation vectors from sensors
        ref_vectors: reference vectors to objects detected by sensors (sun, stars)
        att_init: initial guess of attitude of sat in Euler angles

    Outputs:
        pos_opt: final estimation as quaternion
        dops: dop for quaternion components
        pos_store: All the computed attitude estimates over the iterations
        i: Number of iterations for convergence
    """
    
    # Initialise arrays 
    # store as euler for easy interpretation
    att_store = np.array([(quat2euler(att_init))])

    # Initial iteration params
    i = 0
    max_iter = 100
    tol =  1e-4
    datt = 100

    # initial guess in quaternion
    att_opt = att_init
    # att_opt = att_init

    datts = []
    dops =[]

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



        C_lgcv2body_euler = np.array([[np.cos(psi)*np.cos(theta), np.sin(psi)*np.cos(theta), -np.sin(theta)],
                      [np.cos(psi)*np.sin(theta)*np.sin(phi) - np.sin(psi)*np.cos(phi), np.sin(psi)*np.sin(theta)*np.sin(phi) + np.cos(psi)*np.cos(phi), np.cos(theta)*np.sin(phi)],
                    [np.cos(psi)*np.sin(theta)*np.cos(phi) + np.sin(psi)*np.sin(phi), np.sin(psi)*np.sin(theta)*np.cos(phi) - np.cos(psi)*np.sin(phi), np.cos(theta)*np.cos(phi)]])
      
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

        # Weights matrix based on sensor accuracy
        W = np.eye(vector_obs.size)
        # 36 vector - we have 12 sensors each with x y z accuracy
        # cannot weight euler angles - can only weight sensor accuracy
        # Sun sensor accuracy 0.3
        # star tracker 1
        W[0,0] = 0.3
        W[1,1] = 0.3
        W[2,2] = 0.3

        W[3,3] = 0.3
        W[4,4] = 0.3
        W[5,5] = 0.3

        # # Star tracker z component less accurate
        # W[5,5] = 0.14
        # W[8,8] = 0.14
        # W[11,11] = 0.14
        # W[14,14] = 0.14
        # W[17,17] = 0.14
            
        # Dilution of Precision for w,x,y,z of quaternion     
        dops =  [np.sqrt(np.linalg.inv( H.T @ H)[0,0]), np.sqrt(np.linalg.inv( H.T @ H)[1,1]), np.sqrt(np.linalg.inv( H.T @ H)[2,2]), np.sqrt(np.linalg.inv( H.T @ H)[3,3])]

        # NLLS calculation for change in attitude for next iterations
        datt = np.linalg.inv(H.T @ W @ H) @ H.T @ W @ dy
        # datt is quaternion wxyz

         
        if np.sum(np.abs(datt) > tol):
            # print("START")
            # print(quat2euler(att_opt))
            # print("datt", datt)
            att_opt = att_opt + datt.T
            # print(quat2euler(att_opt))
            att_opt = att_opt/np.linalg.norm(att_opt)
            # print(quat2euler(att_opt))
            datts.append(np.abs(np.linalg.norm(datt)))
            i +=1
         
        
            att_store = np.append(att_store,[quat2euler(att_opt)], axis = 0)
        else:
            att_opt = att_opt/np.linalg.norm(att_opt)
            break
    
    print()
    print(f"Number of iterations to converge: {i}")
    if i == 100:
        not_converged = 1
    else:
        not_converged = 0
    # print("PDOP:",dops)
    return att_opt/np.linalg.norm(att_opt), dops, att_store, i, datts, not_converged

def quaternion_multiply(quaternion1, quaternion0):
    w0, x0, y0, z0 = quaternion0
    w1, x1, y1, z1 = quaternion1
    return np.array([-x1*x0 - y1*y0 - z1*z0 + w1*w0,
                        x1*w0 + y1*z0 - z1*y0 + w1*x0,
                        -x1*z0 + y1*w0 + z1*x0 + w1*y0,
                        x1*y0 - y1*x0 + z1*w0 + w1*z0])


def main():

    # att_init = (np.array([10,32,-45]))
    # initial guess
    att_init = euler2quat(np.array([0,0,0]))
    real_attitude = [10, 30, -45]
    # real_attitude = [ 12.85378879 , 47.59094252,-30.52834642]

    # receive sensor data (stars position and sun position) based on 5 stars detected by star trackers
    # for 22 october

    # time = np.arange(1,10,1)
  

    time_start = 0
    one_period = 300# 60000

    initial_euler = [10, 30, -45]
    initial_quat = attitude_transforms.euler2quat(initial_euler)

    attitude = solve_ivp(propagate_attitude, [time_start, one_period], initial_quat, method='RK45', max_step=10) 
    attitude_poses = attitude.y[:, :]
    times = attitude.t

    # Create a NumPy array of zeros
    eulers = []
    # eulers = np.zeros(len(times), dtype=object)

    # Resize each index to contain a NumPy array of 1x3
    # for i in range(len(times)):
    #     eulers[i] = np.zeros((1, 3))

    # for i in range(0, len(eulers)):
    for i in range(len(times)):
        q = attitude_poses[:, i]
        euler = attitude_transforms.quat2euler(q)
        # print(euler.shape, type(euler))
        # eulers[i] = euler
        eulers.append(euler)
        # print(eulers[i])
        # eulers[i] = attitude_transforms.quat2euler(attitude.y[:, i])
    eulers = np.array(eulers)


    final_estimates = np.zeros((len(times), 3))
    iterations = []
    dops_time = np.zeros((len(times), 4))





    # real_attitude = [10, 50,-45]
    # real_attitude = [15,-45,30]
    # real_attitude = [ 16,   53, -23]

    print((euler2quat(real_attitude)))
    # [ 0.87551721 -0.38865382  0.22162565  0.1824826 ]
    #nice quat ^^

    # [ 0.343162    0.12721205  0.24932435 -0.13266423]
    # [0.14398904 0.47897751 0.83864741 0.21568086]
    
    # print(np.linalg.norm(euler2quat([15,-45,30])))
    # print(np.linalg.norm(euler2quat([10, 50,-45])))
    # print(np.linalg.norm(euler2quat([10, 30,-45])))


    print(quat2euler(euler2quat([70,30,-45])))

    # init guess
    # att_init = euler2quat([60,30,-40])
    

    # eulers = np.array(eulers)
    # # eulers = np.random.normal(0, [100,100,100] , (50,3))
    # print(eulers)

    
    for t in range(len(times)):
        num_stars = 10
        # number of sensor measurements at one time interval
        n = 2

        stars, sun = findStarSunPositions(num_stars*n,times[t])

        # try with just one
        ref_vectors = []
        for i in range(n):
            ref_vectors.append(sun)
        for i in range(n*num_stars):
            ref_vectors.append(stars[:,i])
        ref_vectors = np.array(ref_vectors)
        # ref_vectors = np.array([sun,sun, stars[:,0], stars[:,1],stars[:,2], stars[:,3], stars[:,4], 
        #                         stars[:,5], stars[:,6],stars[:,7], stars[:,8], stars[:,9]])
        # duplicate for the number of measurements we have from each timestep
        # ref_vectors = np.repeat(ref_vectors, n, axis = 0)
        # ref_vectors = np.tile(ref_vectors, n)
        print(ref_vectors)



        # in euler angles
        # sun sensor random error 47 arc seconds = 0.01306 degrees
        # star tracker random error XYZ = 50, 50, 350 microradians = 0.00286, 0.00286, 0.02005 degrees
        # print("time", type(t), times)
        # print()
        # generate in 100 euler angles

        # print("EULERS", type((eulers[t])), eulers[t-1], t)
        print(len(eulers))
        print(len(times))
        print(times[-1])


        euler_sun_sensor_errors = np.random.normal(eulers[t-1], 0.01306 , (n,3)) # 0.01306
        euler_star_tracker_errors =np.random.normal(eulers[t-1], [0.00286,  0.02005, 0.00286,], (n,3)) 


        print("SENSORS",euler_sun_sensor_errors)
        print(euler_star_tracker_errors)
        print("REAL DATA", eulers[t-1])


        observed_vectors = np.zeros((n*(num_stars+1), 3))

        # generate data for sun sensor
        # sun sensor readings are the same 
        count = 0
        # refs are in order sun sun star star
        for angles in euler_sun_sensor_errors:

            quat = euler2quat(angles)
            

            # psi = np.deg2rad(angles[0])
            # theta = np.deg2rad(angles[1])
            # phi = np.deg2rad(angles[2])

            # C_lgcv2body = np.array([[np.cos(psi)*np.cos(theta), np.sin(psi)*np.cos(theta), -np.sin(theta)],
            #             [np.cos(psi)*np.sin(theta)*np.sin(phi) - np.sin(psi)*np.cos(phi), np.sin(psi)*np.sin(theta)*np.sin(phi) + np.cos(psi)*np.cos(phi), np.cos(theta)*np.sin(phi)],
            #             [np.cos(psi)*np.sin(theta)*np.cos(phi) + np.sin(psi)*np.sin(phi), np.sin(psi)*np.sin(theta)*np.cos(phi) - np.cos(psi)*np.sin(phi), np.cos(theta)*np.cos(phi)]])
            
            qw = quat[0]
            qx = quat[1]
            qy = quat[2]
            qz = quat[3]

            C_lgcv2body_quat = np.array([[qw**2 + qx**2 - qy**2 - qz**2, 2*(qx*qy - qw*qz), 2*(qx*qz + qw*qy)],
                                [2*(qx*qy + qw*qz), qw**2 - qx**2 + qy**2 - qz**2, 2*(qy*qz - qw*qx)],
                                [2*(qx*qz - qw*qy), 2*(qy*qz + qw*qx), qw**2 - qx**2 - qy**2 + qz**2]]).T



            observed_vectors[count] = (C_lgcv2body_quat @ ref_vectors[0])
            count = count + 1
        print(observed_vectors)
        # generate data for star_tracker sensor


        # count = 2
        for angles in euler_star_tracker_errors:

            psi = np.deg2rad(angles[0])
            theta = np.deg2rad(angles[1])
            phi = np.deg2rad(angles[2])

            quat = euler2quat(angles)
            qw = quat[0]
            qx = quat[1]
            qy = quat[2]
            qz = quat[3]

            C_lgcv2body_quat = np.array([[qw**2 + qx**2 - qy**2 - qz**2, 2*(qx*qy - qw*qz), 2*(qx*qz + qw*qy)],
                                [2*(qx*qy + qw*qz), qw**2 - qx**2 + qy**2 - qz**2, 2*(qy*qz - qw*qx)],
                                [2*(qx*qz - qw*qy), 2*(qy*qz + qw*qx), qw**2 - qx**2 - qy**2 + qz**2]]).T
            print("count", count)

            observed_vectors[count:  num_stars+count] = ((C_lgcv2body_quat @ ref_vectors[count:num_stars+count].T).T)
            count = count + num_stars
        observed_vectors = np.array(observed_vectors)

        for vector in observed_vectors:
            print(C_lgcv2body_quat.T @ vector)

        # print(ref_vectors)

        # print("OBSERVED VECTORS", observed_vectors)
        # print("REF VECTORS", ref_vectors)
   
        # needs to be in radians
        att_opt, dops, att_store, i, datts, not_converged = nlls_quaternion_weights(observed_vectors, ref_vectors, att_init)
        dops_euler = euler2quat(dops)
        final_estimates[t-1] = (quat2euler(att_opt))
        iterations.append(i)
        # next guess is final estimate
        att_init = att_opt
        
        print("INIT ATT", quat2euler(att_init))
        print("ATT ACTUALLY", eulers[t-1])
        dops_time[t-1] = dops_euler

        # dops times vector error for actual mapping error
        # errors are 1° maybe
        # print("DOP for euler params Z Y X:", quat2euler(dops))

        # get actual error for 
        

        print("Final attitude estimation Z Y X:", quat2euler(att_opt))
        print("DOP for quaternion params W X Y Z:", dops)


    print(final_estimates)
    
    plt.figure()
    plt.plot(range(len(att_store)), np.full((len(att_store)),eulers[-1,0]), label = "True Yaw")
    plt.plot(range(len(att_store)), np.full((len(att_store)),eulers[-1,1]), label = "True Pitch")
    plt.plot(range(len(att_store)), np.full((len(att_store)),eulers[-1,2]), label = "True Roll")
    plt.plot(range(len(att_store)), att_store[:,0], label = "Yaw")
    plt.plot(range(len(att_store)), att_store[:,1], label = "Pitch")
    plt.plot(range(len(att_store)), att_store[:,2], label = "Roll")
    plt.legend(loc = 'upper right')
    plt.xlabel("Iterations")
    plt.ylabel("Value (degrees)")
    plt.savefig("values")

    plt.figure()
    plt.plot(range(i), datts, label = "Quaternion update value")
    plt.legend()
    plt.xlabel("Iterations")
    plt.ylabel("Value")
    plt.savefig("converge")

    # plot_attitude_propagation(times, eulers)
    # eulers_long = np.random.normal(eulers, 0.02, (len(eulers),3))

    
    plot_attitude_propagation(times, eulers-final_estimates)


    # plt.figure()
    # plt.plot(range(len(final_estimates)), final_estimates[:,0] - real_attitude[0], label = "Yaw")
    # plt.plot(range(len(final_estimates)), final_estimates[:,1] - real_attitude[1], label = "Pitch")
    # plt.plot(range(len(final_estimates)), final_estimates[:,2] - real_attitude[2], label = "Roll")
    # plt.legend()
    # plt.ylabel("Final estimates deviation from true")
    # plt.xlabel("Time")
    # plt.savefig("overtime")

    # plt.figure()
    # plt.plot(range(len(dops_time)), dops_time[:,0], label = "W")
    # plt.plot(range(len(dops_time)), dops_time[:,1], label = "X")
    # plt.plot(range(len(dops_time)), dops_time[:,2], label = "Y")
    # plt.plot(range(len(dops_time)), dops_time[:,3], label = "Z")
    # plt.legend()
    # plt.ylabel("DOPs")
    # plt.xlabel("Timesteps")
    # plt.show()

    plt.figure()
    plt.plot(range(len(dops_time)), dops_time[:,0], label = "Yaw")
    plt.plot(range(len(dops_time)), dops_time[:,1], label = "Pitch")
    plt.plot(range(len(dops_time)), dops_time[:,2], label = "Roll")
    plt.legend()
    plt.ylabel("DOPs")
    plt.xlabel("Timesteps")
    plt.show()

    plt.figure()
    plt.plot(range(len(iterations)), iterations, label = "iterations to converge")
    plt.legend()
    plt.ylabel("Iterations to Converge")
    plt.xlabel("Timesteps (over 300s)")
    plt.show()



def testing():

    # att_data = np.array([[10,31,-41],[9,30,-47],[14,32,-44]])

    # say we have 4 sensors this is sensor data
    # this is in body frame
    vector_obs = np.array([[ 0.3,  0.5,  0.3],[-0.08, -0.01,   0.5], [ 0.4, -0.9, -0.4],[ 0.4, -0.9, -0.3]])
    vector_obs = np.array([[-0.3, 0.3, -0.5], [-0.4, -0.2, -0.01], [-0.2, 0.5, 0.9], [-0.2, -0.5, 0.9]])
    # vector_obs = np.array([[-0.0879, 0.5242, -0.6383],[-0.3319, 0.3281, 0.2055], [0.8465, 0.0540, 0.6485]])




    # clean data
    # sensor data in the body frame - corresponds to local frame data
    vector_obs = np.array([[ 0.48,  0.58  ,0.34],[-0.08 ,-0.09  , 0.49]])#[ 0.41 ,-0.90, -0.37]])

    # vector_obs = np.array([[-0.2, 1, -0.3],[-0.4, 0.3, 0.4], [0.5, -0.07, 0.3]])
    # [[ 0.48208398  0.58341093  0.34238388]
    #  [-0.08427267 -0.0996557   0.49291668]
    #  [ 0.41381351 -0.90926659 -0.37681912]
    #  [ 0.41381351 -0.90926659 -0.37681912]]

    # this is in lgcv
    # where does this come from
    # link with position
    # 

    psi = np.deg2rad(real_attitude[0])
    theta = np.deg2rad(real_attitude[1])
    phi = np.deg2rad(real_attitude[2])

    C_lgcv2body = np.array([[np.cos(psi)*np.cos(theta), np.sin(psi)*np.cos(theta), -np.sin(theta)],
                      [np.cos(psi)*np.sin(theta)*np.sin(phi) - np.sin(psi)*np.cos(phi), np.sin(psi)*np.sin(theta)*np.sin(phi) + np.cos(psi)*np.cos(phi), np.cos(theta)*np.sin(phi)],
                    [np.cos(psi)*np.sin(theta)*np.cos(phi) + np.sin(psi)*np.sin(phi), np.sin(psi)*np.sin(theta)*np.cos(phi) - np.cos(psi)*np.sin(phi), np.cos(theta)*np.cos(phi)]])
     
    ref_vectors_lgcv = np.array([[0.2,0.7,-0.4],[0.1,0.3,0.4]])#,[0.7,-0.8,0.1] ])
    print("TEST",C_lgcv2body@ref_vectors_lgcv.T )
    print("TEST",C_lgcv2body@ref_vectors_lgcv[0] )







if __name__ == '__main__':
    main()



