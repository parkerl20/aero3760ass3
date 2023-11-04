import matplotlib.pyplot as plt
import numpy as np
from pos_sun import *
from attitude_propagation import *
from attitude_transforms import *
'''
NLLS for Attitude Determination
'''


def nlls_quaternion_weights(vector_obs, ref_vectors_lgcv, att_init, dy_flag):
    """
    Finds Non Linear Least Squares estimate of attitude in quaternions , DOP
    Inputs:
        
        vector_obs: observation vectors from sensors (sun sensor, star tracker (stars) )
        ref_vectors: reference vectors to objects detected by sensors (sun, stars)
        att_init: initial guess of attitude of sat in Quaternion
        dy_flag: either 1 or -1 depending on convergence - to multiply against dy

    Outputs:
        pos_opt: final estimation as quaternion
        dops: dop for quaternion components
        pos_store: All the computed attitude estimates over the iterations
        i: Number of iterations for convergence
    """
    
    # Initialise arrays 
    # store as euler for easy interpretation - shouldn't affect calculation as this isn't being used
    att_store = np.array([(quat2euler(att_init))])

    # Initial iteration params
    i = 0
    max_iter = 100
    tol =  1e-4
    datt = 100

    # initial guess in quaternion
    att_opt = att_init

    datts = []
    dops =[]

    # until converges
    while (i < max_iter) and (np.sum(np.abs(datt) > tol)):
        
        # for each angle find calculated  based on previous estimate
        qw = att_opt[0]
        qx = att_opt[1]
        qy = att_opt[2]
        qz = att_opt[3]



        # convert LGCV to body frame using quaternion ( transposed to fix maths )
        C_lgcv2body_quat = np.array([[qw**2 + qx**2 - qy**2 - qz**2, 2*(qx*qy - qw*qz), 2*(qx*qz + qw*qy)],
                                [2*(qx*qy + qw*qz), qw**2 - qx**2 + qy**2 - qz**2, 2*(qy*qz - qw*qx)],
                                [2*(qx*qz - qw*qy), 2*(qy*qz + qw*qx), qw**2 - qx**2 - qy**2 + qz**2]]).T



        # Find reference vectors in body frame - expected obs
        y = C_lgcv2body_quat @ ref_vectors_lgcv.T
        # y is vectors y[:,i] is a vector 

        # difference between measured and calculated: actual - expected
        dy = dy_flag*(vector_obs - y.T)

       
        # dy is x y z, x y z, x y z
        # for NLLS formula
        dy = dy.reshape(dy.size)
        
        # # Building H matrix for each vector measurement: Jacobian matrix derivatives for C_lgcv2body matrix
        H = []

        # for each observed vector measurement body frame
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
        # Sun sensor accuracy 0.23
        # star tracker 1
        W[0,0] = 0.23
        W[1,1] = 0.23
        W[2,2] = 0.23

        W[3,3] = 0.23
        W[4,4] = 0.23
        W[5,5] = 0.23

            
        # Dilution of Precision for w,x,y,z of quaternion     
        dops =  [np.sqrt(np.linalg.inv( H.T @ H)[0,0]), np.sqrt(np.linalg.inv( H.T @ H)[1,1]), np.sqrt(np.linalg.inv( H.T @ H)[2,2]), np.sqrt(np.linalg.inv( H.T @ H)[3,3])]

        # NLLS calculation for change in attitude for next iterations
        datt = np.linalg.inv(H.T @ W @ H) @ H.T @ W @ dy
        # datt is for quaternion wxyz

        # if not converged
        if np.sum(np.abs(datt) > tol):
            
            # update attitude - datt is not a quaternion itself
            att_opt = att_opt  + datt.T
            # att_opt = quaternion_multiply(datt/np.linalg.norm(datt), att_opt)
     
            # normalise quaternion
            att_opt = att_opt/np.linalg.norm(att_opt)
            att_opt = att_opt/np.linalg.norm(att_opt)
            
            datts.append(np.abs(np.linalg.norm(datt)))
            i +=1
            # store attitude as euler angles for interpretation
            att_store = np.append(att_store,[quat2euler(att_opt)], axis = 0)

        else:
            break
    
    print(f"\nNumber of iterations to converge: {i}")
    if i == 100:
        not_converged = 1
    else:
        not_converged = 0
    # print("PDOP:",dops)
    return att_opt/np.linalg.norm(att_opt), dops, att_store, i, datts, not_converged


def get_reference_vectors(num_stars, num_sensors, time):
    # number of sensor measurements at one time interval (redundancy number)

    stars, sun = findStarSunPositions(num_stars*num_sensors,time)

    # try with just one
    ref_vectors = []
    for i in range(num_sensors):
        ref_vectors.append(sun)
    for i in range(num_sensors*num_stars):
        ref_vectors.append(stars[:,i])
    ref_vectors = np.array(ref_vectors)
    # ref_vectors = np.array([sun,sun, stars[:,0], stars[:,1],stars[:,2], stars[:,3], stars[:,4], 
    #                         stars[:,5], stars[:,6],stars[:,7], stars[:,8], stars[:,9]])
    # duplicate for the number of measurements we have from each timestep
    # ref_vectors = np.repeat(ref_vectors, n, axis = 0)
    # ref_vectors = np.tile(ref_vectors, n)
    return ref_vectors

# def generate_sun_star_data(t):




def main():

    time_start = 0
    # end of simulation
    one_period = 10000# 60000


    initial_euler = np.array([10, 30, -45])
    initial_quat = euler2quat(initial_euler)
    
    # find real attitude over period based on nadir pointing satellite
    attitude = solve_ivp(propagate_attitude, [time_start, one_period], initial_quat, method='RK45', max_step=10) 
    attitude_poses = attitude.y[:, :]
    times = attitude.t

    # Initialise array of true euler angles over propagation for comparison
    eulers = []
    for i in range(len(times)):
        q = attitude_poses[:, i]
        euler = quat2euler(q)
        eulers.append(euler)
    eulers = np.array(eulers)

    # Initialise arrays to store data from nlls
    final_estimates = np.zeros((len(times), 3))
    iterations = []
    dops_time = np.zeros((len(times), 4))


    for t in range(len(times)):
        # generate reference vectors for sun 20 stars and 2 of each sensor
        num_stars = 20
        num_sensors = 2
        freq = 5 # sensors have 5Hz (5 measurements / second)
        # number of measurements taken each second
        n = num_sensors*freq
        ref_vectors = get_reference_vectors(num_stars, n,times[t])

        # generate euler angles for sun and star data 
        # errors from sensor specs in degrees
        euler_sun_sensor_errors = np.random.normal(eulers[t-1], 0.1 , (n,3)) 
        euler_star_tracker_errors =np.random.normal(eulers[t-1], [0.0023, 0.03, 0.03 ], (n,3)) # yaw, pitch, roll (zyx)

        # quat_sun_sensor_errors = np.random.normal(attitude_poses[:,t], 0.0001, (n,4)) # 0.01306
        # quat_star_tracker_errors =np.random.normal(attitude_poses[:,t], 0.0001, (n,4)) 
        
        
        observed_vectors = np.zeros((n*(num_stars+1), 3))

        # generate data for sun sensor
        # sun sensor readings are the same 
        count = 0
        # refs are in order sun sun star star
        for angles in euler_sun_sensor_errors:

            quat = euler2quat(angles)
            
            qw = quat[0]
            qx = quat[1]
            qy = quat[2]
            qz = quat[3]

            C_lgcv2body_quat = np.array([[qw**2 + qx**2 - qy**2 - qz**2, 2*(qx*qy - qw*qz), 2*(qx*qz + qw*qy)],
                                [2*(qx*qy + qw*qz), qw**2 - qx**2 + qy**2 - qz**2, 2*(qy*qz - qw*qx)],
                                [2*(qx*qz - qw*qy), 2*(qy*qz + qw*qx), qw**2 - qx**2 - qy**2 + qz**2]]).T



            observed_vectors[count] = (C_lgcv2body_quat @ ref_vectors[0])
            count = count + 1

        # generate data for star_tracker sensor
        # count = 2
        for angles in euler_star_tracker_errors:

            quat = euler2quat(angles)
            qw = quat[0]
            qx = quat[1]
            qy = quat[2]
            qz = quat[3]

            C_lgcv2body_quat = np.array([[qw**2 + qx**2 - qy**2 - qz**2, 2*(qx*qy - qw*qz), 2*(qx*qz + qw*qy)],
                                [2*(qx*qy + qw*qz), qw**2 - qx**2 + qy**2 - qz**2, 2*(qy*qz - qw*qx)],
                                [2*(qx*qz - qw*qy), 2*(qy*qz + qw*qx), qw**2 - qx**2 - qy**2 + qz**2]]).T

            observed_vectors[count:  num_stars+count] = ((C_lgcv2body_quat @ ref_vectors[count:num_stars+count].T).T)
            count = count + num_stars
        observed_vectors = np.array(observed_vectors)

        # initial guess in quaternion
        att_init = euler2quat(np.array([10, 30, -45]))  

        # used to stitch together different converging solutions
        dy_flag = 1
        att_opt, dops, att_store, i, datts, not_converged = nlls_quaternion_weights(observed_vectors, ref_vectors, att_init, dy_flag)
        if not_converged:
            dy_flag = -1
            att_opt, dops, att_store, i, datts, not_converged = nlls_quaternion_weights(observed_vectors, ref_vectors, att_init, dy_flag)
        dops_euler = euler2quat(dops)
        final_estimates[t-1] = (quat2euler(att_opt))
        iterations.append(i)
        # next guess is final estimate
        att_init = att_opt
        
        print("Initial guess", quat2euler(att_init))
        print("Actual attitude", eulers[t-1])
        dops_time[t-1] = dops_euler
        # dops times vector error for actual mapping error
        # errors are 1° maybe
        # print("DOP for euler params Z Y X:", quat2euler(dops))

        # get actual error for 
        print("Final attitude estimation", quat2euler(att_opt))

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
    
    plot_attitude_propagation(times, eulers-final_estimates)
    plot_attitude_propagation(times, final_estimates)

    plt.figure()
    plt.plot(range(len(iterations)), iterations, label = "iterations to converge")
    plt.legend()
    plt.ylabel("Iterations to Converge")
    plt.xlabel("Timesteps")
    plt.show()



if __name__ == '__main__':
    main()


