import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from . import attitude_transforms

'''

NLLS for Static Attitude Determination

'''


def propagate_attitude(t, q0):

    '''
        Propagates the attitude from a 
        matrix determined by quaternion angle
        velocities.

        q - np.array (size 4x1)
        p, q, r - angular velocities in Euler
                  angle representation
    '''

    p = 0.0016
    q = 0.0016
    r = 0.0016

    q_matrix = 0.5 * np.array([
        [0, -p, -q, -r],
        [p, 0, r, -q],
        [q, -r, 0, p],
        [r, q, -p, 0]])
    
    q_dot = 0.5 * q_matrix @ q0

    print("Time:", t, "q:", q0)
    
    return q_dot


def plot_attitude_propagation(time, eulers, fig_name):

    z = [sub_array[0] for sub_array in eulers]
    y = [sub_array[1] for sub_array in eulers]
    x = [sub_array[2] for sub_array in eulers]

    # Create a figure and a set of subplots
    fig, axs = plt.subplots(3, 1, figsize=(8, 6))

    # Plot the data
    axs[0].plot(time, x, label='x(t)', color='b')
    axs[1].plot(time, y, label='y(t)', color='r')
    axs[2].plot(time, z, label='z(t)', color='g')

    # Set titles and labels
    axs[0].set_title('X over Time')
    axs[1].set_title('Y over Time')
    axs[2].set_title('Z over Time')

    # Set y-axis label
    axs[0].set_ylabel('Roll (deg)')
    axs[1].set_ylabel('Pitch (deg)')
    axs[2].set_ylabel('Yaw (deg)')

    for ax in axs:
        ax.set_xlabel('Time')
        ax.legend()

    # Adjust layout and display the plot
    plt.tight_layout()
    plt.savefig(fig_name)

    return 0


def main():

    time_start = 0
    one_period = 10000

    initial_euler = [10, 30, -45]
    initial_quat = attitude_transforms.euler2quat(initial_euler)

    attitude = solve_ivp(propagate_attitude, [time_start, one_period], initial_quat, method='RK45', max_step=10) 
    attitude_poses = attitude.y[:, :]
    times = attitude.t

    # Create a NumPy array of zeros
    eulers = np.zeros(len(times), dtype=object)

    # Resize each index to contain a NumPy array of 1x3
    for i in range(len(times)):
        eulers[i] = np.zeros((1, 3))

    for i in range(0, len(eulers)):
        q = attitude_poses[:, i]
        euler = attitude_transforms.quat2euler(q)
        eulers[i] = euler
        # eulers[i] = attitude_transforms.quat2euler(attitude.y[:, i])

    print("\n\n")
    plot_attitude_propagation(times, eulers, "real_att")

    return 0



if __name__ == "__main__":
    main()



# Attitude Truth
'''
Pixel Location Determination
- FOV
- Camera Matrix
- Star Magnitude
Star Tracker Noise Characterisation
- Photon Noise
- Electron Noise
- Analog-to-Digital Noise
Image Processing
- Noise Reduction
- Thresholding and star finding algorithm
- Centroiding (centre of mass method is easier)
Catalog Formation
- Reference catalog
- Map stars onto Earth planet
'''
