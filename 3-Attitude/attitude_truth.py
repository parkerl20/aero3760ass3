import matplotlib.pyplot as plt
import numpy as np

'''

A script that produces an example attitude manoeuvre (analytically) that
acts as our source of "truth" for attitude determination via Kalman Filter.

'''

def generate_truth_orbit():
    
    print("Generating truth orbit...")
    
    # Define the start and end points and the number of values
    start = 0
    end = 2500
    num_values = 20000

    # Generate the linspace
    time_steps = np.linspace(start, end, num_values)
    time_steps = time_steps.reshape(-1, 1)
        
    '''
    Theta_1 = pitch
    Theta_2 = yaw
    Theta_3 = roll
    '''
    
    theta_1 = 10 * np.e**(-0.02 * time_steps) * np.sin(11.4592 * time_steps)
    theta_2 = 10 * np.e**(-0.02 * time_steps) * np.cos(11.4592 * time_steps)
    theta_3 = np.sin(2.8648 * time_steps)
    
    # Create a figure and axes
    fig1, ax = plt.subplots(figsize=(8, 8))
    ax.plot(theta_2, theta_1, label='Archimedean Spiral', color='b')
    ax.set_title('Satellite Manoeuvre')
    ax.set_xlabel('Yaw (degrees)')
    ax.set_ylabel('Pitch (degrees)')
    ax.legend()
    ax.grid(True)
    ax.axis('equal')
    
    # Create a figure with three subplots in one row and three columns
    fig2, axes = plt.subplots(3, 1, figsize=(9, 3))

    axes[0].plot(time_steps, theta_2, label='Yaw (deg)', color='b')
    axes[0].set_title('Plot 1')
    axes[0].legend()

    axes[1].plot(time_steps, theta_1, label='Pitch (deg)', color='g')
    axes[1].set_title('Plot 2')
    axes[1].legend()

    axes[2].plot(time_steps, theta_3, label='Roll (deg)', color='r')
    axes[2].set_title('Plot 3')
    axes[2].legend()

    # Add a title for the entire figure
    fig2.suptitle('Three Subplots')

    plt.show()
    
    return 0