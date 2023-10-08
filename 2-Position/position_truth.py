import matplotlib.pyplot as plt
import numpy as np

'''

A script that produces an orbit (analytically) that
acts as our source of "truth" for position determination via Kalman Filter.

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
    X = X in ECEF Coordinate
    Y = Y in ECEF Coordinates
    Z = Z in ECEF Coordinates
    '''

    '''
    Produce an ECEF plot using functions from assignment 1 and 2
    '''
    
    
    