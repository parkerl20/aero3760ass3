'''

Static Attitude Determination

'''

'''

Preliminary Design Review Functions

'''

def nlls_orbit_determination(attitude_0,t_span,obs):
    '''
    NLLS_ORBIT_DETERMINATION - uses Non-Linear Least Squares to give a set of orbital parameters
    given a set of observations.

    Input: atttitude_0 (numpy.ndarray) = 3 x 1 vector of roll, pitch, yaw for initial observation
           t_span (numpy.ndarray) = n x discrete times the observations were taken at (1xn np.array)
           obs (numpy.ndarray) = n x observations of the satellite position in LGCV at each time in t_span (3xn np.array)

    Output: attitude_matrix (numpy.ndarray) = 3 x 1 matrix of estimated roll, pitch, yaw of attitude
            
    '''
    
    '''
    
    state accepted tolerance
    determine state vector
    calculate H, Jacobians
    assign W, Weightings
    create iteration loop for deltaY and H calculation
        until tolerance is met
    
    return attitude_matrix
    
    '''
    

    
