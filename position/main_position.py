# -------- Add spacesim to script path
import sys
import os

current_file_path = os.path.abspath(__file__)
parent_directory = os.path.dirname(current_file_path)
grandparent_directory = os.path.dirname(parent_directory)
sys.path.append(grandparent_directory)

# --------- other imports
from spacesim import satellite as sat
from spacesim import estimation as est
from spacesim import constants as const
from spacesim import orbit_determination as od
from spacesim import celestial_body as cb

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np


def transition_matrix_func(r: np.ndarray, dt: float) -> np.ndarray:    
    mu = const.MU_EARTH
    r_i, r_j, r_k = r.ravel()
    r_mag = np.linalg.norm(r)
    
    # print(f"r: {r.ravel()}")
    
    # intermediate matrices
    # [
    #     [A, B],
    #     [C, D]
    # ]
   
    A = np.zeros((3, 3))
    B = np.eye(3)
    D = np.zeros((3, 3))
    
    c_00 = -(mu / r_mag**3) + (3 * mu * r_i**2 / r_mag**5)
    c_01 = 3 * mu * r_i * r_j / r_mag**5
    c_02 = 3 * mu * r_i * r_k / r_mag**5

    c_10 = 3 * mu * r_j * r_i / r_mag**5
    c_11 = -(mu / r_mag**3) + (3 * mu * r_j**2 / r_mag**5)
    c_12 = 3 * mu * r_j * r_k / r_mag**5
    
    c_20 = 3 * mu * r_k * r_i / r_mag**5
    c_21 = 3 * mu * r_k * r_j / r_mag**5
    c_22 = -(mu / r_mag**3) + (3 * mu * r_k**2 / r_mag**5)
    
    C = np.array([
        [c_00, c_01, c_02],
        [c_10, c_11, c_12],
        [c_20, c_21, c_22]
    ])
    
    F0_inter = np.concatenate((A, B), axis=1)
    F1_inter = np.concatenate((C, D), axis=1)
    
    F = np.concatenate((F0_inter, F1_inter), axis=0)
    
    trans_mtx = np.eye(6) + (F * dt) + (F @ F * (dt**2 / 2))
    return trans_mtx

def main() -> None:
    # ----------- Parameters
    tle_file = 'rsc/TLE/navstar_43.txt'
    satellite = None
    seed = 22
    propagation_time = 20
    propagation_step = 0.01
    
    # ----------- Setup
    np.random.seed(seed)
    
    # Load the satellite
    with open(tle_file, 'r') as f:
        name = f.readline().strip()
        line1 = f.readline().strip()
        line2 = f.readline().strip()
        
        satellite = sat.Satellite(
            name,
            line1,
            line2,
            cb.earth
        )
    
    # Propagate the satellite
    r_true, v_true, t_obv = satellite.orbit.propagate(
        propagation_time,
        max_step=propagation_step
    )
    
    # Apply gaussian noise to the position
    r_obv = r_true.copy()
    position_noise = 0.5
    r_obv[0] += np.random.normal(0, position_noise, size=r_obv[0].shape)
    r_obv[1] += np.random.normal(0, position_noise, size=r_obv[1].shape)
    r_obv[2] += np.random.normal(0, position_noise, size=r_obv[2].shape)

    
    # Possibility to add time error
    v_1 = od.herrick_gibbs(
        (r_obv[:,0], t_obv[0]),
        (r_obv[:,1], t_obv[1]),
        (r_obv[:,2], t_obv[2]),
        const.MU_EARTH
    )
    
    initial_state = np.concatenate(
        (r_obv[:,1].flatten(), v_1.flatten())
    ).reshape((6, 1))
    
    process_noise = 0.02 * np.eye(6)
    
    ekf = est.ExtendedKalmanFilter(
        transition_matrix_func,
        initial_state,
        process_noise
    )
    
    # remove observation used for IOD
    t_init = t_obv[1]
    r_obv = r_obv[:,2:]
    t_obv = t_obv[2:]
    
    # ----------- Simulate observations
    t_last = t_init
    r_last = r_obv[:,1]
    observation_covariance = position_noise**2 * np.eye(3)
    obvservation_matrix = np.array([
        [1, 0, 0, 0, 0, 0],
        [0, 1, 0, 0 ,0, 0],
        [0, 0, 1, 0, 0, 0]
    ])
    
    innovations = []
    vel_innovations = []
    obv_innovations = []
    
    x_est, y_est, z_est = [], [], []
    
    for i, t in enumerate(t_obv):
        # Make observation
        dt = t - t_last
        r = r_obv[:,i].reshape((3,1))
        
        state_estimate, _, _ = ekf.update(
            r,
            obvservation_matrix,
            observation_covariance,
            f_args=(ekf.curr_state_est[:3], dt)
        )
        
        t_last = t
        
        innovations.append(
            np.linalg.norm(r_true[:,i+2] - state_estimate[:3].ravel())
        )
        
        vel_innovations.append(
            np.linalg.norm(v_true[:,i+2] - state_estimate[3:].ravel())
        )

    fig, ax = plt.subplots()
    
    ax.plot(t_obv, innovations, label='Innovation')
    
    ax.set_title('Innovation vs Time')
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Innovation (m)')
    
    ax.legend()
    ax.grid()    
    
    plt.show()
    
    return
    

if __name__ == '__main__':
    main()