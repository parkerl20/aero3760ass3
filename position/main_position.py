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
import numpy as np


def transition_matrix_func(r: np.ndarray, dt: float) -> np.ndarray:
    mu = const.MU_EARTH
    r_i, r_j, r_k = r.ravel()
    r_mag = np.linalg.norm(r)
    
    # intermediate matrices
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
    tle_file = 'tle.txt'
    satellite = sat.Satellite(tle_file)
    
    propagation_time = const.solar_day_length
    
    # Propagate the satellite
    r_true,_,t = satellite.orbit.propagate(propagation_time)
    
    # Apply gaussian noise to the position
    r_obv = r_true.copy()
    position_noise = 1
    r_obv[0] += np.random.normal(0, position_noise)
    r_obv[1] += np.random.normal(0, position_noise)
    r_obv[2] += np.random.normal(0, position_noise)
    
    # Possibility to add time error
    v_1 = od.gibbs(r[:,0], r[:,1], r[:,2])
    initial_state = np.concatenate((r[:,0], v_1))
    
    # Simulate observations
    ekf = est.ExtendedKalmanFilter(
        transition_matrix_func,
        initial_state,
        np.eye(6)
    )
    
    for r in r_obv:
        estimated_state = ekf.predict(
            r,
            np.eye(3),
            f_args=(r, t)
        )
    
    return
    

if __name__ == '__main__':
    print(const.J2000)