from spacesim import satellite as sat
from spacesim import estimation as est
from spacesim import constants as const
from spacesim import orbit_determination as od
import numpy as np


def transition_matrix_func(r: np.ndarray, t: float) -> np.ndarray:
    pass

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