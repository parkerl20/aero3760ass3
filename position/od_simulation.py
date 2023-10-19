import sys
import os

current_file_path = os.path.abspath(__file__)
parent_directory = os.path.dirname(current_file_path)
grandparent_directory = os.path.dirname(parent_directory)
sys.path.append(grandparent_directory)

from spacesim import estimation as est
from spacesim import satellite as sat
from spacesim import celestial_body as cb
from spacesim import sensor
from spacesim import constants as const
from spacesim import estimation as est

import numpy as np

def od_simulation(
    satellite: sat.RealTimeSatellite,
    r_init: np.ndarray,
    v_init: np.ndarray,
) -> None:
    """Performs orbit determination simulation for a satellite

    Args:
        satellite (sat.RealTimeSatellite): The satellite to perform
            orbit determination on
    """
    # ---------------- Mount sensors to satellite
    # Simulate HG1700 IMU
    imu_noise = 0.065 * np.ones(3,)
    imu = sensor.SatelliteSensor(
        "imu",
        imu_noise,
        imu_simulator,
        frequency=100  
    )
    
    # Simulates Rakon GNSS reciever DB
    gnss_position_noise = 0.7 * np.ones(3,)
    gnss_velocity_noise = 0.03 * np.ones(3,)
    gnss_reciever = sensor.SatelliteSensor(
        "gnss_reciever",
        np.concatenate((gnss_position_noise, gnss_velocity_noise)),
        gnss_reciever_simulator,
        frequency=10
    )
    
    satellite.attach_sensor(imu)
    satellite.attach_sensor(gnss_reciever)
    
    # ------------------ Add Extended Kalman Filter algorithm
    initial_state = np.concatenate(
        (r_init.flatten(), v_init.flatten())
    )
    
    process_noise = 1 * np.eye(6)
    
    ekf = est.ExtendedKalmanFilter(
        EKF_transition_matrix_func,
        initial_state,
        process_noise
    )
    
    ekf_od_algorithm = sat.SatelliteAlgorithm(
        "Orbit Determination EKF",
        ekf,
        EKF_algo_function
    )
    
    satellite.add_algorithm(ekf_od_algorithm)

    return


def imu_simulator(
    imu: sensor.SatelliteSensor,
    satellite: sat.RealTimeSatellite,
    r: np.ndarray,
    v: np.ndarray
) -> bool:
    # Use the two body equation
    mu = satellite.body.gravitational_parameter
    r_mag = np.linalg.norm(r)
    
    a = (-(mu / r_mag**3) * r).flatten()

    # Add noise
    a[0] += np.random.normal(
        0,
        imu.noise_std[0]
    )
    
    a[1] += np.random.normal(
        0,
        imu.noise_std[1]
    )
    
    a[2] += np.random.normal(
        0,
        imu.noise_std[2]
    )
    
    imu.mesurement = a
    return True
    
def gnss_reciever_simulator(
    gnss_reciever: sensor.SatelliteSensor,
    satellite: sat.RealTimeSatellite,
    r: np.ndarray,
    v: np.ndarray
) -> bool:
    r = r.flatten()
    v = v.flatten()
    
    r[0] += np.random.normal(
        0,
        gnss_reciever.noise_std[0]
    )
    
    r[1] += np.random.normal(
        0,
        gnss_reciever.noise_std[1]
    )
    
    r[2] += np.random.normal(
        0,
        gnss_reciever.noise_std[2]
    )
    
    v[0] += np.random.normal(
        0,
        gnss_reciever.noise_std[3]
    )
    
    v[1] += np.random.normal(
        0,
        gnss_reciever.noise_std[4]
    )
    
    v[2] += np.random.normal(
        0,
        gnss_reciever.noise_std[5]
    )
    
    gnss_reciever.mesurement = np.concatenate((r, v))
    return True

def EKF_transition_matrix_func(r: np.ndarray, dt: float) -> np.ndarray:    
    mu = const.MU_EARTH
    r_i, r_j, r_k = r.ravel()
    r_mag = np.linalg.norm(r)
    
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

def EKF_algo_function(
    time: float,
    ekf: sat.SatelliteAlgorithm,
    satellite: sat.RealTimeSatellite,
    sensors: dict[str, sensor.SatelliteSensor]
) -> None:
    """Simulates an Extended Kalman Filter algorithm
    on a satellite in that acts on sensor data.
    
    The Extneded Kalman Filter algorithm function for
    orbit determination recieves data from the IMU and GNSS
    reciever sensors.
    """
    imu = sensors.get("imu", None)
    gnss_reciever = sensors.get("gnss_reciever", None)
    
    # Velocity estimate from previous IMU measurement
    if imu is not None:
        if hasattr(EKF_algo_function, "last_imu"):
            last_v = ekf.algorithm.get_state()[3:].flatten()
            last_imu = EKF_algo_function.last_imu
            dt = time - ekf.t_last
            
            predicted_v = last_v + (last_imu * dt)
            
            imu_H = np.array([
                [0, 0, 0, 1, 0, 0],
                [0, 0, 0, 0, 1, 0],
                [0, 0, 0, 0, 0, 1]
            ])
            
            ekf.algorithm.update(
                predicted_v.reshape(3, 1),
                imu_H,
                np.diag(imu.noise_std**2),
                f_args=(ekf.algorithm.get_state()[:3], dt)
            )
        
        EKF_algo_function.last_imu = imu.get_measurement()
    
    if gnss_reciever is not None:
        measured_state = gnss_reciever.get_measurement().reshape(6, 1)
        gnss_H = np.eye(6)
        
        ekf.algorithm.update(
            measured_state,
            gnss_H,
            np.diag(gnss_reciever.noise_std**2)
        )
    
    return