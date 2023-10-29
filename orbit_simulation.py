from spacesim import estimation as est
from spacesim import satellite as sat
from spacesim import celestial_body as cb
from spacesim import sensor
from spacesim import constants as const
from spacesim import estimation as est
from spacesim import util
from spacesim import orbital_transforms as ot

import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
from util import startup_plotting

def orbit_simulation(
    propagation_time: float,
    a: float,
    e: float,
    i: float,
    RAAN: float,
    arg_p: float,
    true_anom: float,
    epoch: dt.datetime
) -> None:
    """Performs orbit determination simulation for a satellite

    Args:
        satellite (sat.RealTimeSatellite): The satellite to perform
            orbit determination on
    """
    # ----- set up
    np.random.seed(32)
    propagation_length = 20
    
    satellite = sat.RealTimeSatellite(
        a,
        e,
        i,
        RAAN,
        arg_p,
        true_anom,
        cb.earth,
        epoch,
        name="Satellite 1",
        propagation_length=propagation_length,
        propagation_step=propagation_length
    )
    
    propagation_time = satellite.period
    
    
    # ---------------- Mount sensors to satellite
    # Simulates Rakon GNSS reciever DB
    gnss_position_noise = 0.7 * np.ones(3,)
    gnss_velocity_noise = 0.03 * np.ones(3,)
    gnss_reciever = sensor.SatelliteSensor(
        "gnss_reciever",
        np.concatenate((gnss_position_noise, gnss_velocity_noise)),
        gnss_reciever_simulator,
        frequency=20
    )
    
    satellite.attach_sensor(gnss_reciever)
    
    # ------------------ Add algorithms
    initial_state = np.concatenate(
        (
            satellite.init_r_eci.flatten(),
            satellite.init_v_eci.flatten()
        )
    ).reshape(6, 1)
    
    process_noise = 0.15 * np.eye(6)
    
    ekf = est.ExtendedKalmanFilter(
        EKF_transition_matrix_func,
        initial_state,
        process_noise
    )
    
    ekf_od_algorithm = sat.SatelliteAlgorithm(
        "OD EKF",
        ekf,
        EKF_algo_function,
        util.DictLogger()
    )
    
    satellite.add_algorithm(ekf_od_algorithm)
    
    # Simulate the satellite
    r_residuals = [[], [], []]
    v_residuals = [[], [], []]
    raw_r = [[], [], []]
    time_steps = []    
    gnss_times = []
    
    
    for r_true, v_true, t in satellite:
        if t > propagation_time:
            break
        
        od_ekf = satellite.algorithms["OD EKF"]
        
        dt = t - od_ekf.t_last
        ekf_state = od_ekf.algorithm.predict_state(
            f_args=(
                od_ekf.algorithm.get_state()[:3],
                dt
            )
        )
        
        # Orbit determination residuals
        ekf_r = ekf_state[:3].flatten()
        ekf_v = ekf_state[3:].flatten()
        
        x, y, z = ekf_r
        
        # r_residuals[0].append(r_true[0] - ekf_r[0])
        # r_residuals[1].append(r_true[1] - ekf_r[1])
        # r_residuals[2].append(r_true[2] - ekf_r[2])
        
        in_track, cross_track, radial = ot.ECI_to_mapping_error(
            ekf_r,
            r_true,
            v_true
        )
        
        r_residuals[0].append(in_track)
        r_residuals[1].append(cross_track)
        r_residuals[2].append(radial)
        
        v_residuals[0].append(v_true[0] - ekf_v[0])
        v_residuals[1].append(v_true[1] - ekf_v[1])
        v_residuals[2].append(v_true[2] - ekf_v[2])
        
        time_steps.append(t)
     
    create_od_results(
        r_residuals,
        v_residuals,
        np.array(satellite.algorithms["OD EKF"].logger.get_log("r_residual")).T,
        time_steps
    )

    return

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
    sensors: dict[str, sensor.SatelliteSensor],
    logger: util.DictLogger = None
) -> None:
    """Simulates an Extended Kalman Filter algorithm
    on a satellite in that acts on sensor data.
    
    The Extneded Kalman Filter algorithm function for
    orbit determination recieves data from the GNSS
    reciever sensor.
    """
    gnss_reciever = sensors.get("gnss_reciever", None)
    
    dt = time - ekf.t_last
      
    if gnss_reciever is not None:
        measured_state = gnss_reciever.get_measurement().reshape(6, 1)
        gnss_H = np.eye(6)
        
        gnss_noise = np.diag(gnss_reciever.noise_std**2) * 1
        
        ekf.algorithm.update(
            measured_state,
            gnss_H,
            gnss_noise,
            f_args=(ekf.algorithm.get_state()[:3], dt)
        )
    
    # Log results
    if logger is not None:
        log_EKF_algo(
            logger,
            time,
            ekf.algorithm.get_state()[:3].flatten(),
            ekf.algorithm.get_state()[3:].flatten(),
            satellite.current_r_eci.flatten(),
            satellite.current_v_eci.flatten(),
        )
    
    return

def log_EKF_algo(
    logger: util.DictLogger,
    t: float,
    r: np.ndarray,
    v: np.ndarray,
    r_true: np.ndarray,
    v_true: np.ndarray,
) -> None:
    """Logging function for the OD EKF algorithm
    """
    if len(logger.get_log_keys()) == 0:
        logger.add_log("t", [])
        logger.add_log("r_residual", [])
        logger.add_log("v_residual", [])
    
    logger.get_log("t").append(t)
    logger.get_log("r_residual").append(r - r_true)
    logger.get_log("v_residual").append(v - v_true)
    
    return

def create_od_results(
    r_residuals: list[list[float]],
    v_residuals: list[list[float]],
    r_log: np.ndarray,
    time_steps: list[float]
) -> None:
    """Creates plots for the residuals in position and velocity
    """
    startup_plotting()
    print(f"x - mean: {np.mean(r_log[0])}\tstd: {np.std(r_log[0])}")
    print(f"y - mean: {np.mean(r_log[1])}\tstd: {np.std(r_log[1])}")
    print(f"z - mean: {np.mean(r_log[2])}\tstd: {np.std(r_log[2])}")
    
    # Plot the residuals
    
    # Plot all on one graph
    # r_fig, r_ax = plt.subplots()
    
    # r_ax.plot(time_steps, r_residuals[0], label="x")
    # r_ax.plot(time_steps, r_residuals[1], label="y")
    # r_ax.plot(time_steps, r_residuals[2], label="z")
    
    # r_ax.set_title("Residuals in position")
    # r_ax.set_xlabel("Time (s)")
    # r_ax.set_ylabel("Residual (m)")
    
    # r_ax.legend()
    # r_fig.tight_layout()
    
    # r_fig.savefig("4-Plots/od_r_residuals.png")
    
    rx_fig, rx_ax = plt.subplots()
    rx_ax.plot(time_steps, r_residuals[0])
    
    
    rx_ax.set_title("EKF In-Track Error")
    rx_ax.set_xlabel("Time (s)")
    rx_ax.set_ylabel("Residual (m)")
    
    rx_fig.tight_layout()
    
    rx_fig.savefig("figures/od_EKF_in_track.png")
    
    ry_fig, ry_ax = plt.subplots()
    
    ry_ax.plot(time_steps, r_residuals[1])
    
    ry_ax.set_title("EKF Cross-Track Error")
    ry_ax.set_xlabel("Time (s)")
    ry_ax.set_ylabel("Residual (m)")
    
    ry_fig.tight_layout()
    ry_fig.savefig("figures/od_EKF_cross_track.png")
    
    rz_fig, rz_ax = plt.subplots()
    
    rz_ax.plot(time_steps, r_residuals[2])
    
    rz_ax.set_title("EKF Radial Error")
    rz_ax.set_xlabel("Time (s)")
    rz_ax.set_ylabel("Residual (m)")
    
    rz_fig.tight_layout()
    rz_fig.savefig("figures/od_EKF_radial_error.png")
    
    plt.show() 
    return

def sun_sensor_simulator(
    sun_sensor: sensor.SatelliteSensor,
    satellite: sat.RealTimeSatellite,
    r: np.ndarray,
    v: np.ndarray
) -> bool:
    """Simulates the sun sensor on a satellite"""
    # Get ECI position estimation
    position_ekf = satellite.algorithms["OD EKF"]
    
    ekf_dt = satellite.current_time - position_ekf.t_last
    eci_est = position_ekf.algorithm.predict_state(
        f_args=(
            position_ekf.algorithm.get_state()[:3],
            ekf_dt
        )
    )
    
    current_epoch = satellite.epoch + dt.timedelta(seconds=satellite.current_time)

def star_tracker_simulator(
    star_tracker: sensor.SatelliteSensor,
    satellite: sat.RealTimeSatellite,
    r: np.ndarray,
    v: np.ndarray
) -> bool:
    """Simulates the star tracker on a satellite"""
    pass
