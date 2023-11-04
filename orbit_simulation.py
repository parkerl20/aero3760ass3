import math
from spacesim import estimation as est
from spacesim import satellite as sat
from spacesim import celestial_body as cb
from spacesim import sensor
from spacesim import constants as const
from spacesim import estimation as est
from spacesim import util
from spacesim import orbital_transforms as ot
from spacesim import orbit as orb

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
    
    # Mapping params
    elevation_deg = 60
    
    
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
    
    star_tracker = sensor.SatelliteSensor(
        "star_tracker",
        np.zeros(3,),       # place holder
        star_tracker_simulator,
        frequency=5
    )
    
    sun_sensor = sensor.SatelliteSensor(
        "sun_sensor",
        np.zeros(3,),       # place holder
        sun_sensor_simulator,
        frequency=5
    )
    
    satellite.attach_sensor(gnss_reciever)
    # satellite.attach_sensor(star_tracker)
    # satellite.attach_sensor(sun_sensor)
    
    # ------------------ Add algorithms
    initial_state = np.concatenate(
        (
            satellite.init_r_eci.flatten(),
            satellite.init_v_eci.flatten()
        )
    ).reshape(-1, 1)
    
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
    
    # Attitude NLLS
    attitude_nlls_algo = sat.SatelliteAlgorithm(
        "Attitude NLLS",
        None,
        attitude_NLLS_algo_function
    )
    
    
    satellite.add_algorithm(ekf_od_algorithm)
    satellite.add_algorithm(attitude_nlls_algo)
    
    # Simulate the satellite
    r_residuals = [[], [], []]
    v_residuals = [[], [], []]
    att_residuals = [[], [], [], []]
    raw_r = [[], [], []]
    time_steps = []    
    gnss_times = []
    mappings = [[], [], [], [], [], []]
    
    
    for r_true, v_true, attitude_true, t in satellite:
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
        
        # Attitude residuals
        att_estimation = None       # NLLS estmation
        
        # r_residuals[0].append(r_true[0] - ekf_r[0])
        # r_residuals[1].append(r_true[1] - ekf_r[1])
        # r_residuals[2].append(r_true[2] - ekf_r[2])
        
        # in_track, cross_track, radial = ot.ECI_to_mapping_error(
        #     ekf_r,
        #     r_true,
        #     v_true
        # )
        
        # TODO: Add when attitude is done
        in_track, cross_track, radial, azimuth, nadir, rms = mapping_budget(
            ekf_r,
            r_true,
            v_true,
            attitude_true,  # needs to be converted to euler (currently in quarternions)
            att_estimation
        )
        
        r_residuals[0].append(in_track)
        r_residuals[1].append(cross_track)
        r_residuals[2].append(radial)
        
        v_residuals[0].append(v_true[0] - ekf_v[0])
        v_residuals[1].append(v_true[1] - ekf_v[1])
        v_residuals[2].append(v_true[2] - ekf_v[2])
        
        time_steps.append(t)
    
        # print("Current attitude: ", satellite.current_attitude)
        att_residuals[0].append(satellite.current_attitude[0])
        att_residuals[1].append(satellite.current_attitude[1])
        att_residuals[2].append(satellite.current_attitude[2])
        att_residuals[3].append(satellite.current_attitude[3])

        mappings[0].append(in_track)
        mappings[1].append(cross_track)
        mappings[2].append(radial)
        mappings[3].append(azimuth)
        mappings[4].append(nadir)
        mappings[5].append(rms)
        
        
    create_od_results(
        r_residuals,
        v_residuals,
        np.array(satellite.algorithms["OD EKF"].logger.get_log("r_residual")).T,
        time_steps
    )

    create_att_results(
        att_residuals,
        time_steps
    )
    
    # Plot mapping errors on seperate axes
    print(f"In track mean: {np.mean(mappings[0])}\tstd: {np.std(mappings[0])}")
    print(f"Cross track mean: {np.mean(mappings[1])}\tstd: {np.std(mappings[1])}")
    print(f"Radial mean: {np.mean(mappings[2])}\tstd: {np.std(mappings[2])}")
    
    in_track_fig, in_track_ax = plt.subplots()
    in_track_ax.plot(time_steps, mappings[0])
    
    in_track_ax.set_title("In-Track Mapping Error")
    in_track_ax.set_xlabel("Time (s)")
    in_track_ax.set_ylabel("Error (m)")
    
    in_track_fig.tight_layout()
    in_track_fig.savefig("figures/od_mapping_in_track.png")
    
    cross_track_fig, cross_track_ax = plt.subplots()
    cross_track_ax.plot(time_steps, mappings[1])
    
    cross_track_ax.set_title("Cross-Track Mapping Error")
    cross_track_ax.set_xlabel("Time (s)")
    cross_track_ax.set_ylabel("Error (m)")
    
    cross_track_fig.tight_layout()
    cross_track_fig.savefig("figures/od_mapping_cross_track.png")
    
    radial_fig, radial_ax = plt.subplots()
    radial_ax.plot(time_steps, mappings[2])
    
    radial_ax.set_title("Radial Mapping Error")
    radial_ax.set_xlabel("Time (s)")
    radial_ax.set_ylabel("Error (m)")
    
    radial_fig.tight_layout()
    radial_fig.savefig("figures/od_mapping_radial.png")

    azimuth_fig, azimuth_ax = plt.subplots()
    azimuth_ax.plot(time_steps, mappings[3])
    
    azimuth_ax.set_title("Azimuth Mapping Error")
    azimuth_ax.set_xlabel("Time (s)")
    azimuth_ax.set_ylabel("Error (m)")
    
    azimuth_fig.tight_layout()
    azimuth_fig.savefig("figures/att_mapping_azimuth.png")

    nadir_fig, nadir_ax = plt.subplots()
    nadir_ax.plot(time_steps, mappings[4])
    
    nadir_ax.set_title("Nadir Mapping Error")
    nadir_ax.set_xlabel("Time (s)")
    nadir_ax.set_ylabel("Error (m)")
    
    nadir_fig.tight_layout()
    nadir_fig.savefig("figures/att_mapping_nadir.png")

    rms_fig, rms_ax = plt.subplots()
    rms_ax.plot(time_steps, mappings[5])
    
    rms_ax.set_title("RMS Mapping Error")
    rms_ax.set_xlabel("Time (s)")
    rms_ax.set_ylabel("Error (m)")
    
    rms_fig.tight_layout()
    rms_fig.savefig("figures/att_mapping_rms.png")
    
    plt.show()
    
    return

def orbit_dynamics(
    t: float,
    state_vector: np.ndarray, 
    orbit: orb.Orbit
) -> np.ndarray:
    mu = orbit.body.gravitational_parameter
    
    r = state_vector[0:3]
    v = state_vector[3:6]
    
    r_mag = np.linalg.norm(r)
    
    q0 = state_vector[6:10].reshape(4,1)     # Quaternion
    
    # Two body equation
    r_dot = v
    v_dot = -(mu / r_mag**3) * r
    
    # Quaternion propagation
    p = 0.0016
    q = 0.0016
    r = 0.0016

    q_matrix = 0.5 * np.array([
        [0, -p, -q, -r],
        [p, 0, r, -q],
        [q, -r, 0, p],
        [r, q, -p, 0]])
    
    q_dot = (q_matrix @ q0).flatten()
    
    return [*r_dot, *v_dot, *q_dot]
    
    

def gnss_reciever_simulator(
    gnss_reciever: sensor.SatelliteSensor,
    satellite: sat.RealTimeSatellite,
    r: np.ndarray,
    v: np.ndarray
) -> bool:
    r = r.flatten()
    v = v.flatten()
    
    # Getting current attitude
    attitude = satellite.current_attitude
    
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

def create_att_results( # -----------------------------------------------------------------------------------------------------
    att_residuals: list[list[float]],
    time_steps: list[float]
) -> None:
    """
    """
    startup_plotting()
    # print("Current Attitude:", att_residuals)
    # print("Size of time:", len(time_steps))
    # print("Size of att_residuals[0]:", len(att_residuals[0]))

    attw_fig, attw_ax = plt.subplots()
    attw_ax.plot(time_steps[1:], att_residuals[0][1:])
    
    attw_ax.set_title("Quaternion[0]")
    attw_ax.set_xlabel("Time (s)")
    attw_ax.set_ylabel("Quaternion Angle (units)")
    
    attw_fig.tight_layout()
    
    attw_fig.savefig("figures/at_quaternion[0].png")


    attx_fig, attx_ax = plt.subplots()
    attx_ax.plot(time_steps[1:], att_residuals[1][1:])
    
    attx_ax.set_title("Quaternion[1]")
    attx_ax.set_xlabel("Time (s)")
    attx_ax.set_ylabel("Quaternion Angle (units)")
    
    attx_fig.tight_layout()
    
    attx_fig.savefig("figures/at_quaternion[1].png")


    atty_fig, atty_ax = plt.subplots()
    atty_ax.plot(time_steps[1:], att_residuals[2][1:])
    
    atty_ax.set_title("Quaternion[2]")
    atty_ax.set_xlabel("Time (s)")
    atty_ax.set_ylabel("Quaternion Angle (units)")
    
    atty_fig.tight_layout()
    
    atty_fig.savefig("figures/at_quaternion[2].png")


    attz_fig, attz_ax = plt.subplots()
    attz_ax.plot(time_steps[1:], att_residuals[3][1:])
    
    attz_ax.set_title("Quaternion[3]")
    attz_ax.set_xlabel("Time (s)")
    attz_ax.set_ylabel("Quaternion Angle (units)")
    
    attz_fig.tight_layout()
    
    attz_fig.savefig("figures/at_quaternion[3].png")

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
    
    # plt.show() 
    return

def sun_sensor_simulator(
    sun_sensor: sensor.SatelliteSensor,
    satellite: sat.RealTimeSatellite,
    r: np.ndarray,
    v: np.ndarray
) -> bool:
    """Simulates the sun sensor on a satellite"""
    # Get ECI position estimation
    attitude_quart = satellite.current_attitude
    
    # 0.1 deg error std
    
    # Convert to euler angles
    euler = ot.quat2euler(attitude_quart.flatten())
    euler += np.random.normal(0, 0.1, 3)
    
    sun_sensor.mesurement = euler
    return True

def star_tracker_simulator(
    star_tracker: sensor.SatelliteSensor,
    satellite: sat.RealTimeSatellite,
    r: np.ndarray,
    v: np.ndarray
) -> bool:
    """Simulates the star tracker on a satellite"""
    attitude_quart = satellite.current_attitude
    
    yaw_noise = 0.03
    pitch_noise = 0.002
    roll_noise = 0.002
    
    # Convert to euler angles
    euler = ot.quat2euler(attitude_quart.flatten())
    
    # Apply noise
    euler[0] += np.random.normal(0, yaw_noise)
    euler[1] += np.random.normal(0, pitch_noise)
    euler[2] += np.random.normal(0, roll_noise)

    star_tracker.mesurement = euler
    
    return True

def calculate_rms(
    data: np.array
) -> float:
    """
    Calculate the root mean square of a list of values.
    """
    n = len(data)
    if n < 1:
        raise ValueError("Input list must contain at least one element.")

    sum_of_squares = sum(x ** 2 for x in data)
    rms = math.sqrt(sum_of_squares / n)
    return rms

def mapping_budget(
    r_eci, 
    r_true, 
    v_true, 
    att_true, 
    att_estimate
):
    
    # print("att_true:", att_true)

    delta_I, delta_C, delta_R = ot.ECI_to_mapping_error(r_eci, r_true, v_true)
    delta_azimuth, delta_elevation = ot.ECI_to_azimuth_error(att_true, att_estimate)


    R_E = const.R_EARTH     # m
    H = np.linalg.norm(r_eci) - R_E
    R_T = R_E
    R_S = R_E + H

    elevation_deg = 60
    elevation_rad = np.deg2rad(elevation_deg)

    delta_phi_deg = 0.0005
    delta_phi_rad = np.deg2rad(delta_phi_deg)

    delta_eta_deg = 0.0005
    delta_eta_rad = np.deg2rad(delta_eta_deg)

    p_rad = np.arcsin(R_E/(R_E+H))
    p_deg = np.rad2deg(p_rad)

    eta_rad = np.arcsin(np.cos(elevation_rad) * np.sin(p_rad))
    eta_deg = np.rad2deg(eta_rad)

    lam_deg = 90 - elevation_deg - p_deg
    lam_rad = np.deg2rad(lam_deg)

    phi_deg = 30
    phi_rad = np.deg2rad(phi_deg)

    H_mapping = np.arcsin(np.sin(lam_rad)*np.sin(phi_rad))
    G_mapping = np.arcsin(np.sin(lam_rad)*np.cos(phi_rad))
 
    intrack_error = delta_I * R_T / R_S * np.cos(H_mapping)
    crosstrack_error = delta_C * R_T / R_S * np.cos(G_mapping)
    radial_error = delta_R * np.sin(eta_rad) / np.sin(elevation_rad)

    division = 0.0451649 # for elevation 60 deg (worst case)
    D = R_E * division

    azimuth_error = delta_azimuth * D * np.sin(eta_rad)
    nadir_error = delta_elevation * D / np.sin(elevation_rad)

    data = [azimuth_error, nadir_error, intrack_error, crosstrack_error, radial_error]

    rms = calculate_rms(data)

    return intrack_error, crosstrack_error, radial_error, azimuth_error, nadir_error, rms

def attitude_NLLS_algo_function(
    time: float,
    nlls: sat.SatelliteAlgorithm,
    satellite: sat.RealTimeSatellite,
    sensors: dict[str, sensor.SatelliteSensor],
    logger: util.DictLogger = None
) -> None:
    # nlls.algorithm == None
    # returns None if star tracker not in sensors dictionary
    star_tracker = sensors.get("star_tracker", None)
    sun_sensor = sensors.get("sun_sensor", None)
    
    # One way
    # star_tracker_attitude = None
    # sun_sensor_attitude = None
    
    if star_tracker is not None:
        # do start tracker stuff
        star_tracker_attitude = star_tracker.get_measurement()
        pass
    
    if sun_sensor is not None:
        # do sun sensor stuff
        sun_sensor_attitude = sun_sensor.get_measurement()
        pass
    
    pass