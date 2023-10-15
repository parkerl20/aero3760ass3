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
from spacesim import orbit as orb

from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import MultipleLocator
from pyvista import examples as pv_ex
import pyvista as pv
import matplotlib.pyplot as plt
import numpy as np
import datetime


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

def perturbation_dynamics(t: float, state_vector: np.ndarray, orbit: orb.Orbit) -> np.ndarray:
    # Satellite dynamics that accounts for perturbations
    def perturb_acceleration(r: np.ndarray) -> np.ndarray:
        # Inner function to calculate perturbation accleration
        x, y, z = r
        r_mag = np.sqrt(r.dot(r))
        
        P = - (3 * mu * J2 * R**2) / (2 * r_mag**5)
        a_x = (1 - 3 * z**2 / r_mag**2) * x
        a_y = (1 - 3 * z**2 / r_mag**2) * y
        a_z = (3 - 3 * z**2 / r_mag**2) * z
        
        return P * np.array([a_x, a_y, a_z], dtype=np.float64)
    
    mu = orbit.body.gravitational_parameter
    J2 = orbit.body.J2
    R = orbit.body.radius
    
    r = state_vector[0:3]
    v = state_vector[3:6]
    
    r_mag = np.sqrt(r.dot(r))
    
    # Dynamics
    P = perturb_acceleration(r)
    
    r_dot = v
    v_dot = -(mu / (r_mag ** 3)) * r + P    # two-body equation
    
    return [*r_dot, *v_dot]

def main() -> None:
    # ----------- Parameters
    seed = 38
    propagation_time = 20
    propagation_step = 0.01
    
    # Orbit params
    a = 6932386.57371916
    e = 0
    i = -33
    arg_p = 0
    RAAN = 238.82
    theta = 0
    epoch = datetime.datetime(2023, 1, 1)
    
    # ----------- Setup
    np.random.seed(seed)
    earth_mesh = pv_ex.planets.load_earth(radius=const.R_EARTH)
    earth_tex = pv_ex.load_globe_texture()
    
    np.set_printoptions(linewidth=np.inf)
    
    earth = cb.CelestialBody(
        "Earth",
        const.M_EARTH,
        const.R_EARTH,
        J2 = const.J2_EARTH,
        colour="blue",
        body_mesh=earth_mesh,
        body_texture=earth_tex,
    )
    
    orbit = orb.Orbit(
        a,
        e,
        i,
        RAAN,
        arg_p,
        theta,
        earth,
        epoch,
        "SPY",
        # orb_dyn=perturbation_dynamics,
        colour="red"
    )
    
    # Propagate the satellite
    r_true, v_true, t_obv = orbit.propagate(
        propagation_time,
        max_step=propagation_step
    )
    
    # Apply gaussian noise to the position
    r_obv = r_true.copy()
    position_noise = 0.87
    r_obv[0] += np.random.normal(0, position_noise, size=r_obv[0].shape)
    r_obv[1] += np.random.normal(0, position_noise, size=r_obv[1].shape)
    r_obv[2] += np.random.normal(0, position_noise, size=r_obv[2].shape)

    
    # Possibility to add time error
    obv_index = 40
    
    v_1 = od.herrick_gibbs(
        (r_obv[:,0], t_obv[0]),
        (r_obv[:,obv_index], t_obv[obv_index]),
        (r_obv[:,80], t_obv[80]),
        const.MU_EARTH
    )
    
    # print(f"predict diff: {v_1.flatten() - v_true[:,obv_index].flatten()}")
    
    initial_state = np.concatenate(
        (r_obv[:,obv_index].flatten(), v_1.flatten())
    ).reshape((6, 1))
    
    process_noise = 5 * np.eye(6)
    
    ekf = est.ExtendedKalmanFilter(
        transition_matrix_func,
        initial_state,
        process_noise
    )
    
    # remove observation used for IOD
    t_init = t_obv[obv_index]
    r_obv = r_obv[:,obv_index+1:]
    t_obv = t_obv[obv_index+1:]
    
    # ----------- Simulate observations
    t_last = t_init
    # r_last = r_obv[:,1]
    observation_covariance = position_noise**2 * np.eye(3)
    obvservation_matrix = np.array([
        [1, 0, 0, 0, 0, 0],
        [0, 1, 0, 0 ,0, 0],
        [0, 0, 1, 0, 0, 0]
    ])
    
    innovations = []
    vel_innovations = []
    predicted = []
    innovation_sum = np.zeros(3,)
    innovation_averages = []
    x_innovation = []
    y_innovation = []
    z_innovation = []
    
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
            np.linalg.norm(r_true[:,i+  obv_index + 1] - state_estimate[:3].ravel())
        )
        
        vel_innovations.append(
            np.linalg.norm(v_true[:,i + obv_index + 1] - state_estimate[3:].ravel())
        )
        
        innovation_sum += r_true[:,i + obv_index + 1] - state_estimate[:3].ravel()
        innovation_averages.append(
            np.linalg.norm(innovation_sum / (i + 1))
        )
        
        x_innovation.append(
            r_true[0,i + obv_index + 1] - state_estimate[0]
        )
        
        y_innovation.append(
            r_true[1,i + obv_index + 1] - state_estimate[1]
        )
        
        z_innovation.append(
            r_true[2,i + obv_index + 1] - state_estimate[2]
        )
        
        predicted.append(state_estimate[:3].ravel())
    
    # predicted = np.array(predicted).T
    x_fig, x_ax = plt.subplots()
    y_fig, y_ax = plt.subplots()
    z_fig, z_ax = plt.subplots()
    
    x_ax.plot(t_obv[500:1000], x_innovation[500:1000], label='Innovation')
    y_ax.plot(t_obv[500:1000], y_innovation[500:1000], label='Innovation')
    z_ax.plot(t_obv[500:1000], z_innovation[500:1000], label='Innovation')
    
    x_ax.set_title('X - Innovation vs Time')
    y_ax.set_title('Y - Innovation vs Time')
    z_ax.set_title('Z - Innovation vs Time')
    
    x_ax.set_xlabel('Time (s)')
    y_ax.set_xlabel('Time (s)')
    z_ax.set_xlabel('Time (s)')
    
    x_ax.set_ylabel('Innovation (m)')
    y_ax.set_ylabel('Innovation (m)')
    z_ax.set_ylabel('Innovation (m)')
    
    x_ax.grid()
    y_ax.grid()
    z_ax.grid()
    
    x_fig.tight_layout()
    y_fig.tight_layout()
    z_fig.tight_layout()
    
    x_fig.savefig("position/plots/x_innovation.png")
    y_fig.savefig("position/plots/y_innovation.png")
    z_fig.savefig("position/plots/z_innovation.png")
    
    vel_fig, vel_ax = plt.subplots()
    vel_ax.plot(t_obv, vel_innovations, label='Innovation')
    vel_ax.set_title('Velocity Innovation vs Time')
    vel_ax.set_xlabel('Time (s)')
    vel_ax.set_ylabel('Innovation (m/s)')
    
    vel_ax.grid()
    vel_fig.tight_layout()
    vel_fig.savefig("position/plots/vel_innovation.png")

    fig, ax = plt.subplots()
    
    ax.plot(t_obv, innovations, label='Innovation')
    ax.plot(t_obv, innovation_averages, label='Innovation Average')
    
    ax.set_title('Innovation vs Time')
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Innovation (m)')
    
    ax.legend()
    ax.grid()
    fig.tight_layout()
    fig.savefig("position/plots/innovation.png")
    
    # Plot predicted and true orbits in 3D
    fig_3d = plt.figure()
    ax_3d = fig_3d.add_subplot(111, projection='3d')
    
    start_i = 0
    end_i = 5
    
    ax_3d.plot(
        r_true[0][start_i+2:end_i+2],
        r_true[1][start_i+2:end_i+2],
        r_true[2][start_i+2:end_i+2], 
        label='True Orbit'
    )
    
    ax_3d.plot(
        predicted[0][start_i:end_i],
        predicted[1][start_i:end_i],
        predicted[2][start_i:end_i],
        label='Predicted Orbit'
    )
    
    ax_3d.legend()
    fig_3d.tight_layout()
    
    plt.show()
    
    return
    

if __name__ == '__main__':
    main()