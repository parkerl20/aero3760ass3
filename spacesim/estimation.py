from spacesim import constants as const
import numpy as np


def linear_least_squares(
    y: np.ndarray,
    H: np.ndarray,
    W: np.ndarray = None
) -> np.ndarray:
    """Performs linear least squares estimation.

    Args:
        y (np.ndarray): An (n x 3) vector of observations.
        H (np.ndarray): The system model
        W (np.ndarray, optional): The weighting matrix. If not provided, the
            identity matrix is used.

    Returns:
        np.ndarray: The estimated state vector.
    """
    if W is None:
        W = np.eye(H.shape[0])
        
    return np.linalg.inv(H.T @ W @ H) @ H.T @ W @ y


class OrbitKalmanFilter():
    """An Kalman filter for a Keplerian orbit.
    """

    def __init__(self, init_r, init_v, init_t, mu: float = const.MU_EARTH) -> None:
        self.init_r = init_r.flatten()
        self.init_v = init_v.flatten()
        self.t = init_t
        
        self.current_r = self.init_r
        self.current_v = self.init_v
        self.mu = mu
        
        self.F = OrbitKalmanFilter.state_matrix(self.init_r, self.mu)
        self.current_x = np.concatenate((self.init_r, self.init_v)).reshape((6, 1))
    

        #------------- Kalman filter matrices
        # Measurement matrix
        self.H = np.array([
            [1, 0, 0, 0, 0, 0],
            [0, 1, 0, 0 ,0 ,0],
            [0, 0, 1, 0 ,0 ,0]
        ])
        
        # Variance matrices
        self.P = np.eye(6)
        self.Q = 0.1 * np.eye(6)
        
        self.R = 30 * np.eye(3)
        
        # Noise matrices - not modeling noise yet
        self.w = 0 * np.eye(6)
        self.v = 10 * np.eye(6)
        
        return
    
    def predict(self, r: np.ndarray, t: float) -> np.ndarray:
        """Predicts the state vector at the next time step.

        Args:
            r (np.ndarray): The observation vector.
            t (float): The time of the observation.

        Returns:
            np.ndarray: A (6 x 1) vector of the predicted state.
        """
        dt = t - self.t
        self.t = t
        
        phi = np.eye(6) + self.F * dt + self.F @ self.F * (dt**2 / 2)
        X_bar = phi @ self.current_x                # predicted state
        P_bar = phi @ self.P @ phi.T + self.Q       # predicted covariance
        
        y_tilde = r - self.H @ X_bar
        S = self.H @ P_bar @ self.H.T + self.R      # innovation covariance
        K = P_bar @ self.H.T @ np.linalg.inv(S)     # Kalman gain
        delta_x = K @ y_tilde                       # state correction
        
        # Update states
        self.current_x = X_bar + delta_x
        self.P = (np.eye(6) - K @ self.H) @ P_bar
        
        self.F = OrbitKalmanFilter.state_matrix(self.current_x[:3, 0], self.mu)
        
        return self.current_x.flatten()
        

    @staticmethod
    def state_matrix(r: np.ndarray, mu: float = const.MU_EARTH) -> np.ndarray:
        """_summary_

        Args:
            r (np.ndarray): A position vector
            mu (float, optional): Gravitational parameter of the primary body.
                Defaults to const.MU_EARTH.

        Returns:
            np.ndarray: A (6 x 6) matrix
        """
        r_i, r_j, r_k = r
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
        
        return np.concatenate((F0_inter, F1_inter), axis=0)
        
    
    @staticmethod
    def __universal_conic_section_orbit(
        t: float, 
        init_pos: np.ndarray, 
        init_vel: np.ndarray,
        mu: float = const.MU_EARTH
    ) -> tuple[np.ndarray, np.ndarray]:
        """
        UNIVERSAL_CONIC_SECTION_ORBIT - determines the position and velocity of a
        satellite given the orbital parameters and an input time using a
        universal conic section solution to Kepler's equations (from initial
        position and velocity)

        Inputs:
            t (float): time since epoch in seconds
            init_pos (numpy.ndarray): Satellite position at epoch (3x1 numpy array)
            init_vel (numpy.ndarray): Satellite velocity at epoch (3x1 numpy array)

        Outputs:
            Position_ECI (numpy.ndarray): ECI frame position at time t (3, numpy array)
            Velocity_ECI (numpy.ndarray): ECI frame velocity at time t (3, numpy array)

        See Week 6 slides on Orbital Determination using NLLS for details.
        """

        # Analytically solve change in eccentric anomaly equation for position and
        # velocity at time t
        a_semi = 1 / (2 / np.linalg.norm(init_pos) - np.linalg.norm(init_vel) ** 2 / mu)
        phi = OrbitKalmanFilter.__deltaE_solve(a_semi, mu, init_pos, init_vel, t)

        # Compute intermediate variables
        f = 1 - a_semi * (1 - np.cos(phi)) / np.linalg.norm(init_pos)
        g = t - (a_semi ** (3 / 2)) * (phi - np.sin(phi)) / np.sqrt(mu)
        r = a_semi * (1 - (1 - (np.linalg.norm(init_pos) / a_semi)) * np.cos(phi)) + \
            np.dot(init_pos.T, init_vel) * np.sqrt(a_semi / mu) * np.sin(phi)
        fdot = -np.sqrt(mu * a_semi) * np.sin(phi) / (r * np.linalg.norm(init_pos))
        gdot = 1 - a_semi * (1 - np.cos(phi)) / r

        # Compute position and velocity at time t
        position_eci = f * init_pos + g * init_vel
        velocity_eci = fdot * init_pos + gdot * init_vel

        return position_eci, velocity_eci

    @staticmethod
    def __deltaE_solve(a, mu, r, rdot, delta_t):
        # DELTAE_SOLVE - solve equation for change in eccentric anomaly using a
        # Newton/Raphson method

        tol = 1e-12
        phik = 0
        del_phi = 1
        i = 0
        while del_phi > tol:
            phi_old = phik
            B = (a ** (3 / 2)) / np.sqrt(mu)
            C = 1 - np.linalg.norm(r) / a
            D = np.dot(r, rdot) / np.sqrt(mu * a)
            fn = B * (phik - C * np.sin(phik) + D * (1 - np.cos(phik))) - delta_t
            fndot = B - B * C * np.cos(phik) + B * D * np.sin(phik)
            phik = phik - (fn / fndot)
            del_phi = np.abs(phi_old - phik) % (2 * np.pi)
            i += 1

        return phik