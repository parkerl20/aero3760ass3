from spacesim import constants as const
from typing import Callable
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

class ExtendedKalmanFilter():
    """A class implementing a general Extended Kalman Filter.
    
    A basic Kalman filter can also be implemented by having the 
    transition matrix function return a constant matrix.
    
    Examples:
    ```python 
        def transition_matrix_func(time: float, r: np.ndarry) -> np.ndarray:
            # A function that returns a time-varying transition matrix
            x, y, z = r     # unpack into position components
            # model matrix
            F = np.array([
                [0, 0, 0, 1, 0, 0],
                [0, 0, 0, 0, 1, 0],
                [0, 0, 0, 0, 0, 1],
                [x, 0, 0, 0, 0, 0],
                [0, y, 0, 0, 0, 0],
                [0, 0, z, 0, 0, 0]
            ])
            
            trans_mtx = np.eye(6) + F * time
            return trans_mtx
        
        ekf = ExtendedKalmanFilter(
            transition_matrix_func,
            np.zeros(6,1),
            np.eye(6),
            np.eye(3)
        )
        
        measured_positions = [[1, 2, 3], [2, 3, 4], [3, 4, 5], ...]
        time_stamps = [0, 1, 2, ...]
        
        # Estimate the state of a system at each time step
        for r, t in zip(measured_positions, time_stamps):
            # observation matrix
            H = np.array([
                [1, 0, 0, 0, 0, 0],
                [0, 1, 0, 0, 0, 0],
                [0, 0, 1, 0, 0, 0]
            ])
            
            state_est, innovation, est_covar = ekf.predict(
                r,
                H,
                f_args=(t, r)
            )

    ```

    """
    def __init__(
        self,
        transition_matrix_func: Callable[[any], np.ndarray],
        initial_state: np.ndarray,
        process_covariance: np.ndarray,
        observation_covariance: np.ndarray
    ) -> None:
        
        self.transition_matrix_func = transition_matrix_func
        self.curr_state_est = initial_state             # state estimate
        self.est_covar = np.eye(initial_state.shape[0]) # state estimate covariance
        self.process_covar = process_covariance         # process covariance
        self.obs_covar = observation_covariance     	# observation covariance
        
        return
    
    def predict(
        self,
        measurement: np.ndarray,
        observation_matrix: np.ndarray,
        *,
        f_args: tuple = ()
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Predicts the next state of the system 

        Args:
            measurement (np.ndarray): The measured state of the system
            observation_matrix (np.ndarray): A matrix mapping from 
                observation space to state space.
            f_args (tuple): Additional arguments to pass to the transition 
                matrix function.

        Returns:
			tuple[np.ndarray, np.ndarray, np.ndarray]: The predicted state,
				the innovation, and the predicted state covariance.
        """
        trans_mtx = self.transition_matrix_func(*f_args)					
        
        state_est = trans_mtx @ self.curr_state_est
        est_covar = trans_mtx @ self.est_covar @ trans_mtx.T + self.process_covar
        
        innovation = measurement - observation_matrix @ state_est		# state residual
        innovation_covar = observation_matrix @ est_covar @ observation_matrix.T + self.obs_covar
        kalman_gain = est_covar @ observation_matrix.T @ np.linalg.inv(innovation_covar)
        
        self.curr_state_est = state_est + kalman_gain @ innovation
        self.est_covar = (np.eye(self.est_covar.shape[0]) - kalman_gain @ observation_matrix) @ est_covar
        
        post_fit_residual = measurement - observation_matrix @ self.curr_state_est
        
        return self.curr_state_est, post_fit_residual, self.est_covar