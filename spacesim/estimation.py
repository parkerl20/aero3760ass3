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