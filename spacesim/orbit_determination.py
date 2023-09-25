import numpy as np
import math


def gibbs(
    r1: np.ndarray,
    r2: np.ndarray,
    r3:np.ndarray,
    mu: float,
    tolerance: float = 1e-2
) -> np.ndarray:
    """The gibbs method of orbit determination.

    Args:
        r1 (np.ndarray): A (3 x 1) position vector.
        r2 (np.ndarray): A (3 x 1) position vector.
        r3 (np.ndarray): A (3 x 1) position vector.
        mu (float): The gravitational parameter of the primary body.

    Returns:
        np.ndarray: The (3 x 1) velocity vector at r2.
    """
    r1 = r1.ravel()
    r2 = r2.ravel()
    r3 = r3.ravel()
    
    r1_mag = np.linalg.norm(r1)
    r2_mag = np.linalg.norm(r2)
    r3_mag = np.linalg.norm(r3)
    
    #----- Co-planar check
    C_23 = np.cross(r2, r3) / np.linalg.norm(np.cross(r2, r3))
    u_1 = r1 / r1_mag
    
    if np.abs(np.dot(C_23, u_1)) > tolerance:
        return np.zeros((3, 1))
    
    N = r1_mag * np.cross(r2, r3) + r2_mag * np.cross(r3, r1) + r3_mag * np.cross(r1, r2)
    N_mag = np.linalg.norm(N)
    
    D = np.cross(r1, r2) + np.cross(r2, r3) + np.cross(r3, r1)
    D_mag = np.linalg.norm(D)
    
    S = (r2_mag - r3_mag) * r1 + (r3_mag - r1_mag) * r2 + (r1_mag - r2_mag) * r3
    
    v2 = math.sqrt(mu / (N_mag * D_mag)) * (np.cross(D, r2) / r2_mag + S)
    
    return v2.reshape(3, 1)

def herrick_gibbs(
        obv1: tuple[np.ndarray, float],
        obv2: tuple[np.ndarray, float],
        obv3: tuple[np.ndarray, float],
        mu: float
    ) -> np.ndarray:
    """The herrick gibbs method of orbit determination.

    Args:
        obv1 (tuple[np.ndarray, float]): The eci position vector and time of the first observation.
        obv2 (tuple[np.ndarray, float]): The eci position vector and time of the second observation.
        obv3 (tuple[np.ndarray, float]): The eci position vector and time of the third observation.
        mu (float): The gravitational parameter of the primary body.

    Returns:
        np.ndarray: The (3 x 1) velocity vector at obv2.
    """
    # TODO: Check vectors are co-planar
    r1, t1 = obv1
    r2, t2 = obv2
    r3, t3 = obv3
    
    r1_mag = np.linalg.norm(r1)
    r2_mag = np.linalg.norm(r2)
    r3_mag = np.linalg.norm(r3)
    
    # Time deltas
    d_t32 = t3 - t2
    d_t21 = t2 - t1
    d_t31 = t3 - t1
    
    A = d_t32 * (1 / (d_t21 * d_t31) + mu / (12 * r1_mag**3))
    B = (d_t32 - d_t21) * (1 / (d_t21 * d_t32) + mu / (12 * r2_mag**3))
    C = d_t21 * (1 / (d_t31 * d_t32) + mu / (12 * r3_mag**3))
    
    v2 = -A * r1 + B * r2 + C * r3
    
    return v2.reshape(3,1)
    
    

if __name__ == "__main__":
    from spacesim import constants as const
    r1 = np.array([[-4036.39024968, -1283.00774237, -5371.32144987]], dtype=np.float64).T
    r2 = np.array([[-3870.33066615, -1688.0824012,  -5322.30251991] ], dtype=np.float64).T
    r3 = np.array([[-3722.74717789, -2048.92475439, -5310.92844261]], dtype=np.float64).T
    
    t_1 = 2.9412000000000003
    t_2 = 3.0 
    t_3 = 3.0588

    
    obv1 = (r1, t_1)
    obv2 = (r2, t_2)
    obv3 = (r3, t_3)
    
    # v2 = gibbs(r1, r2, r3, const.MU_EARTH)
    v2 = herrick_gibbs(obv1, obv2, obv3, const.MU_EARTH)
    print(v2)
    