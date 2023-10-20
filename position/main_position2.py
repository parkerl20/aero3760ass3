import sys
import os

current_file_path = os.path.abspath(__file__)
parent_directory = os.path.dirname(current_file_path)
grandparent_directory = os.path.dirname(parent_directory)
sys.path.append(grandparent_directory)


from spacesim import constants as const
from spacesim import orbit_determination as od
from spacesim import orbital_transforms as ot
from spacesim import satellite as sat
from spacesim import celestial_body as cb

import numpy as np
import datetime as dt
import matplotlib.pyplot as plt

from .od_simulation import od_simulation

def main_position(
    r: np.ndarray,
    t: np.ndarray,
    epoch: dt.datetime
) -> None:
    """Main function for position determination.
    
    Simulates precise orbit determination for a single satellite

    Args:
        r (np.ndarray): The first three position vectors
            of the satellites in the ECI frame as column vectors
        t (np.ndarray): The times at which the position vectors
            were recorded.
        epoch (dt.datetime): The epoch of the second position vector
    """
    # Parameters
    propagation_time = 1000
    propagation_length = 10         # Prop time per iteration on sat
    max_propagation_step = 1        # Max step size within iteraction interval
    
    # Initial orbit determination
    obv_0 = r[:,0], t[0]
    obv_1 = r[:,1], t[1]
    obv_2 = r[:,2], t[2]
    
    v_1 = od.herrick_gibbs(obv_0, obv_1, obv_2, const.MU_EARTH)
    
    a,e,i,RAAN, arg_p, true_anom = ot.ECI_to_elements(
        r[:,1], v_1, const.MU_EARTH
    ).flatten()
    
    # Temporary values
    # a = 6932386.765062842
    # e = 0.021
    # i = 33
    # RAAN = 58.82
    # arg_p = 180
    # true_anom = 0
    # epoch = dt.datetime(2023, 1, 1)
    
    # Set up satellite
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
        propagation_step=max_propagation_step
    )
    
    od_simulation(
        propagation_time,
        satellite,
    )
    
    return

if __name__ == "__main__":
    main_position(
        None,
        None,
        None
    )