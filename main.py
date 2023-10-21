"""
This will be the main code file that is used to run the entire program with all crucial plots. 
In particular, this helps with integration as results can be sent from section to section 
easily, like ECI data from orbitt getting sent to gee to plot swathe width 
"""

from orbit.mainOrbit import mainOrbit
from orbit_simulation import orbit_simulation
from gee.mainGee import mainGee
import datetime as dt


def main():
    # Runs the main orbit code with results being r, v, t of the 4 satellites
    results = mainOrbit()
    
    # ---------------- Satellite simulation
    # Temporary values
    a = 6932386.765062842
    e = 0.021
    i = 33
    RAAN = 58.82
    arg_p = 180
    true_anom = 0
    epoch = dt.datetime(2023, 1, 1)
    simulation_time = 1000

    orbit_simulation(
        simulation_time,
        a,
        e,
        i,
        RAAN,
        arg_p,
        true_anom,
        epoch
    )

    # Runs the main GEE code opening a map in your browser
    # mainGee()

    return 0


if __name__ == "__main__":
    main()