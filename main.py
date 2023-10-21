"""
This will be the main code file that is used to run the entire program with all crucial plots. 
In particular, this helps with integration as results can be sent from section to section 
easily, like ECI data from orbitt getting sent to gee to plot swathe width 
"""
from spacesim import orbital_transforms as ot
from spacesim import constants as const
from orbit.mainOrbit import mainOrbit
from orbit_simulation import orbit_simulation
from gee.mainGee import mainGee
import datetime as dt


def main():
    # Runs the main orbit code with results being r, v, t of the 4 satellites
    results = mainOrbit(show_results=1)
    
    # ---------------- Satellite simulation
    r0 = results[0]['r'][:,0]
    v0 = results[0]['v'][:,0]
    
    a,e,i,RAAN, arg_p, true_anom = ot.ECI_to_elements(r0, v0, const.MU_EARTH).flatten()

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
    # ---------------- Orbit determination setup
    # First three position vectors, used for IOD
    # r_obv = results[0]['r'][:]
    # t_obv = results[0]['t'][:]
    # epoch = dt.datetime(2023, 1, 1)

    # main_position(r_obv, t_obv, epoch)

    # Runs the main GEE code opening a map in your browser
    # mainGee(results)

    return 0


if __name__ == "__main__":
    main()