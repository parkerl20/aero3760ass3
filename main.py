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
from attitude.attitude_nlls import run_attitude_determ
from remote_sensing import data_accuracy


def main():
    # Uncomment the second line to run through the calculations instead of just results
    satellite_simulation(run_sim=0)
    # satellite_simulation(run_sim=1)

    # length_flag = 1 for short version, length_flag = 2 for long version
    run_attitude_determ(length_flag = 1)
    # run_attitude_determ(length_flag = 2))


    return 0


def satellite_simulation(run_sim):
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
    
    mapping_error = 2 # From attitude
    mainGee(results, mapping_error, run_sim)
    
    data_accuracy.main()

if __name__ == "__main__":
    main()