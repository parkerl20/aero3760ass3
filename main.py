"""
This will be the main code file that is used to run the entire program with all crucial plots. 
In particular, this helps with integration as results can be sent from section to section 
easily, like ECI data from orbitt getting sent to gee to plot swathe width 
"""

from orbit.mainOrbit import mainOrbit
from position.main_position2 import main_position
from gee.mainGee import mainGee
import datetime as dt


def main():
    # Runs the main orbit code with results being r, v, t of the 4 satellites
    results = mainOrbit(show_results=1)
    
    # ---------------- Orbit determination setup
    # First three position vectors, used for IOD
    # r_obv = results[0]['r'][:]
    # t_obv = results[0]['t'][:]
    # epoch = dt.datetime(2023, 1, 1)

    # main_position(r_obv, t_obv, epoch)

    # Runs the main GEE code opening a map in your browser
    mainGee(results)

    return 0


if __name__ == "__main__":
    main()