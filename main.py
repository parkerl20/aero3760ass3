"""
This will be the main code file that is used to run the entire program with all crucial plots. 
In particular, this helps with integration as results can be sent from section to section 
easily, like ECI data from orbitt getting sent to gee to plot swathe width 
"""

from orbit.mainOrbit import mainOrbit
from gee.mainGee import mainGee


def main():
    # Runs the main orbit code with results being r, v, t of the 4 satellites
    results = mainOrbit()

    # Runs the main GEE code opening a map in your browser
    mainGee()

    return 0


if __name__ == "__main__":
    main()