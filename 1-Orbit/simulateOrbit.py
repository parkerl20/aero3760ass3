# -------- Add spacesim to script path
import sys
import os

current_file_path = os.path.abspath(__file__)
parent_directory = os.path.dirname(current_file_path)
grandparent_directory = os.path.dirname(parent_directory)
sys.path.append(grandparent_directory)

from spacesim import celestial_body as cb
from spacesim import satellite as sat
from spacesim import orbital_system as orbsys
from spacesim import orbital_system_plotter as osplt
from spacesim import constants as const
from spacesim import ground_station as gs
from spacesim import orbital_transforms as ot
from spacesim import orbit as orb
from datetime import datetime

from pyvista import examples as pv_ex
import numpy as np
import matplotlib.pyplot as plt

def simulateOrbit(a, e, i, rt_asc, arg_p, theta) -> None:
    #---------------- Create earth object
    earth_mesh = pv_ex.planets.load_earth(radius=const.R_EARTH)
    earth_tex = pv_ex.load_globe_texture()
    
    earth = cb.CelestialBody(
        "Earth",
        const.M_EARTH,
        const.R_EARTH,
        J2 = const.J2_EARTH,
        colour="blue",
        body_mesh=earth_mesh,
        body_texture=earth_tex,
    )

    # Specify the epoch (start time) for the orbit
    epoch = datetime(2023, 1, 1)
    
    # Create orbit
    orbit = orb.Orbit(a, e, i, rt_asc, arg_p, theta, earth, epoch, "1 satellite configuration")
    
    # Params
    propagation_time = 60 * 60 * 24

    # Results
    result_r, result_v, result_t = orbit.propagate(propagation_time)

    # Adding the orbit to a system
    earth_orbital_system = orbsys.OrbitalSystem(earth)
    earth_orbital_system.add_orbit(orbit)
    
    # System plotter
    plotter = osplt.SystemPlotter(earth_orbital_system)

    # Plot groundtrack
    gt_fig, gt_ax = plotter.groundtrack(
        propagation_time,
        map_img="./rsc/bluemarble.jpg"
    )
    plotter.plot3d(propagation_time).show()
    
    plt.show()

    return result_r, result_v, result_t