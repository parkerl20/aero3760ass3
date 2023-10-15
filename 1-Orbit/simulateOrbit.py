# -------- Add spacesim to script path
import sys
import os

current_file_path = os.path.abspath(__file__)
parent_directory = os.path.dirname(current_file_path)
grandparent_directory = os.path.dirname(parent_directory)
sys.path.append(grandparent_directory)

from spacesim import celestial_body as cb
from spacesim import orbital_system as orbsys
from spacesim import orbital_system_plotter as osplt
from spacesim import constants as const
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
    orbit1 = orb.Orbit(a, e, i, rt_asc, arg_p, theta, earth, epoch, "1st Satellite")
    orbit2 = orb.Orbit(a, e, i, rt_asc+90.0, arg_p, theta, earth, epoch, "2nd Satellite")
    orbit3 = orb.Orbit(a, e, i, rt_asc+180.0, arg_p, theta, earth, epoch, "3rd Satellite")
    orbit4 = orb.Orbit(a, e, i, rt_asc+270.0, arg_p, theta, earth, epoch, "4th Satellite")
    
    # Params
    propagation_time = 60 * 60 * 60

    # Results
    # result_r1, result_v1, result_t1 = orbit1.propagate(propagation_time)

    # Adding the orbit to a system
    earth_orbital_system = orbsys.OrbitalSystem(earth)
    earth_orbital_system.add_orbit(orbit1)
    earth_orbital_system.add_orbit(orbit2)
    earth_orbital_system.add_orbit(orbit3)
    earth_orbital_system.add_orbit(orbit4)
    
    # System plotter
    plotter = osplt.SystemPlotter(earth_orbital_system)

    # Plot groundtrack
    gt_fig, gt_ax = plotter.groundtrack(
        propagation_time,
        map_img="./rsc/bluemarble.jpg"
    )

    # Plot 3D orbit
    earth_pl = plotter.plot3d(
        propagation_time
    )

    # Add space background
    cubemap = pv_ex.download_cubemap_space_16k()
    earth_pl.add_actor(cubemap.to_skybox())
    earth_pl.show()
    
    plt.show()

    return 0 #result_r, result_v, result_t