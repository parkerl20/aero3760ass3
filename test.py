from spacesim import constants as const
from scipy import constants as spconst
from spacesim import orbit
from spacesim import celestial_body as cb
from spacesim import orbital_system as orbsys
from spacesim import orbital_system_plotter as osplt
import numpy as np
from pyvista import examples as pv_ex
import pyvista as pv
import matplotlib.pyplot as plt
import datetime as dt

def main() -> None:
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
    
    d_RAAN = 60
    shift_RAAN = -40
    h = 658.899283   # km
    a = const.R_EARTH + h * 1000
    e = 0
    i_1 = 90
    i_2 = i_1
    omega = 40
    date = dt.datetime(2023, 10, 1)
    theta_1 = 90
    theta_2 = theta_1
    
    orbit_1 = orbit.Orbit(
        a,
        e,
        i_1,
        shift_RAAN,
        omega,
        theta_1,
        earth,
        date,
        "orbit_1",
        colour="red"
    )
    
    orbit_2 = orbit.Orbit(
        a,
        e,
        i_2,
        shift_RAAN + d_RAAN,
        omega,
        theta_2,
        earth,
        date,
        "orbit_2",
        colour="lime"
    )
    
    orbit_3 = orbit.Orbit(
        a,
        e,
        i_2,
        shift_RAAN + 2 * d_RAAN,
        omega,
        theta_2,
        earth,
        date,
        "orbit_3",
        colour="yellow"
    )
    
    orbit_system = orbsys.OrbitalSystem(earth)
    orbit_system.add_orbit(orbit_1)
    orbit_system.add_orbit(orbit_2)
    orbit_system.add_orbit(orbit_3)
    
    k = 1
    
    delta = (spconst.pi / 3) / const.ROT_V_EARTH
    
    prop_time = k * delta
    t_start = prop_time - (spconst.pi / 8) / const.ROT_V_EARTH
    # prop_time = 1 * 24 * 60 * 60
    # t_start = 0
    
    orbit_system_plotter = osplt.SystemPlotter(orbit_system)
    orbit_system_plotter.groundtrack(
        prop_time,
        t_start=t_start,
        map_img="./rsc/bluemarble.jpg"
    )
    
    plt.show()
    
    # Add space background
    # earth_pl = orbit_system_plotter.plot3d(
    #     orbit_1.period
    # )
    # cubemap = pv_ex.download_cubemap_space_16k()
    # earth_pl.add_actor(cubemap.to_skybox())
    
    
    # earth_pl.show()
    return

if __name__ == "__main__":
    # main()
    # T = 5760
    
    t = (2 *spconst.pi) / const.ROT_V_EARTH
    T = t / 15
    print(T)
    
    a = (T**2 * const.MU_EARTH / (4 * spconst.pi**2))**(1/3)
    h = (a - const.R_EARTH) / 1000
    
    print(f"Altitude: {h} km")