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
    
    d_RAAN = 90
    shift_RAAN = -5
    h = 274
    a = const.R_EARTH + h * 1000
    e = 0
    i = -35
    omega = 0
    
    orbit_1 = orbit.Orbit(
        a,
        e,
        i,
        shift_RAAN,
        omega,
        0,
        earth,
        dt.datetime.now(),
        "orbit_1",
        colour="red"
    )
    
    orbit_2 = orbit.Orbit(
        a,
        e,
        i,
        shift_RAAN + d_RAAN,
        omega,
        0,
        earth,
        dt.datetime.now(),
        "orbit_2",
        colour="lime"
    )
    
    orbit_system = orbsys.OrbitalSystem(earth)
    orbit_system.add_orbit(orbit_1)
    orbit_system.add_orbit(orbit_2)
    
    k_orbit = 12
    
    prop_time = k_orbit * orbit_1.period
    t_start = (k_orbit - 1) * orbit_1.period
    
    orbit_system_plotter = osplt.SystemPlotter(orbit_system)
    orbit_system_plotter.groundtrack(
        prop_time,
        t_start=t_start,
        map_img="./rsc/bluemarble.jpg"
    )
    
    plt.show()
    return

if __name__ == "__main__":
    main()
    # T = 5400
    # a = (T**2 * const.MU_EARTH / (4 * spconst.pi**2))**(1/3)
    # h = (a - const.R_EARTH) / 1000
    
    # print(f"Altitude: {h} km")