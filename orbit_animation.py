from spacesim import satellite as sat
from spacesim import orbital_system as osys
from spacesim import orbital_system_plotter as osplt
from spacesim import celestial_body as cb
from spacesim import constants as const

from pyvista import examples as pv_ex
import pyvista as pv
import datetime as dt

def main() -> None:
    a = 6932386.57371916
    e = 0
    i = -33
    arg_p = 0
    RAAN = 238.82
    theta = 0
    epoch = dt.datetime(2023, 1, 1)
    
    propagation_time = 6 * 60 * 60
    propagation_length = 40
    propagation_step = propagation_length / 2
    
    realtime_sat_1 = sat.RealTimeSatellite(
        a,
        e,
        i,
        RAAN,
        arg_p,
        theta,
        cb.earth,
        epoch,
        name="Real Time",
        colour="red",
        propagation_length=propagation_length,
        propagation_step=propagation_step
    )
    
    realtime_sat_2 = sat.RealTimeSatellite(
        a,
        e,
        i,
        RAAN + 90,
        arg_p,
        theta,
        cb.earth,
        epoch,
        name="Real Time",
        colour="lime",
        propagation_length=propagation_length,
        propagation_step=propagation_step
    )
    
    realtime_sat_3 = sat.RealTimeSatellite(
        a,
        e,
        i,
        RAAN + 180,
        arg_p,
        theta,
        cb.earth,
        epoch,
        name="Real Time",
        colour="yellow",
        propagation_length=propagation_length,
        propagation_step=propagation_step
    )
    
    realtime_sat_4 = sat.RealTimeSatellite(
        a,
        e,
        i,
        RAAN + 270,
        arg_p,
        theta,
        cb.earth,
        epoch,
        name="Real Time",
        colour="blue",
        propagation_length=propagation_length,
        propagation_step=propagation_step
    )
    
    earth_mesh = pv_ex.planets.load_earth(radius=const.R_EARTH)
    earth_tex = pv_ex.load_globe_texture()
    earth = cb.earth
    earth.body_mesh = earth_mesh
    earth.body_texture = earth_tex
    
    
    earth_orbital_system = osys.OrbitalSystem(cb.earth)
    earth_orbital_system.add_orbit(realtime_sat_1)
    earth_orbital_system.add_orbit(realtime_sat_2)
    earth_orbital_system.add_orbit(realtime_sat_3)
    earth_orbital_system.add_orbit(realtime_sat_4)
    
    os_plotter = osplt.SystemPlotter(earth_orbital_system)
    
    plotter = os_plotter.animate3d(
        propagation_time,
        "orbit_animation_black.gif",
        fps=10,
        fade_out=True,
        animation_step=3
    )
    
    realtime_sat_1.reset()
    realtime_sat_2.reset()
    realtime_sat_3.reset()
    realtime_sat_4.reset()
    
    os_plotter.animate_groundtrack(
        propagation_time,
        "orbit_animation_groundtrack2.gif",
        fps=10,
        map_img="./rsc/bluemarble.jpg",
        fade_out=True,
        fade_out_length=40,
        animation_step=3
    )
    
    return

if __name__ == "__main__":
    main()