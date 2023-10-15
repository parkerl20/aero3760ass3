"""Previous code to read TLE's, plot orbits and plot ground tracks

"""

from spacesim import celestial_body as cb
from spacesim import satellite as sat
from spacesim import orbital_system as orbsys
from spacesim import orbital_system_plotter as osplt
from spacesim import constants as const
from spacesim import ground_station as gs
from spacesim import orbital_transforms as ot
from spacesim import orbit as orb

from pyvista import examples as pv_ex
import numpy as np
import matplotlib.pyplot as plt

def perturbation_dynamics(t: float, state_vector: np.ndarray, orbit: orb.Orbit) -> np.ndarray:
    # Satellite dynamics that accounts for perturbations
    def perturb_acceleration(r: np.ndarray) -> np.ndarray:
        # Inner function to calculate perturbation accleration
        x, y, z = r
        r_mag = np.sqrt(r.dot(r))
        
        P = - (3 * mu * J2 * R**2) / (2 * r_mag**5)
        a_x = (1 - 3 * z**2 / r_mag**2) * x
        a_y = (1 - 3 * z**2 / r_mag**2) * y
        a_z = (3 - 3 * z**2 / r_mag**2) * z
        
        return P * np.array([a_x, a_y, a_z], dtype=np.float64)
    
    
    mu = orbit.body.gravitational_parameter
    J2 = orbit.body.J2
    R = orbit.body.radius
    
    r = state_vector[0:3]
    v = state_vector[3:6]
    
    r_mag = np.sqrt(r.dot(r))
    
    # Dynamics
    P = perturb_acceleration(r)
    
    r_dot = v
    v_dot = -(mu / (r_mag ** 3)) * r + P    # two-body equation

    return [*r_dot, *v_dot]

def main() -> None:
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
    
    earth_orbital_system = orbsys.OrbitalSystem(earth)
    
    # Params
    propagation_time = 60 * 60 * 24
    
    #---------------- Create single satellite
    tle_file = "./rsc/TLE/navstar_43.txt"
    satellite = None
    
    with open(tle_file, "r") as f:
        name = f.readline().strip()
        tle_line_1 = f.readline().strip()
        tle_line_2 = f.readline().strip()
        
        satellite = sat.Satellite(name, tle_line_1, tle_line_2, earth)
        print(satellite)
        
    earth_orbital_system.add_orbit(satellite.orbit)    
    
    
    # ---------------- Ground station demo
    sydney_location = (-33.8688, 151.2093, 0)
    sydney = gs.OrbitObservatory("Sydney", sydney_location)
    
    sydney_obs = sydney.observe_orbit(
        satellite.orbit,
        propagation_time,
        analytical=True
    )

    start = sydney_obs.visibility_period[0][0]
    end = sydney_obs.visibility_period[0][1]
    print(end - start)
    
    #---------------- Create multiple satellites
    tle_files = [
        "./rsc/TLE/navstar_43.txt",
        "./rsc/TLE/navstar_47.txt",
        "./rsc/TLE/navstar_51.txt",
        "./rsc/TLE/navstar_52.txt",
    ]
    
    satellites: list[sat.Satellite] = []
    
    for tle in tle_files:
        with open(tle, "r") as f:
            name = f.readline().strip()
            tle_line_1 = f.readline().strip()
            tle_line_2 = f.readline().strip()
            
            satellite = sat.Satellite(name, tle_line_1, tle_line_2, earth)
            satellites.append(satellite)
            
            earth_orbital_system.add_orbit(satellite.orbit)
    
    
    # ---------------- Plot orbits
    system_plotter = osplt.SystemPlotter(earth_orbital_system)
    
    gt_fig, gt_ax = system_plotter.groundtrack(
        propagation_time,
        map_img="./rsc/bluemarble.jpg"
    )
    
    plt.show()
    
    earth_pl = system_plotter.plot3d(
        propagation_time
    )
    
    # Add space background
    cubemap = pv_ex.download_cubemap_space_16k()
    earth_pl.add_actor(cubemap.to_skybox())
    
    
    earth_pl.show()
    

    return
    


if __name__ == "__main__":
    main()