from spacesim import orbit as ob
from spacesim import celestial_body as cb
from spacesim import orbital_system
from spacesim import orbital_system_plotter as osplt
from spacesim import ground_station as gs 
from matplotlib.figure import Figure
from matplotlib.axes import Axes
import matplotlib.pyplot as plt
import numpy as np

import datetime as dt
import matplotlib.pyplot as plt
from pyvista import examples as pv_ex

def main() -> None:
    earth = cb.earth
    my_orbit = ob.Orbit(
        7212100,
        0,
        38,
        0,
        0,
        0,
        earth,
        dt.datetime.now(),
        colour="red"
    )
    earth_system = orbital_system.OrbitalSystem(earth)
    earth_system.add_orbit(my_orbit)
    orb_plotter = osplt.SystemPlotter(earth_system)
    fig, ax = orb_plotter.groundtrack(86400, map_img="./rsc/bluemarble.jpg")


    # orbit = orb_plotter.plot3d(86400)

    # # Add space background
    # cubemap = pv_ex.download_cubemap_space_16k()
    # orbit.add_actor(cubemap.to_skybox())
    
    
    # orbit.show()
    # orbit.show()
    # orb_plotter.
    # plt.show()


    # Groundstation
    marker_lat, marker_lon = -33.89610713284069, 151.19635000670468
    ax.plot(marker_lon, marker_lat, 'go', markersize=11)  

    # plt.show()

    joe = gs.OrbitObservatory("Joe", (-33.89610713284069, 151.19635000670468, 37))
    observe = joe.observe_orbit(my_orbit, 86400)
    visible = observe[4]
    print(f"Total visible time: {calculate_total_time(visible)}")

    # observation = ground_station.observe_orbit(
    #     satellite.orbit,
    #     obs_end,
    #     t_start=obs_start,
    #     analytical=True,
    # )
    # pass_index = 0
    # elev = observe.elevation[pass_index]
    # azim = np.radians(observe.azimuth[pass_index])
    
    # polar_plot(
    #     azim,
    #     elev,
    #     # polar_ax,
    #     arrow_every=100
    # )
    fig, polar_ax = plt.subplots(subplot_kw={"projection": "polar"})
    print(visible)


    polar_ax.set_title('Fly-overs of DEBRA in a day', va='top', fontsize=32)
    for i, period in enumerate(visible):
        elev = observe.elevation[i]
        azim = np.radians(observe.azimuth[i])
        line, = polar_ax.plot(azim, 90 - np.array(elev), label=f"Flyover {i}")  # Plot each flyover and capture the line object
        add_arrow(line, size=20)  # Add arrows to each line

    # Legend and tick size
    polar_ax.legend()
    polar_ax.tick_params(axis='both', labelsize=28)  # Adjust font size for ticks    
    polar_ax.legend(fontsize=21)    

    plt.show()
    return


def calculate_total_time(datetime_pairs):
    total_seconds = sum((end - start).total_seconds() for start, end in datetime_pairs)
    return total_seconds / 60


def add_arrow(line, colour: str = None, size: int = 10, every: int = 10):
    """Adds  arrows to a line on a plot.
    
    Source: https://stackoverflow.com/questions/34017866/arrow-on-a-line-plot-with-matplotlib

    Args:
        ax (Axes): The axes to add the arrows to.
        line (_type_): The line object most recently plotted.
        colour (str, optional): Colour of arrows. Defaults to None.
        size (int, optional): Size of arrows. Defaults to 12.
        every (int, optional): How many points to skip between arrows. Defaults to 10.
    """
    if colour is None:
        colour = line.get_color()
    
    xdata = line.get_xdata()
    ydata = line.get_ydata()
    
    # Calculate points to add arrows to
    points: list[tuple] = []
    dxdy: list[tuple] = []

    for i in range(every, len(xdata) - 1, every):
        points.append((xdata[i], ydata[i]))
        dxdy.append((xdata[i+1], ydata[i+1]))
    
    for p, d in zip(points, dxdy):
        line.axes.annotate('',
            xytext=(p[0], p[1]),
            xy=(d[0], d[1]),
            arrowprops=dict(arrowstyle="-|>", color=colour),
            size=size
        )
    
    return

def polar_plot(theta: np.ndarray, r: np.ndarray, ax: Axes = None, **kwargs) -> tuple[Figure, Axes]:
    """_summary_

    Args:
        theta (np.ndarray): _description_
        r (np.ndarray): _description_
        ax (Axes, optional): Provided axes to plot on.
        If None, a new figure and axes will be created. Defaults to None.

    Returns:
        tuple[Figure, Axes]: _description_
    """
    r_centre = kwargs.setdefault("r_centre", 90)
    r_max = kwargs.setdefault("r_max", 0)
    r_step = kwargs.setdefault("r_step", 10)
    label = kwargs.setdefault("label", None)
    
    arrow_every = kwargs.setdefault("arrow_every", 10)
    arrow_size = kwargs.setdefault("arrow_size", 10)
    
    fig = None
    
    if ax is None:
        fig, ax = plt.subplots(subplot_kw={"projection": "polar"})
        ax.set_rlim(r_centre, r_max, r_step)
    
    lines = ax.plot(theta, r, label=label, linewidth=1.5)
    add_arrow(lines[0], size=arrow_size, every=arrow_every)
    # plt.show()
    return fig, ax

if __name__ == "__main__":
    main()