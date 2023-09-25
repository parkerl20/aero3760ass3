"""WARNING: This module is deprecated
"""

from spacesim import orbit as orb
from spacesim import orbital_transforms as ot
from spacesim import celestial_body as cb

from matplotlib.figure import Figure
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.axes import Axes
from matplotlib import patches
import matplotlib.pyplot as plt
import numpy as np

import pyvista as pv
import pyvista.examples as pv_ex


class OrbitPlotter():
    """A class to plot orbits.
    
    Attributes:
        orbit (Orbit): The orbit to be plotted.
        background (tuple[float]): The background colour of the plot. Defaults to black.
    """
    
    def __init__(self, orbit: orb.Orbit) -> None:
        self.orbit = orbit
        self.background = (0,0,0,0.9)     # black
    
    def orbit3D(self, t: float, colour: str = "red", use_km: bool = False,
                t_start: float = 0, plot_body: bool = True) -> tuple[Figure, Axes3D]:
        """Plots the orbit of a satellite in 3D in the ECI frame.

        Args:
            t (float): The end time of the orbit propergation being plot in seconds.
            use_km (bool, optional): When `True` the orbit is presented in km. Defaults to False.
            t_start (float, optional): The start time of the orbit propergation. Defaults to 0.
            plot_body (bool): Determines whether the body being orbitted is plotted
            Defaults to True.
            plot_options (dict, optional): A dictionary of options to customise the plot. Defaults to None.
            plot_earth (bool, optional): Use hardcoded earth function to plot.

        Returns:
            tuple[Figure, Axes3D]: The figure and 3D axes objects with the orbit plotted.
        """
        r, _, _ = self.orbit.propagate(t, t_start, use_km)
        
        fig = plt.figure("3D Orbit Plot")
        ax: Axes3D = fig.add_subplot(111, projection='3d')
        
        if plot_body:
            ax = self.__plot_body3D(ax, use_km)
            
        ax.plot3D(r[0], r[1], r[2], color=colour, zorder=10, label=self.orbit.name)
        ax.set_aspect("equal")
        
        unit = "km" if use_km else "m"
        
        if self.orbit.name is not None:
            ax.set_title(f"{self.orbit.name} Orbit in ECI Frame")
            ax.legend(loc="upper right")
        
        ax.set_xlabel(f"x ({unit})")
        ax.set_ylabel(f"y ({unit})")
        ax.set_zlabel(f"z ({unit})")
        
        ax = OrbitPlotter.__set_background_colour3D(ax, self.background)
        fig.tight_layout()
        
        return fig, ax
          
    def __plot_body3D(self, ax: Axes3D, use_km: bool = False) -> Axes3D:
        """Plots a body on an 3D axes. Assumed to be a planet.

        Args:
            ax (Axes3D): A 3D axes object.

        Returns:
            Axes3D: An axes object with the planet plotted.
        """
        sphere_precision = 40
        theta = np.linspace(0, 2 * np.pi, sphere_precision)
        phi = np.linspace(0, np.pi, sphere_precision)
        
        THETA, PHI = np.meshgrid(theta, phi)
        
        R = self.orbit.body.radius
        R /= 1000 if use_km else 1
        
        x = R * np.cos(THETA) * np.sin(PHI)
        y = R * np.sin(THETA) * np.sin(PHI)
        z = R * np.cos(PHI)
        
        colour = self.orbit.body.colour
        if colour == "blue":
            cmap = "Blues"
        
        ax.plot_surface(x, y, z, cmap=cmap, zorder=2)
        # Fake point to allow legend label
        ax.scatter(0, 0, 0, color=colour, edgecolors=None, label=self.orbit.body.name)
        
        return ax
    
    # Hard coded function that plots a 3D earth
    def __plot_earth(self, ax: Axes3D, use_km: bool = False) -> Axes3D:
        image_file = "./rsc/earth.jpg"
        img = plt.imread(image_file)

        # define a grid matching the map size, subsample along with pixels
        theta = np.linspace(0, np.pi, img.shape[0])
        phi = np.linspace(0, 2*np.pi, img.shape[1])

        count = 180 # keep 180 points along theta and phi
        theta_inds = np.linspace(0, img.shape[0] - 1, count).round().astype(int)
        phi_inds = np.linspace(0, img.shape[1] - 1, count).round().astype(int)
        theta = theta[theta_inds]
        phi = phi[phi_inds]
        img = img[np.ix_(theta_inds, phi_inds)]

        theta,phi = np.meshgrid(theta, phi)
        R = self.orbit.body.radius
        R /= 1000 if use_km else 1
        

        # sphere
        x = R * np.sin(theta) * np.cos(phi)
        y = R * np.sin(theta) * np.sin(phi)
        z = R * np.cos(theta)

        # create 3d Axes
        ax.plot_surface(x.T, y.T, z.T, facecolors=img/255, cstride=1, rstride=1) # we've already pruned ourselves

        # make the plot more spherical
        ax.axis('scaled')
        return ax
    
    def orbit3D_nice(self, t: float, colour: str = "red", use_km: bool = False,
                     t_start: float = 0, plot_body: bool = True) -> pv.Plotter:
        """Creates a 3D plot of the orbit in PyVista

        Args:
            t (float): The end time of the orbit propergation.
            colour (str, optional): Trajectory colour. Defaults to "red".
            use_km (bool, optional): Whether to use km or m. Defaults to False.
            t_start (float, optional): The start time of the orbit propergation. Defaults to 0.
            plot_body (bool, optional): Whether to plot the orbit's body. Defaults to True.

        Returns:
            tuple[Figure, Axes3D]: _description_
        """
        r,_,_ = self.orbit.propagate(t, t_start, use_km)
        
        pl = pv.Plotter()
        
        cubemap = pv_ex.download_cubemap_space_16k()
        pl.add_actor(cubemap.to_skybox())
        
        trajectory = pv.Spline(r.T)   # Splines for 3D lines
        pl.add_mesh(trajectory, color=colour, line_width=4, label=self.orbit.name)
        
        if plot_body:
            # Use Earth in example package, currently hard coded
            R = self.orbit.body.radius
            R /= 1000 if use_km else 1
            earth = pv_ex.planets.load_earth(radius=R)
            earth_texture = pv_ex.load_globe_texture()
            
            pl.add_mesh(earth, texture=earth_texture, smooth_shading=True)
        
        # pl.add_camera_orientation_widget()      # adds tool to help orientation
        pl.add_axes()
        pl.add_legend(face="r")
        
        return pl
        
    
    def orbit2D(self, t: float, colour: str = "red", use_km: bool = False,
                   t_start: float = 0, plot_body: bool = True) -> tuple[Figure, Axes]:
        """_summary_

        Args:
            t (float): _description_
            colour (str, optional): _description_. Defaults to "red".
            use_km (bool, optional): _description_. Defaults to False.
            t_start (float, optional): _description_. Defaults to 0.
            plot_body (bool, optional): _description_. Defaults to True.

        Returns:
            tuple[Figure, Axes]: The figure and 2D axes objects with the orbit plotted.
        """
        r, _, _= self.orbit.propagate(t, t_start, use_km, "perifocal")
        
        fig = plt.figure("2D Orbit Plot")
        ax = fig.add_subplot()
        ax.set_aspect("equal", adjustable="box")
        
        R = self.orbit.body.radius
        R /= 1000 if use_km else 1
        
        if plot_body:
            body_colour = OrbitPlotter.__nicer_colour(self.orbit.body.colour)
            
            body_circle = patches.Circle((0,0), R, alpha = 0.8, color=body_colour)
            ax.add_patch(body_circle)
            ax.scatter(0,0, color=body_colour, alpha=0.8, edgecolors=None, label=self.orbit.body.name)
            
        ax.plot(r[0], r[1], color=colour, label=self.orbit.name)
        
        if self.orbit.name is not None:
            ax.set_title(f"{self.orbit.name} Orbit in Perifocal Frame",
                         fontname="Times New Roman", fontweight="bold", fontsize=16)
            ax.legend(loc="upper right")
        
        unit = "km" if use_km else "m"
        ax.set_xlabel(f"p ({unit})", fontname="Times New Roman", fontweight="bold", fontsize=12)
        ax.set_ylabel(f"q ({unit})", fontname="Times New Roman", fontweight="bold", fontsize=12)
        
        ax.grid()
        ax.set_aspect("equal", adjustable="box")
        ax.set_facecolor("black")
        fig.tight_layout()
        
        return fig, ax
    
    @staticmethod
    def __set_background_colour3D(ax: Axes3D, rgba: tuple[float]) -> Axes3D:
        """Sets the background of the provided 3D axes.

        Args:
            ax (Axes3D): A 3D axes object.
            rgba (tuple[float], optional): The colour of the background.
        Returns:
            Axes3D: An axes object with a black background.
        """
        ax.w_xaxis.set_pane_color(rgba)
        ax.w_yaxis.set_pane_color(rgba)
        ax.w_zaxis.set_pane_color(rgba)
        return ax
            
    def groundtrack(self, t: int, t_start: int = 0, map_img: str = None,
                    colour: str = "red") -> tuple[Figure, Axes3D]:
        """_summary_

        Args:
            epoch (dt.datetime): _description_
            map_img (str, optional): _description_. Defaults to None.

        Returns:
            tuple[Figure, Axes3D]: _description_
        """
        r_eci, _, time_eci = self.orbit.propagate(t, t_start)
        lat, lng = [], []
        
        # Convert eci to latitude, longitude
        for i in range(len(time_eci)):
            r_ecef = ot.ECI_to_ECEF(r_eci[:,i], time_eci[i], self.orbit.epoch)
            lat_i, lng_i, _ = ot.ECEF_to_LLH(r_ecef).flatten()
            
            # Remove discontinuities
            if i > 0 and abs(lng[-1] - lng_i) > 180:
                lng.append(np.nan)
                lat.append(np.nan)
            
            lng.append(lng_i)
            lat.append(lat_i)
        
        # Plot groundtrack
        fig = plt.figure("Groundtrack")
        ax = fig.add_subplot()
        
        if map_img is not None:
            img = plt.imread(map_img)
            ax.imshow(img, extent=[-180, 180, -90, 90])
            
        ax.plot(lng, lat, color='red', label=self.orbit.name, linewidth=0.3)
        
        x_ticks = [-180 + x for x in range(0, 380, 30)]
        y_ticks = [-90 + x for x in range(0, 195, 30)]
        
        ax.set_xlabel("Longitude", fontname="Times New Roman", fontweight="bold", fontsize=12)
        ax.set_ylabel("Latitude", fontname="Times New Roman", fontweight="bold", fontsize=12)
        ax.set_xticks(x_ticks)
        ax.set_yticks(y_ticks)
        
        if self.orbit.name is not None:
            ax.set_title(f"{self.orbit.name} Groundtrack",
                         fontname="Times New Roman", fontweight="bold", fontsize=16)
            ax.legend(loc="upper right")
        
        ax.grid()
        fig.tight_layout()
        
        return fig, ax

    @staticmethod
    def __nicer_colour(colour: str) -> str:
        """Provides better alternatives for standard matplotlib colours.

        Args:
            colour (str): matplotlib colour.

        Returns:
            str: Hex code for the colour.
        """
        colour_hex = "#ffffff"
        
        if colour == "blue":
            colour_hex ="#4682b4"       # Light blue
        
        return colour_hex
    
    @staticmethod
    def __filter_visible(r: np.ndarray, elev: float, azim: float) -> np.ndarray:
        """
        DEPRECATED: Used for plotting with matplotlib in 3D.
        
        Filters the points out that are behind a 3D sphere given the viewing
        elevation and azimuth.
        
        Source: https://stackoverflow.com/questions/41699494/how-to-obscure-a-line-behind-a-surface-plot-in-matplotlib

        Args:
            r (np.ndarray): A (3 x n) matrix of position vectors.
            elev (float): The viewing elevation in degrees.
            azim (float): The viewing azimuth in degrees.

        Returns:
            np.ndarray: A (3 x k) matrix of position vectors where k <= n.
        """
        r_copy = r.copy()
        # convert to radians
        a = azim * np.pi / 180.0 - np.pi
        e = elev * np.pi / 180.0 - np.pi / 2.
        # normal vector to viewing plane
        X = np.array([ np.sin(e) * np.cos(a), np.sin(e) * np.sin(a), np.cos(e)])
        cond = (X @ r) < 0     # filter
        
        r_copy[:, cond] = np.nan
        
        return r_copy