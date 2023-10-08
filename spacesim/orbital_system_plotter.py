from spacesim import orbital_system as os
from spacesim import celestial_body as cb
from spacesim import orbit as orb
from spacesim import orbital_transforms as ot

from matplotlib.figure import Figure
from matplotlib.axes import Axes
from matplotlib import patches
import numpy as np
import matplotlib.pyplot as plt
import pyvista as pv
import warnings


class SystemPlotter():
    def __init__(self, system: os.OrbitalSystem) -> None:
        self.system = system
    
    def plot2d(
        self,
        t: float,
        *,
        t_start: float = 0,
        use_km: bool = False
    ) -> tuple[Figure, Axes]:
        """Plots the orbital system in 2D.

        Args:
            t (float): The propagation time from each orbits epoch to plot.
            t_start (float, optional): The start time of the plot. Defaults to 0.
            use_km (bool, optional): Whether to use km or m as the unit of distance.
            Defaults to metres.
        
        Returns:
            tuple[Figure, Axes]: The figure and axes of the plot.
        """
        if len(self.system.orbits) > 1:
            warnings.warn("Plotting in 2D is not ideal for a system with\
                multiple orbits if orbits are on different orbital planes.")
        
        coord_frame = "perifocal"
        
        # Plot body
        system_fig, system_ax = plt.subplots()
        
        body_radius = self.system.primary_body.radius
        body_radius = body_radius / 1000 if use_km else body_radius

        body_circle = patches.Circle((0,0), body_radius, color=self.system.primary_body.colour)
        system_ax.add_patch(body_circle)
        
        # Plot orbits
        for orbit in self.system.orbits:
            r, _, _ = orbit.propagate(t, t_start, use_km, coord_frame, max_step=20)
            system_ax.plot(r[0], r[1], color=orbit.body.colour, label=orbit.body.name)
        
        # Plot formatting
        system_ax.set_title(f"{self.system.name} Orbital System")
        
        unit = "km" if use_km else "m"
        system_ax.set_xlabel(f"p ({unit})", fontname="Times New Roman", fontweight="bold", fontsize=12)
        system_ax.set_ylabel(f"q ({unit})", fontname="Times New Roman", fontweight="bold", fontsize=12)
        
        system_ax.grid()
        system_ax.set_aspect("equal", adjustable="box")
        system_ax.set_facecolor("black")
        system_fig.tight_layout()
        
        return system_fig, system_ax
    
    def plot3d(
        self,
        t: float,
        *,
        t_start: float = 0,
        use_km: bool = False,
        max_step: float = 20,
        analytical: bool = False
    ) -> pv.Plotter:
        """Plots the orbital system in 3D.

        Args:
            t (float): The propagation time from each orbits epoch to plot.
            t_start (float, optional): The start time of the plot. Defaults to 0.
            use_km (bool, optional): Whether to use km or m as the unit of distance.
        """
        system_pl = pv.Plotter()
        
        # Plot body
        system_pl.add_mesh(self.system.primary_body.body_mesh, texture=self.system.primary_body.body_texture)
        
        # Plot orbits
        for orbit in self.system.orbits:
            r, _, _ = orbit.propagate(
                t,
                t_start=t_start,
                use_km=use_km,
                max_step=max_step,
                analytical=analytical
            )
            
            trajectory = pv.Spline(r.T)   # Splines for 3D lines
            system_pl.add_mesh(trajectory, color=orbit.colour, line_width=2, label=orbit.name)
        
        return system_pl
    
    def groundtrack(
        self,
        t: float,
        *,
        t_start: float = 0,
        map_img: str = None,
    ) -> None:
        """Plots the groundtrack of the orbital system.

        Args:
            t (float): The propagation time from each orbits epoch to plot.
            t_start (float, optional): The start time of the plot. Defaults to 0.
            map_img (str, optional): The path to the map image to use. Defaults to None.
        """
        # Set up figure
        track_fig = plt.figure("Groundtracks")
        track_ax = track_fig.add_subplot()
        
        if map_img is not None:
            img = plt.imread(map_img)
            track_ax.imshow(img, extent=[-180, 180, -90, 90])
        
        x_ticks = [-180 + x for x in range(0, 380, 30)]
        y_ticks = [-90 + x for x in range(0, 195, 30)]
        
        track_ax.set_xlabel("Longitude", fontname="Times New Roman", fontweight="bold", fontsize=12)
        track_ax.set_ylabel("Latitude", fontname="Times New Roman", fontweight="bold", fontsize=12)
        track_ax.set_xticks(x_ticks)
        track_ax.set_yticks(y_ticks)
            
        track_ax.grid()

        # Add orbits
        for orbit in self.system.orbits:
            r_eci, _, t_eci = orbit.propagate(t, t_start=t_start)
            lat, lng = [], []
            
            # Convert eci to lat long
            for i in range(len(t_eci)):
                r_ecef = ot.ECI_to_ECEF(r_eci[:, i], t_eci[i], orbit.epoch)
                lat_i, lng_i, _ = ot.ECEF_to_LLH(r_ecef).flatten()
                
                if i > 0 and abs(lng[-1] - lng_i) > 180:
                    lng.append(np.nan)
                    lat.append(np.nan)
                
                lng.append(lng_i)
                lat.append(lat_i)
            
            track_ax.plot(lng, lat, label=orbit.name, linewidth=0.5)
        
        track_ax.legend()
        track_fig.tight_layout()
        
        return track_fig, track_ax