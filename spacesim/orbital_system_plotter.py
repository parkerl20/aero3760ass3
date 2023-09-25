from spacesim import orbital_system as os
from spacesim import celestial_body as cb
from spacesim import orbit as orb

from matplotlib.figure import Figure
from matplotlib.axes import Axes
from matplotlib import patches
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
        max_step: float = 20
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
            r, _, _ = orbit.propagate(t, t_start, use_km, max_step=max_step)
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
        pass