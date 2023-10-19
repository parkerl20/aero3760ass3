from spacesim import orbital_system as os
from spacesim import celestial_body as cb
from spacesim import orbit as orb
from spacesim import orbital_transforms as ot
from spacesim import constants as const
from spacesim import time_util as tu

from matplotlib.figure import Figure
from matplotlib.axes import Axes
from matplotlib import patches
import matplotlib.animation as animation
import matplotlib.pyplot as plt
from pyvista import examples as pv_ex
import datetime as dt
import pyvista as pv
import numpy as np
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
    
    def animate3d(
        self,
        propagation_time: float,
        gif_name: str,
        *,
        fps: int = 30,
        fade_out: bool = False,
        fade_out_length: int = 30,
        animation_step: int = 1
    ) -> None:
        """Provides an animation of the orbital system in 3D.
        
        Hard coded for earth as the primary body.

        Args:
            t (float): The propagation time from each orbits epoch to plot.
            gif_name (str): The name of the gif to save.
        """
        plotter = pv.Plotter(notebook=False, off_screen=True)
        plotter.background_color = 'black'
        cubemap = pv_ex.download_cubemap_space_16k()
        
        earth = pv_ex.planets.load_earth(radius=const.R_EARTH)
        earth_texture = pv_ex.load_globe_texture()
        
        plotter.add_mesh(earth, texture=earth_texture, smooth_shading=True)
        
        # Adjust camera view
        plotter.camera_position = "yz"
        plotter.camera.zoom(1.5)
        original_azimuth = plotter.camera.azimuth + 140
        plotter.camera.azimuth = original_azimuth
        plotter.camera.elevation = -23
        
        earth_rot_offset = 180
        
        # Create gif
        plotter.open_gif(
            gif_name,
            framerate=fps
        )
        
        running = True
        points = [[] for _ in range(len(self.system.orbits))]
        current_datetime = self.system.orbits[0].epoch
        
        opacity_values = np.linspace(0.2, 1.0, fade_out_length)
        
        while running:   
            plotter.clear()
            # plotter.add_actor(cubemap.to_skybox())  
                  
            # Plot orbits
            t_current = -1
            for i, orbit in enumerate(self.system.orbits):
                for _ in range(animation_step):
                    r, _, t = next(orbit)
                    t_current = t
                    points[i].append(r)
                    
                    if fade_out and len(points[i]) > fade_out_length:
                        points[i].pop(0)
                
                r_pos = np.array(points[i])
                trajectory = pv.lines_from_points(r_pos)
                trajectory["opacity"] = opacity_values[-len(r_pos):]
                
                if fade_out:
                    plotter.add_mesh(
                        trajectory,
                        cmap=[orbit.colour],
                        scalars=np.arange(len(r_pos)),
                        opacity="opacity",
                        line_width=3,
                        label=orbit.name,
                        show_scalar_bar=False
                    )
                else:
                    plotter.add_mesh(
                        trajectory,
                        color=orbit.colour,
                        line_width=3,
                        label=orbit.name
                    )

            if t_current > propagation_time:
                running = False
                break
            
            # Rotate earth
            current_datetime = self.system.orbits[0].epoch + dt.timedelta(seconds=int(t_current))
            ERA = np.degrees(tu.gmst(current_datetime) * const.ROT_V_EARTH)
            print(f"t: {t_current}\tERA: {ERA:.2f}")
            
            plotter.add_mesh(earth.rotate_z(ERA + earth_rot_offset), texture=earth_texture, smooth_shading=True)
            plotter.camera.azimuth = original_azimuth + ERA
            
            plotter.write_frame()
        
        plotter.close()
        return
    
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
            
            track_ax.plot(
                lng,
                lat,
                label=orbit.name,
                color=orbit.colour,
                linewidth=0.5
            )
        
        track_ax.legend()
        track_fig.tight_layout()
        
        return track_fig, track_ax
    
    def animate_groundtrack(
        self,
        propagation_time: float,
        gif_name: str,
        *,
        t_start: float = 0,
        map_img: str = None,
        fps: int = 30,
        fade_out: bool = False,
        fade_out_length: int = 30,
        animation_step: int = 1
    ) -> None:
        
        track_fig, track_ax = plt.subplots()
        
        writer = animation.PillowWriter(fps=fps)
        
        if map_img is not None:
            img = plt.imread(map_img)
            track_ax.imshow(img, extent=[-180, 180, -90, 90])
        
        x_ticks = [-180 + x for x in range(0, 380, 30)]
        y_ticks = [-90 + x for x in range(0, 195, 30)]
        
        track_ax.set_xlabel("Longitude", fontname="Times New Roman", fontweight="bold", fontsize=12)
        track_ax.set_ylabel("Latitude", fontname="Times New Roman", fontweight="bold", fontsize=12)
        track_ax.set_xticks(x_ticks)
        track_ax.set_yticks(y_ticks)
        track_fig.tight_layout()
            
        track_ax.grid()
        
        # Create gif
        running = True
        latitudes = [[] for _ in range(len(self.system.orbits))]
        longitudes = [[] for _ in range(len(self.system.orbits))]
        
        with writer.saving(track_fig, gif_name, dpi=300):
            while running:
                lines = []      # store lines to remove later
                
                t_last = -1
                # Plot orbits
                for i, orbit in enumerate(self.system.orbits):
                    for _ in range(animation_step):
                        r, _, t = next(orbit)
                        
                        t_last = t
                        
                        if t > propagation_time:
                            running = False
                            break

                        # Convert eci to lat long
                        r_ecef = ot.ECI_to_ECEF(r, t, orbit.epoch)
                        lat, lng, _ = ot.ECEF_to_LLH(r_ecef).flatten()
                        
                        # Check for discontinuity
                        if len(longitudes[i]) > 0 and abs(longitudes[i][-1] - lng) > 180:
                            longitudes[i].append(np.nan)
                            latitudes[i].append(np.nan)
                        
                        
                        longitudes[i].append(lng)
                        latitudes[i].append(lat)
                        
                        if fade_out and len(longitudes[i]) > fade_out_length:
                            longitudes[i].pop(0)
                            latitudes[i].pop(0)
                            
                    # possibly multiple lines created
                    lines.append(
                        track_ax.plot(
                            longitudes[i],
                            latitudes[i],
                            label=orbit.name,
                            color=orbit.colour,
                            linewidth=1
                        )[0]
                    )
                print(t_last)
                
                # track_ax.legend()                
                # Write frame
                writer.grab_frame()
                
                # Clean up frame
                for line in lines:
                    line.remove()
                
        # writer.finish()
        plt.close(track_fig)
        return
                