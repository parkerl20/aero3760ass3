from spacesim import constants as const
import scipy.constants as spconst
import pyvista as pv


class CelestialBody():
    """A class to represent a celestial body.
    
    Attributes:
        name (str): The name of the celestial body.
        mass (float): The mass of the celestial body in kilograms.
        radius (float): The radius of the celestial body in meters.
        J2 (float): The J2 parameter of the celestial body.
        colour (str): The colour of the celestial body used for 2D plotting.
        body_mesh (pv.PolyData): A mesh of the celestial body used for 3D plotting.
        body_texture (pv.Texture): The texture of the 3D mesh.
    """
    
    def __init__(
        self,
        name: str,
        mass: float,
        radius: float = None,
        J2: float = None,
        colour: str = None,
        body_mesh: pv.PolyData = None,
        body_texture: pv.Texture = None
    ) -> None:
        """Initializes a CelestialBody object.

        Args:
            name (str): The name of the celestial body.
            mass (float): The mass of the celestial body in kilograms.
            radius (float, optional): The radius of the celestial body 
            in meters. Defaults to None.
            J2 (float, optional): The J2 parameter of the celestial body.
            Defaults to None.
            colour (str, optional): The colour of the celestial body used
            for 2D plotting. Defaults to None.
            body_mesh (pv.PolyData, optional): A mesh of the celestial
            body used for 3D plotting. Defaults to None.
            body_texture (pv.Texture, optional): The texture of the 3D mesh.
            Defaults to None.
        """
        self.name: str = name
        self.mass: float = mass
        self.radius: float = radius
        self.J2: float = J2
        self.gravitational_parameter: float = spconst.G * mass      # m^3 / s^2
        
        # Plotting attributes
        self.colour: str = colour
        self.body_mesh: pv.PolyData = body_mesh
        self.body_texture: pv.Texture = body_texture
        
        return

# Default celestial bodies
earth = CelestialBody(
    "Earth",
    const.M_EARTH,
    const.R_EARTH,
    J2 = const.J2_EARTH,
    colour="blue"
)