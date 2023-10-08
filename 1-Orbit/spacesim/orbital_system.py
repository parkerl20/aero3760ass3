from spacesim import celestial_body as cb
from spacesim import orbit as orb


class OrbitalSystem():
    """A class to represent an orbital system.
    
    Attributes:
        primary_body (cb.CelestialBody): The primary body of the orbital system.
        orbits (list[orb.Orbit]): A list of orbits in the orbital system.
    """
    def __init__(self, primary_body: cb.CelestialBody) -> None:
        self.primary_body: cb.CelestialBody = primary_body
        self.orbits: list[orb.Orbit] = []
        
        return
    
    def add_orbit(self, orbit: orb.Orbit) -> None:
        """Adds an orbit to the orbital system.

        Args:
            orbit (orb.Orbit): The orbit to be added.
        """
        self.orbits.append(orbit)
        return