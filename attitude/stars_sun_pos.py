import numpy as np
from matplotlib import pyplot as plt

from skyfield.api import Star, load
from skyfield.data import hipparcos


# Set the location of our satellite, in heliocentric rectangular coordinates
location = np.array([-3.36867753, -10.1950307, -36.55292215])


def build_stereographic_projection(location):
    """Compute *x* and *y* coordinates at which to plot the positions."""

    p = location

    print("p", p)
    u = p / np.linalg.norm(p)
    if len(u.shape) > 1:
        c = u.mean(axis=1)
        c = c / np.linalg.norm(c)
    else:
        c = u
    x_c, y_c, z_c = c

    def project(position):
        print("Project(position)", position)
        p = position.position.au
        u = p / np.sqrt((p * p).sum(axis=0))
        x, y, z = u

        t0 = 1/np.sqrt(x_c**2 + y_c**2)
        t1 = x*x_c
        t2 = np.sqrt(-z_c**2 + 1)
        t3 = t0*t2
        t4 = y*y_c
        t5 = 1/(t1*t3 + t3*t4 + z*z_c + 1)
        t6 = t0*z_c

        return t0*t5*(x*y_c - x_c*y), -t5*(t1*t6 - t2*z + t4*t6)

    return project

def findStarSunPositions(num_stars, day, hour):
    """
    finds sun and num_stars position in metres at given day and hour (of october 2023)

    """
    # Julian Date Time

    ts = load.timescale()
    t= ts.utc(2023, 10, day, hour)
    print("t:", t)

    # An ephemeris from the JPL provides Sun and Earth positions.

    eph = load('de421.bsp')
    sun = eph['sun']
    earth = eph['earth']

    # The Hipparcos mission provides our star catalog.

    with load.open(hipparcos.URL) as f:
        stars = hipparcos.load_dataframe(f)

    # LOCATION DETERMINED HERE --------------------------------------------------------
    # from position determination

    projection = build_stereographic_projection(location)
    field_of_view_degrees = 45.0
    limiting_magnitude = 7.0

    # Now that we have constructed our projection, compute the x and y
    # coordinates that each star will have on the plot.

    star_positions = earth.at(t).observe(Star.from_dataframe(stars))
    sun_position = earth.at(t).observe(sun)
    # print("THE SUN", sun_position.position)


    stars['x'], stars['y'] = projection(star_positions)

    # Create a True/False mask marking the stars bright enough to be
    # included in our plot.  And go ahead and compute how large their
    # markers will be on the plot.

    bright_stars = (stars.magnitude <= limiting_magnitude)

    stars_measured = star_positions[bright_stars][0:num_stars+1].position.m


    return stars_measured,sun_position.position.m

# stars, sun = (findStarSunPositions(5, 22, 1))
