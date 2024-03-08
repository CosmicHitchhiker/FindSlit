import numpy as np
from matplotlib.patches import Polygon
from astropy import units as u

class mySlit(Polygon):
    """
    Create a patch representing a latitude-longitude quadrangle rotated by
    given angle.

    See astropy.visualization.wcsaxes.Quadrangle for details.

    Parameters
    ----------
    anchor : tuple or `~astropy.units.Quantity` ['angle']
        Center of slit.
        This can be either a tuple of two `~astropy.units.Quantity` objects, or
        a single `~astropy.units.Quantity` array with two elements.
    width : `~astropy.units.Quantity` ['angle'], default=1.5*u.arcsec
        The width of the slit.
    height : `~astropy.units.Quantity` ['angle'], default=3*u.arcmin
        The lenght of the slit.
    theta : `~astropy.units.Quantity` ['angle'], default=0*u.deg
        PA of the slit.
    resolution : int, optional
        The number of points that make up each side of the quadrangle -
        increase this to get a smoother quadrangle.
    vertex_unit : `~astropy.units.Unit` ['angle'], default=u.deg
        The units in which the resulting polygon should be defined - this
        should match the unit that the transformation (e.g. the WCS
        transformation) expects as input.

    Notes
    -----
    Additional keyword arguments are passed to `~matplotlib.patches.Polygon`
    """

    def __init__(self, anchor, width=1.5 * u.arcsec, height_up=3 * u.arcmin,
                 height_down=0 * u.arcmin,
                 theta=0 * u.deg, resolution=100, vertex_unit=u.deg, **kwargs):

        # Extract longitude/latitude, either from a tuple of two quantities, or
        # a single 2-element Quantity.
        lon_c, lat_c = u.Quantity(anchor).to_value(vertex_unit)
        center = np.array([[lon_c, lat_c]])

        theta = u.Quantity(theta).to_value(u.rad)

        # Convert the quadrangle dimensions to the appropriate units
        width = width.to_value(vertex_unit)
        height_up = height_up.to_value(vertex_unit)
        height_down = height_down.to_value(vertex_unit)

        # Corner coordinates
        longitude = lon_c - width * 0.5
        latitude = lat_c - height_down

        # Create progressions in longitude and latitude
        lon_seq = longitude + np.linspace(0, width, resolution + 1)
        lat_seq = latitude + np.linspace(0, height_down + height_up,
                                         resolution + 1)

        # Trace the path of the quadrangle
        lon = np.concatenate([lon_seq[:-1],
                              np.repeat(lon_seq[-1], resolution),
                              np.flip(lon_seq[1:]),
                              np.repeat(lon_seq[0], resolution)])
        lat = np.concatenate([np.repeat(lat_seq[0], resolution),
                              lat_seq[:-1],
                              np.repeat(lat_seq[-1], resolution),
                              np.flip(lat_seq[1:])])

        # Create polygon vertices
        vertices = np.array([lon, lat])

        # Rotation matrix
        rot_matrix = np.array([[np.cos(theta), +np.sin(theta)],
                               [-np.sin(theta), np.cos(theta)]])

        # Rotate Quadrangle
        vertices = vertices - center.T
        vertices = rot_matrix @ vertices
        vertices = vertices + center.T
        vertices = vertices.T

        super().__init__(vertices, **kwargs)