
# 3rd party imports
import numpy as np

# GMAP imports
import GMAP.src.tools.constants as GM_con
import GMAP.src.tools.math_functions as GM_mf
import GMAP.src.tools.print_tools as GM_pt


class NeighborMap:
    """Contains a single neighbormapping

    This mapping is a grid. When provided with two values (one for each
    direction), a bilinearly interpolated value is obtained.
    After making an instance of this class, the main use is for the
    instance to have the function get_delta called (which is in practice
    the application of the map).

    Parameters
    ----------
    fname : str or `pathlib.Path`
        the name of the file that contains the map information
    map_ : :class:`~GMAP.src.tools.MapReader.Map`
        The object that stores everything the program currently knows
        about this map.

    Attributes
    ----------
    dim : int
        The amount of rows and columns to expect from the file
    space : float or int
        The amount of degrees between two columns/rows. Should equal
        360/(dim-1)
    data : `np.ndarray`
        An array of shape (dim, dim) that contains the actual data
        points.
    """

    def __init__(self, fname, map_):
        self.dim = 13
        self.space = 30
        try:
            self.data = np.loadtxt(fname)
        except Exception as ex:
            GM_pt.Printer.warning(
                "\nCould not interpret the data "
                f" in the file {fname}. Please make sure "
                "the file was not changed since downloading, contains "
                "13 lines of 13 numbers, and no characters that cannot be "
                "interpreted as floats.",
                "map_AmideBB_3", exception=ex
            )
            map_.success = False
            return
        if self.data.shape != (13, 13):
            GM_pt.Printer.warning(
                "\nCould not interpret the data "
                f" in the file {fname}. Please make sure "
                "the file was not changed since downloading, and contains "
                "13 lines of 13 numbers.",
                "map_AmideBB_3"
            )
            map_.success = False

    def process_angle(self, angle):
        """Takes an angle in radians, and converts in to degrees, and
        also the column corresponding to the map.

        Parameters
        ----------
        angle : float
            The angle (in radians) to convert

        Returns
        -------
        angle : float
            The same angle, but in degrees
        angle_N : int
            The column that forms the 'left' boundary for interpolation.
            Found through int((angle + 180) // self.space)
        """

        angle *= GM_con.rad2deg
        angle_N = int((angle + 180) // self.space)
        if angle_N == (self.dim - 1):
            angle_N = self.dim - 2
        return angle, angle_N

    def get_delta(self, Nosc, Cosc, system):
        """Calculate the shift for this pair of oscillators.

        Parameters
        ----------
        Nosc, Cosc : :class:`~GMAP.src.tools.SystemReader.Oscillator`
            The oscillators for which the shift should be calculated.
        system : :class:`~GMAP.src.tools.SystemReader.System`
            The object that stores everything the program currently
            knows about the system being treated (names, numbers, types,
            masses, charges of all atoms, for example)

        Returns
        -------
        delta : float
            The amount that should be added to the frequency to
            compensate for interactions between these two oscillators.
        """

        # calculate the ramachandran angles
        Nbpos = Nosc.positions_box
        Cbpos = Cosc.positions_box
        phi_ang = GM_mf.dihedral_boxcoords(  # around NCA, angle of two Cs
            Nbpos[0], Nbpos[3], Nbpos[5], Cbpos[0], system.boxvects)
        phi_ang, phi_N = self.process_angle(phi_ang)
        psi_ang = GM_mf.dihedral_boxcoords(  # around CAC, angle of two Ns
            Nbpos[3], Nbpos[5], Cbpos[0], Cbpos[3], system.boxvects)
        psi_ang, psi_N = self.process_angle(psi_ang)

        if (0 <= phi_N < (self.dim - 1)) and (0 <= psi_N < (self.dim - 1)):
            y1 = self.data[psi_N, phi_N]
            y2 = self.data[psi_N+1, phi_N]
            y3 = self.data[psi_N+1, phi_N+1]
            y4 = self.data[psi_N, phi_N+1]

            u = (phi_ang % self.space) / self.space
            t = (psi_ang % self.space) / self.space

            # bilinear interpolation!
            delta = (1-u)*(1-t)*y1 + (1-u)*t*y2 + u*t*y3 + u*(1-t)*y4

        else:
            GM_pt.Printer.warning(
                "Ill defined ramachandran angles found for residue "
                f"{Nosc.resnames[1]}{Nosc.resnums[1]}. The nearest "
                "neighbour shift will be set to zero.", "map_AmideBB_5"
            )
            delta = 0
        return delta


def read_maps(map_):
    """Find and read all files with neighbor maps.
    Creates a NeighborMap instance for each file, and adds them to the
    map for later use.

    Parameters
    ----------
    map_ : :class:`~GMAP.src.tools.MapReader.Map`
        The object that stores everything the program currently knows
        about this map.
    """

    mapdir = map_.directory / "neighbor_maps"
    map_.neighbormaps = {
        fname.stem: NeighborMap(fname, map_) for fname in mapdir.iterdir()}
