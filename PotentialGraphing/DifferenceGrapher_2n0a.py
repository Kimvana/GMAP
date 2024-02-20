import numpy as np
import matplotlib.pyplot as plt
import MDAnalysis as MDA


def read_files(folder, spheresize):
    """
    Turns a .txt file of potentials into a list of floats of those
    potentials. The .txt file must be of the format: 1 pot1 pot2 ...
    potn \n

    :param folder: Name of the folder in which the .txt file resides
        and the prefix of the .txt file.
    :type folder: str
    :param spheresize: Size of the spheresize used to compute the
        potentials in the .txt file.
    :type spheresize: int
    :return: List of the potentials, ordered in ascending atom number.
    :rtype: list[float]
    """

    spheresize_pot_list = open(f"{folder}{spheresize}.txt").readline()
    pot_list = spheresize_pot_list.split(" ")[1:-1]
    return [float(pot) for pot in pot_list]


def system_arranger(system, potentials):
    """
    This function takes a system made by MDA and the potentials and
    combines them into one list.

    :param system: Contains a list of tuples, where each tuple
        represents an atom. The first item in each tuple is the name
        of the atom, which is not per se its type. The second item
        in each tuple is the number of the residue in which the atom
        resides.
    :type system: list[tuple[str,str]]
    :param potentials: Contains a list of floats, where each float
        represents the potential at the point of an atom compared
        to a position infinitely far away from the atom.
    :type potentials: list[float]
    :return: Contains a list of lists where each sublist represents
        a residue. The sublists contain tuples which represent the
        atoms in the residue. The first item in each tuple is the
        name of the atom and the second item is the potential at
        the position of atom compared to a position infinitely far
        away from the atom.
    :rtype: list[list[tuple[str, float]]]
    """

    arranged_system = []
    cur_res = 0  # The residue currently being processed
    for res, pot in zip(system, potentials):
        # Start working on a new residue if the atom from MDA no
        # longer belongs to the current residue
        if cur_res != res[1]:
            arranged_system.append([])
            cur_res += 1
        # Add the atom to the current residue
        arranged_system[-1].append([res[0], pot])
    return arranged_system


def potential_difference(arranged_system, relevant_atoms):
    """
    Calculates the absolute potential difference between relevant
    atoms if they are in the same residue and returns the average.

    :param arranged_system: Contains a list of lists where each
        sublist represents a residue. The sublists contain tuples
        which represent the atoms in the residue. The first item
        in each tuple is the name of the atom and the second item
        is the potential at the position of atom compared to a
        position infinitely far away from the atom.
    :type arranged_system: list[list[tuple[str, float]]]
    :param relevant_atoms: A list that contains all the atom names.
        The potentials between atoms with any of these names will
        be calculated if they are in the same residue.
    :type relevant_atoms: list[str]
    :return: The average absolute potential between all relevant atoms.
    :rtype: int
    """

    total_potential_difference = 0
    total_number_of_potentials = 0
    for residue in arranged_system:
        potentials = []  # All of the potentials between atoms
        for atom in residue:
            if atom[0] in relevant_atoms and not np.isnan(atom[1]):
                potentials.append(atom[1])
        while potentials:
            for potential in potentials[1:]:
                total_potential_difference += abs(potentials[0] - potential)
                total_number_of_potentials += 1
            potentials = potentials[1:]
    average_potential = total_potential_difference / total_number_of_potentials
    return average_potential


def potential_difference_grapher(topfile, trjfile, potential_folders,
                                 radii, relevant_atoms):
    """
    Creates .png files displaying the potentials versus spheresize.
    The potentials are the averages of the absolute potential
    between relevant atoms that are in the same residue.

    :param topfile: Name of the topology file and its path.
    :type topfile: str
    :param trjfile: Name of the trajectory file and its path.
    :type trjfile: str
    :param potential_folders: List containing the names of the folders
        where the potentials are stored. The folders are where the
        .txt files and the prefix of the .txt files reside.
    :type potential_folders: list[str]
    :param radii: List containing the spheresizes. Each spheresize
        corresponds to a .txt file.
    :type radii: list[int]
    :param relevant_atoms: A list that contains all the atom names.
        The potentials between atoms with any of these names will
        be calculated if they are in the same residue.
    :type relevant_atoms: list[str]
    :return: None
    """
    universe = MDA.Universe(topfile, trjfile)
    system = [*zip(universe.atoms.names, universe.atoms.resnums)]
    for folder in potential_folders:
        pot_out = []
        for spheresize in radii:
            potentials = read_files(folder, spheresize)
            arranged_system = system_arranger(system, potentials)
            pot_out.append(potential_difference(arranged_system,
                                                relevant_atoms))
        plt.plot(radii, pot_out)
        plt.xlabel("Spheresize (A)")
        plt.ylabel("Potential (V)")
    plt.legend(["ma", "mm"])
    plt.savefig("2N0A Graph")
    plt.clf


if __name__ == "__main__":

    topfile = "2n0a.tpr"
    trjfile = "2n0a.xtc"

    potential_folders = ["2n0a_ma/subbox_ma_", "2n0a_perres/perres_mm_"]
    radii = np.arange(10, 59)  # Used spheresizes
    relevant_atoms = ["O", "N"]

    potential_difference_grapher(topfile, trjfile, potential_folders,
                                 radii, relevant_atoms)
