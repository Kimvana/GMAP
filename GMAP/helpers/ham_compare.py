"""
The purpose of this script is to compare two hamiltonians to see whether
they're equal, and if not, indicate where / by how much.

Something similar is now implemented for the tests of the singles maps
AmideBB and AmideSC.
"""

import sys
from pathlib import Path

import numpy as np


def makeham(fname, n_singles):
    with open(fname, "r") as fhand:
        flatham = fhand.readline()
        flatham = flatham.split()
        # we have the final +1 as the frame number is still in front!
        if len(flatham) != (n_singles * (n_singles + 1)) / 2 + 1:
            print("warning! provided hamiltonian length doesn't match!")
            sys.exit()
        flatham = flatham[1:]
        squareham = np.zeros((n_singles, n_singles))
        squareham[np.triu_indices(n_singles, k=0)] = flatham
        squareham = squareham + squareham.T - np.diag(np.diag(squareham))
    return squareham


callargs = sys.argv

if len(callargs) < 2:
    fname1 = Path("../../hamiltonian.txt")
else:
    fname1 = Path(callargs[1])

if len(callargs) < 3:
    fname2 = Path(
        "D:/Data/PhD/AIM installable/AIM-version-1.0-installable/"
        "2024-10-18_16-50-03_AIM_V1-0-2_Hamiltonian.txt")
else:
    fname2 = Path(callargs[2])

ham1 = makeham(fname1, 128)
ham2 = makeham(fname2, 128)

ene1 = np.diag(ham1)
ene2 = np.diag(ham2)

if np.all(np.round(ene1, 3) == np.round(ene2, 3)):
    print("The two hamiltonians have the same energies")
else:
    print("There are the following differences between the two:")
    for ix, (en1, en2) in enumerate(zip(ene1, ene2)):
        if abs(en1 - en2) > 0.001:
            print(ix, en1, en2)
