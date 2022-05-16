# standard lib imports
import sys


def help():
    print(
        "\n\n"
        "Usage:\n"
        "\n"
        "    GMP GEM demo\n"
        "Launches GEM in demo-mode. Performs a basic calculation to "
        "demonstrate basic use and to verify the program is installed "
        "correctly.\n"
        "\n"
        "    GMP GEM [name of input file]\n"
        "Performs a run of GEM using the parameters specified in the included "
        "file.\n"
        "\n\n"
        "The purpose of GEM is to take an MD trajectory and compute the "
        "time-dependent Hamiltonian to be used in electronic spectral "
        "calculations. Instructions on how to deal with specific "
        "chromophores have to be included in the corresponding .emap file.\n"
        "\n"
        "For more information, check the manual on N/A.\n"
    )


def main(callcommand):
    print("entered main of GEM - yet to be constructed")


if __name__ == "__main__":
    callcommand = sys.argv
    if len(callcommand) == 1:
        help()
    else:
        main(callcommand)
