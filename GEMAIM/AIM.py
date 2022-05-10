# standard lib imports
import sys


def help():
    print(
        "\n\n"
        "Usage:\n"
        "\n"
        "    GMP AIM demo\n"
        "Launches AIM in demo-mode. Performs a basic calculation to "
        "demonstrate basic use and to verify the program is installed "
        "correctly.\n"
        "\n"
        "    GMP AIM [name of input file]\n"
        "Performs a run of AIM using the parameters specified in the included "
        "file.\n"
        "\n\n"
        "The purpose of AIM is to take an MD trajectory and compute the "
        "time-dependent Hamiltonian to be used in infrared spectral "
        "calculations. Natively, only (protein) amide groups are supported, "
        "but AIM allows the user to specify their own oscillating groups.\n"
        "\n"
        "For more information, check the manual on github.com/Kimvana/AIM.\n"
    )


def main(callcommand):
    print("entered main of AIM - yet to be constructed")


if __name__ == "__main__":
    callcommand = sys.argv
    if len(callcommand) == 1:
        help()
    else:
        main(callcommand)
