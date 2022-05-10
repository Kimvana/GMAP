# my lib imports
import GEMAIM


def help():
    main(["help"])


def main(callcommand):
    if len(callcommand) == 1:
        print(
            "\n\n\n"
            "Welcome to the GEMAIM package. Currently, this is a work in "
            "progress. The following options are available:\n\n\n"
            "    GMP\n"
            "prints this help.\n"
            "\n"
            "    GMP help\n"
            "prints this help.\n"
            "\n"
            "    GMP help [program]\n"
            "prints the help for the requested program\n"
            "\n"
            "    GMP AIM\n"
            "This program will create the required files for calculating "
            "infrared spectra.\n"
            "\n"
            "    GMP GEM\n"
            "This program will create the required files for calculating "
            "electronic spectra.\n"
        )
        return

    if hasattr(GEMAIM, callcommand[1]):
        progchoice = getattr(GEMAIM, callcommand[1])
        progchoice.help()

    else:
        print(
            "Program wasn't recognized. Please type the following "
            "for more information:\nGMP"
        )
