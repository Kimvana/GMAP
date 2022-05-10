# standard lib imports
import sys

# my lib imports
import GEMAIM


def main():
    callcommand = sys.argv
    if len(callcommand) == 1:
        callcommand.append("help")
    if len(callcommand) == 2:
        callcommand.append("help")

    choice = callcommand[1]
    if choice not in GEMAIM.alltools:
        print(
            "Choice of program wasn't recognized. Please type the following "
            "for more information on how to use this package:\nGMP"
        )
        return

    target_module = getattr(GEMAIM, choice)
    if callcommand[2] == "help":
        getattr(target_module, "help")()
    else:
        getattr(target_module, "main")(callcommand[1:])


if __name__ == "__main__":
    main()
