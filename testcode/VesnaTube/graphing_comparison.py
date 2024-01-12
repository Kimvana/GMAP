import matplotlib.pyplot as plt
import numpy as np

#some numbers are hardcoded because of predefined radii that wouldnt normally be
#instances of 100 would usually be "maxDist2"


def smoothnessPlotter(radii, smoothFactors, type):
    print(type)

    match type:
        case "linear":
            appliedFactors = list(map(lambda x : 1/(10*(1 - x)), smoothFactors))
            factorCurves = list(map(lambda x : weightFactorLinear(radii, x), appliedFactors))
            # print(len(factorCurves))
            # print(len(factorCurves[0]))
        case "quadratic":
            appliedFactors = list(map(lambda x : 1/(100 * (1 - x ** 2)), smoothFactors))
            factorCurves = list(map(lambda x : weightFactorQuadratic(radii, x ), appliedFactors))
        case _: 
            print("Smoothing type not available")
    return factorCurves
            


def weightFactorLinear(radii, appliedFactor):
    return list(map(lambda x : min((10 - x ) * appliedFactor, 1) , radii))

    # print(weightFactor)

    # plt.plot(radii, weightFactor)
    # plt.show()


def weightFactorQuadratic(radii, appliedFactor):
    return list(map(lambda x : min((100 - x ** 2) * appliedFactor, 1), radii))

def graphBuilder(radii, smoothFactors, types):
    print(types)
    for i in types:
        factorCurves = smoothnessPlotter(radii, smoothFactors, i)
        print(factorCurves)
        for j in factorCurves:
            plt.plot(radii,j)
    

def start():
    radii  = [r/100 for r in range(1, 1000)] #0.01 to 10
    smoothFactors = [0.1,0.5,0.9] #cant be 0
    types = ["linear", "quadratic"]

    graphBuilder(radii, smoothFactors, types)

    plt.show()

start()