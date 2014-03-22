#!/usr/bin/python
import scipy.optimize as so
import sys
import itertools
import numpy


muB = 6.717138840811653E5
xr =  9.974686590633269
Vzero = 3.326
xLimZero = [10, 100]
xLimBottom = [0, 9]
hbar = 7.638233053909147
mass = 723.4530251050578

KStdPoints = [2.5, 2.73, 2.85, 9.5]


def PotV(cB, p, Bp, RLOffset, x):    
    return (RLOffset + p*Vzero * (1.0-1.0/(1.0+pow(x/xr, 2) ) ) - muB * cB * Bp * x)


def GetZero(cB, p, Bp, RLOffset):
    return so.brentq(lambda x: PotV(cB, p, Bp, RLOffset, x), xLimZero[0], xLimZero[1])


def GetBottom(cB, p, Bp, RLOffset):
    return so.fminbound(lambda x: PotV(cB, p, Bp, RLOffset, x), xLimBottom[0], xLimBottom[1])

def GetKPoints(RLOffset):
    E = numpy.multiply(KStdPoints, KStdPoints) * pow(hbar, 2) / (2.0*mass)
    E += (RLOffset - 0.5)
    return numpy.sqrt(2*mass*E)/hbar

def main():
    cB = float(sys.argv[1])
    p = 0.6388322
    Bp = 18.91987E-8
    Offset = [0.5]
    for RLOffset in Offset:
        KP = GetKPoints(RLOffset)
        print RLOffset, GetZero(cB, p, Bp, RLOffset), GetBottom(cB, p, Bp, RLOffset), KP[0], KP[1], KP[2], KP[3]



if __name__ == "__main__":
    """Launch."""
    main()
