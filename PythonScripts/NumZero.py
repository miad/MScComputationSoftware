#!/usr/bin/python
import scipy.optimize as so

RLOffset = 0.5
muB = 6.717138840811653E5
xr =  9.974686590633269
Vzero = 3.326
xLim = [10, 20]



def PotV(cB, p, Bp, x):    
    return (RLOffset + p*Vzero * (1.0-1.0/(1.0+pow(x/xr, 2) ) ) - muB * cB * Bp * x)


def GetZero(cB, p, Bp):
    return so.brentq(lambda x: PotV(cB, p, Bp, x), xLim[0], xLim[1])



#def main():
    #print GetZero(1.00457, 0.63811, 18.92E-8)


#if __name__ == "__main__":
#    """Launch."""
#    main()
