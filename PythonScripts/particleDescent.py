#!/usr/bin/python
"""
This script will modify a config file for execution, and run it.
For usage, run the script with the -h flag.
"""
import sys, argparse, subprocess, os, errno, numpy
from multiprocessing import Pool
import itertools
import tempfile

cbValues=[1.00457, 1.00311, 0.99968, 0.98989]

expectedRates = [35.25, 30.12, 21.76, 8.28]
expectedSigma = [3.57, 2.81, 1.12, 0.49]

dPV = 1E-5
dBV = 1E-5

hbar = 7.638233053909147
mass = 723.4530251050578

inFileContents = []

tempCounter = 0

def mkdir_p(path):
    """Create directory if does not exist."""
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise
                     

def generateConfig(CBval, Bval, Pval):
    configFD, configFilename = tempfile.mkstemp()
    os.close(configFD)
    outputFD, outputFilename = tempfile.mkstemp()
    os.close(outputFD)

    with open(configFilename, "w") as f:
        for l in inFileContents:
            currLine = l
            currLine = currLine.replace("_PV_", "%.10f" % Pval)
            currLine = currLine.replace("_BP_", "%.10fE-8" % Bval)
            currLine = currLine.replace("_CB_", "%.10f" % CBval)
            currLine = currLine.replace("_INTERESTING_POINTS_FILE_", outputFilename)
            f.write(currLine)

    return configFilename, outputFilename

def computeRun(configFile):
    """Run the simulation on a specific macro file."""
    try:
    #sys.stdout.write("Running with config file " + configFile + "\n")
    #with open("/dev/null", "w") as f:
        with open(os.devnull, "w") as f:
            subprocess.call(['./Compute',  '--configFile', configFile], stdout=f, stderr=f)
    except:
        print sys.stdout.write("Error: " + str(sys.exc_info()[0]));
        raise

def K2E(k):
    return k*k*hbar*hbar/(2*mass)

def cleanup(configFile, outputFile):
    os.remove(configFile)
    os.remove(outputFile)
    

def ComputeEnergy(CBval, Bval, Pval):
    configFile, outputFile = generateConfig(CBval, Bval, Pval)
    computeRun(configFile)
    E = []
    with open(outputFile, "r") as f:
        for line in f:
            re, im = [float(x) for x in line.split()]
            E.append(K2E(re + 1j*im))
    cleanup(configFile, outputFile)
    if len(E) != 1:
        print "Warning: found " + str(len(E)) + " energies for cB = ", CBval, "B = ", Bval, " P = ", Pval
        return 1E99

    return -2*E[0].imag/hbar*1.0E6

def ComputeEnergyMapper(myTriplet):
    return ComputeEnergy(myTriplet[0], myTriplet[1], myTriplet[2])
                      

def Func(Bval, Pval):
    pool = Pool(4)
    rates = pool.map(ComputeEnergyMapper, (zip(cbValues, [Bval] * 4, [Pval] * 4) ) )
    pool.close()
    pool.join()
    ans = sum(pow((x-y)/z,2) for (x, y, z) in zip(rates, expectedRates, expectedSigma))
    return ans

def Minimize(direction, params):
    Bval = params[0]
    Pval = params[1]
    Bdir = direction[0]
    Pdir = direction[1]
    x = 0.0
    xold = 1.0
    dx = numpy.sqrt(dBV*dBV + dPV * dPV)*1E-4

    #Do Newton-Rhapson
    while abs(x-xold) > 1E-9:
        xold = x
        f = Func(Bval + Bdir * x, Pval + Pdir * x)
        fPlus = Func(Bval + Bdir*(x + dx), Pval + Pdir * (x + dx)) 
        fMinus =  Func(Bval + Bdir*(x - dx), Pval + Pdir * (x - dx))

        fPrime = (fPlus - fMinus) / (2 * dx)
        fBis = (fPlus + fMinus - 2 * f)/(dx*dx)

        if fBis != 0:
            x = x - fPrime / fBis
        else:
            break
    return [Bval + Bdir * x, Pval + Pdir * x]

    
    
    

def PerformStep(params):
    Bval = params[0]
    Pval = params[1]
    f = Func(Bval, Pval)
    dfdB = ( Func(Bval + dBV, Pval) - Func(Bval - dBV, Pval) ) / (2 * dBV)
    dfdP = ( Func(Bval, Pval + dPV) - Func(Bval, Pval - dPV) ) / (2 * dPV)
    return  Minimize([-dfdB, -dfdP], params), f


def GetParser():
    parser = argparse.ArgumentParser(description='Generate scriptfiles for the hbar_gshfs simulation.', epilog='For any questions, contact Rikard')
    parser.add_argument('infile',metavar='infile',type=str, help='Base (main) macro file to build new macro files.')
    parser.add_argument('resultfile',metavar='resultfile',type=str, help='File to dump result to')
    parser.add_argument('startB', metavar='startB', type=float, help='Initial B-value')
    parser.add_argument('startP', metavar='startP', type=float, help='Initial P-value')
    return parser

def ReadInfile(inFile):
    with open(inFile,"r") as f:
        for line in f:
            inFileContents.append(line)

def main():
    parser = GetParser()
    args = parser.parse_args()
    ReadInfile(args.infile)

    #startValues = [18.98, 0.64]
    currentValues = [args.startB, args.startP]
    currentResult = 1E99
    nextResult = 1E10
    
    #May there be stepping here.
    with open(args.resultfile,"w", 0) as f:
        while abs(currentResult - nextResult) > 1E-9:
            currentResult = nextResult
            f.write("%.10f" % currentValues[0] + " %.10f" % currentValues[1] + " %.10f" % currentResult + "\n")
            f.flush() #Since we'll probably abort the process at some point anyway.
            print currentValues, " ", currentResult
            currentValues, nextResult = PerformStep(currentValues)
        

if __name__ == "__main__":
    """Launch."""
    main()
