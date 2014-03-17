#!/usr/bin/python
"""
This script will modify a config file for execution, and run it.
For usage, run the script with the -h flag.
"""
import sys, argparse, subprocess, os, errno, numpy
from multiprocessing import Pool
import itertools
import NumZero

NumberOfThreads=8

#pValues = numpy.arange(0.632, 0.644, 1.2E-4)
pValues = [0.6388322]
bValues = [18.91987]
#bValues = numpy.arange(18.86, 19.04, 1.8E-3)
cbValues=[1.00457, 1.00311, 0.99968, 0.98989]

allValues = list(itertools.product(pValues, bValues, cbValues))

rfoNames = ["KCurve", "KFound", "Potential", "PotentialPrecision", "InterestingPoints", "EnergyFile"]

def mkdir_p(path):
    """Create directory if does not exist."""
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise

def configFileNames(configDir):
    names = []
    for val in allValues:
        names.append(configDir + "/run_" + "%.5f" % val[0] + "_%.5f" % val[1] + "_%.5f" % val[2] + ".conf")
    return names

def outputDirNames(outputDir):
    names = []
    for val in allValues:
        names.append(outputDir + "/run_" + "%.5f" % val[0] + "_%.5f" % val[1] + "_%.5f" % val[2])
    return names
                     

                     
#Get content of config file
def generateConfigFiles(inFile, configDir, outputDir):
    oNames = outputDirNames(outputDir)
    cfNames = configFileNames(configDir)

    lines = []
    with open(inFile,"r") as f:
        for line in f:
            lines.append(line)

    for nameOfFile, val, dir0 in zip(cfNames, allValues, oNames):
        with open(nameOfFile,"w") as f:
            for l in lines:
                currLine = l
                for rfoName in rfoNames:
                    if rfoName + " = \"" in currLine:
                        currLine = rfoName + " = " + "\"" + str(dir0) + "/" + rfoName + ".dat" + "\"" + "\n"
                currLine = currLine.replace("_PV_", "%.10f" % val[0])
                currLine = currLine.replace("_BP_", "%.10fE-8" % val[1])
                currLine = currLine.replace("_CB_", "%.10f" % val[2])
                currLine = currLine.replace("_VRC_", "%.10f" % NumZero.GetZero(val[2], val[0], val[1]*1E-8))

                f.write(currLine )

    return cfNames, oNames

def computeRun(configFile):
    """Run the simulation on a specific macro file."""
    try:
    #sys.stdout.write("Running with config file " + configFile + "\n")
    #with open("/dev/null", "w") as f:
        with open(os.devnull, "w") as f:
            sys.stdout.write("Computing using file " + configFile + ".\n")
            subprocess.call(['./Compute',  '--configFile', configFile], stdout=f, stderr=f)
    except:
        print sys.stdout.write("Error: " + str(sys.exc_info()[0]));
        raise

def main():
    """Main method."""
    parser = argparse.ArgumentParser(description='Generate scriptfiles for the hbar_gshfs simulation.', epilog='For any questions, contact Rikard')
    parser.add_argument('infile',metavar='infile',type=str, help='Base (main) macro file to build new macro files.')
    parser.add_argument('configdir',metavar='configdir',type=str, help='Directory to store processed config files in.')
    parser.add_argument('outdir',metavar='outdir',type=str, help='Root file output directory.')
    parser.add_argument('-run', metavar='',type=bool, nargs='?', const=True, default=False, help="Specify this flag to actually launch the G4 processes in parallell.")
    parser.add_argument('-rm',metavar='',type=bool, nargs='?', const=True, default=False, help="Specify this flag to remove the macro files again after completion (usually used with the -run flag).")

    args = parser.parse_args()
    mkdir_p(args.configdir)
    sys.stdout.write("Generating " + str(len(allValues)) + " config files...\n")
    configFiles, outputDirectoriesToCreate = generateConfigFiles(args.infile, args.configdir, args.outdir)
    
    if args.run:
        mkdir_p(args.outdir)
        for dirac in outputDirectoriesToCreate:
            mkdir_p(dirac)

        pool=Pool(NumberOfThreads)
        for cFile in configFiles:
            pool.apply_async(computeRun,args=(cFile,))

        pool.close();
        pool.join();

    if args.rm:
        for f in configFiles:
            os.remove(f)
        #os.removedirs(args.configdir)

if __name__ == "__main__":
    """Launch."""
    main()

