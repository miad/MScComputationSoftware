#!/usr/bin/python
"""
This script will modify a config file for execution, and run it.
For usage, run the script with the -h flag.
"""
import sys, argparse, subprocess, os, errno, numpy
from multiprocessing import Pool


potDepths = numpy.arange(-5.0, 5.0, 0.05)
#for n in potDepths:
#    n = round(n, 3)

rfoNames = ["KCurve", "KFound", "Potential", "PotentialPrecision", "InterestingPoints", "Wavefunctions", "EnergyFile", "ProductWavefunction"]

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
    for d in potDepths:
        names.append(configDir + "/run_" + "%.2f" % d + ".conf")
    return names

def outputDirNames(outputDir):
    names = []
    for d in potDepths:
        names.append(outputDir + "/output_" + "%.2f" % d )
    return names
                     

                     
#Get content of config file
def generateConfigFiles(inFile, configDir, outputDir):
    oNames = outputDirNames(outputDir)
    cfNames = configFileNames(configDir)

    lines = []
    with open(inFile,"r") as f:
        for line in f:
            lines.append(line)

    for nameOfFile, depth, dir0 in zip(cfNames, potDepths, oNames):
        with open(nameOfFile,"w") as f:
            for l in lines:
                currLine = l
                for rfoName in rfoNames:
                    if rfoName + " = \"" in currLine:
                        currLine = rfoName + " = " + "\"" + str(dir0) + "/" + rfoName + ".dat" + "\"" + "\n"
                minusOneOverg = -1.0 / depth
                currLine = currLine.replace("DEPTH", "%.10f" % minusOneOverg)
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
    sys.stdout.write("Generating config files...\n")
    configFiles, outputDirectoriesToCreate = generateConfigFiles(args.infile, args.configdir, args.outdir)
    
    if args.run:
        mkdir_p(args.outdir)
        for dirac in outputDirectoriesToCreate:
            mkdir_p(dirac)

        pool=Pool(12)
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

