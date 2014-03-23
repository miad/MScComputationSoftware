#!/usr/bin/python
import numpy
import sys

hbar = 7.638233053909147;
mass = 723.4530251050578;

import sys

def main():
    if len(sys.argv) != 3:
        print "Usage: ./", sys.argv[0], " infile outfile"
        exit(1)
    outF = open(sys.argv[2], "w")
    for line in open(sys.argv[1], "r"):
        lt=line.strip().split('\t')
        nb=float(lt[0])+1j*float(lt[1])
        conv=numpy.sqrt(2*mass*nb)/hbar
        outF.write("%13.10e" % conv.real + " " + "%13.12e" % conv.imag + "\n");
    outF.close()

if __name__ == "__main__":
    main()
