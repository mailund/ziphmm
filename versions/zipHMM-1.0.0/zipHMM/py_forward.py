from pysrf import *

import sys
import os
from subprocess import call


def printUsage(cmd):
    print "Usage: python ", cmd, "<hmm filename <seq filename>"

if __name__ == "__main__":
    if not len(sys.argv) == 3:
        printUsage(sys.argv[0])
        sys.exit(-1)

    hmmFilename = sys.argv[1]
    seqFilename = sys.argv[2]

    (pi, A, B) = readHMM(hmmFilename)
    nStates = pi.getHeight()
    nObservables = B.getWidth()

    f = Forwarder(seqFilename = seqFilename, nStates = nStates, nObservables = nObservables)
    loglik = f.forward(pi, A, B)
    print nStates, nObservables, loglik

