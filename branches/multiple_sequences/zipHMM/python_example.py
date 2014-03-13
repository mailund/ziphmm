from pyZipHMM import *

f1 = Forwarder.fromSequence(seqFilename = "../example.seq", alphabetSize = 3, minNoEvals = 10)
f1.writeToDirectory("../example_out")

## Do something crazy

f2 = Forwarder.fromDirectory(directory = "../example_out")
pi, A, B = readHMM("../example.hmm")

print "loglikelihood: ", f2.forward(pi, A, B)
# print "loglikelihood: ", f2.ptforward(pi, A, B) # parallelized version

pdPath, pdTable = posteriorDecoding("../example.seq", pi, A, B)
print "posterior path[0:10]:", pdPath[0:10]
print "posterior table[:0:7]:"
for r in xrange(pdTable.getHeight()):
    for c in xrange(7):
        print ("%.5f" % pdTable[r,c]), "\t",
    print
    

viterbiPath, viterbi_ll  = viterbi("../example.seq", pi, A, B)
print "viterbi log likelihood:", viterbi_ll
print "viterbi path[0:10]:", viterbiPath[0:10]


## pi, A and B can also be created using
# pi = Matrix(nStates, 0)
# pi[0,0] = 0.5
# ... and so on.
