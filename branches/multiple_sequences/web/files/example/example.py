from pyZipHMM import *

f = Forwarder.fromSequence(seqFilename = "example.seq", alphabetSize = 3, minNoEvals = 10)
f.writeToDirectory("example_out")

pi, A, B = readHMM("example.hmm")

print "loglikelihood: ", f.forward(pi, A, B)


pdPath, pdTable = posteriorDecoding("example.seq", pi, A, B)
print "posterior path[0:10]:", pdPath[0:10]
print "posterior table[:0:10]:"
for r in xrange(pdTable.getHeight()):
    for c in xrange(10):
        print ("%.5f" % pdTable[r,c]), "\t",
    print
    

viterbiPath, viterbi_ll  = viterbi("example.seq", pi, A, B)
print "viterbi log likelihood:", viterbi_ll
print "viterbi path[0:10]:", viterbiPath[0:10]
