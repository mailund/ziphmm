from pyZipHMM import *

f = Forwarder.fromSequenceDirectory(dirname = "sequences", alphabetSize = 3, minNoEvals = 10)
f.writeToDirectory("example_out")

pi, A, B = readHMM("example.hmm")

print "loglikelihood: ", f.forward(pi, A, B)
# print "loglikelihood: ", f.ptforward(pi, A, B)
# print "loglikelihood: ", f.mrforward(pi, A, B)
