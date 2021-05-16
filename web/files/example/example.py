from pyZipHMM import *

f = Forwarder.fromSequence(seqFilename = "example.seq", alphabetSize = 3, minNoEvals = 10)
f.writeToDirectory(b"example_out")

pi, A, B = readHMM(b"example.hmm")

print("loglikelihood: ", f.forward(pi, A, B))


pdPath, pdTable = posteriorDecoding(b"example.seq", pi, A, B)
print("posterior path[0:10]:", pdPath[0:10])
print("posterior table[:0:10]:")
for r in range(pdTable.getHeight()):
    for c in range(10):
        print(("%.5f" % pdTable[r,c]), "\t", end=' ')
    print()
    

viterbiPath, viterbi_ll  = viterbi(b"example.seq", pi, A, B)
print("viterbi log likelihood:", viterbi_ll)
print("viterbi path[0:10]:", viterbiPath[0:10])
