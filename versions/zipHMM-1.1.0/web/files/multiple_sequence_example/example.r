library('rZipHMM')

f = Forwarder$new()
f$readSeqDirectory(dirname = "sequences", alphabetSize = 3, minNoEvals = 10)
f$writeToDirectory("example_out")

hmm = readHMM("example.hmm")

cat("loglikelihood:", f$forward(hmm), "\n")
# cat("loglikelihood:", f$ptforward(hmm), "\n") # parallelized version
# cat("loglikelihood:", f$mrforward(hmm), "\n") # parallelized version
