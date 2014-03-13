library('rZipHMM')

f = Forwarder$new()
f$readSeq(seqFilename = "example.seq", alphabetSize = 3, minNoEvals = 10)
f$writeToDirectory("example_out")

hmm = readHMM("example.hmm")

cat("loglikelihood:", f$forward(hmm), "\n")
# cat("loglikelihood:", f$ptforward(hmm), "\n") # parallelized version


pd = posteriorDecoding("example.seq", hmm)
cat("posterior path [1:10]:", pd$path[1:10], "\n")
cat("posterior table [1:10,:]:\n")
pd$table[,1:10]


v = viterbi("example.seq", hmm)
cat("Viterbi log-likelihood:", v$loglik, "\n")
cat("Viterbi path [1:10]:\n")
v$path[1:10]
