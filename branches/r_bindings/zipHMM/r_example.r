source('rZipHMM.r')

f1 = Forwarder$new()
f1$readSeq(seqFilename = "../example.seq", alphabetSize = 3, minNoEvals = 10)

f1$writeToDirectory("../example_out")

## Do something crazy

f2 = Forwarder$new()
f2$readFromDirectory(directory = "../example_out")
hmm = readHMM("../example.hmm")

ll = f2$forward(hmm)
#ll = f2$ptforward(hmm) # parallelized version
cat("loglikelihood:", ll, "\n")


pd = posteriorDecoding("../example.seq", hmm)
cat("posterior path [1:10]:", pd$path[1:10], "\n")
cat("posterior table [1:8,:]:\n")
pd$table[,1:8]


v = viterbi("../example.seq", hmm)
cat("Viterbi log-likelihood:", v$loglik, "\n")
cat("Viterbi path [1:10]:\n")
v$path[1:10]
