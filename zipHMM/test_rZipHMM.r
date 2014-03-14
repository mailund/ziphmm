source('rZipHMM.r')

spec = readHMMspec("../test_data/test2.hmm")

stopifnot(spec$noStates == 3)
stopifnot(spec$alphabetSize == 4)

hmm = readHMM("../test_data/test2.hmm")
noStates = hmm$noStates
alphabetSize = hmm$alphabetSize

stopifnot(noStates == 3)
stopifnot(alphabetSize == 4)

writeHMM(hmm, "testtest.hmm")

v = viterbi("../test_data/test2.seq", hmm)
stopifnot(all.equal(v$loglik, -37.15332, tolerance=0.001))

p = posteriorDecoding("../test_data/test2.seq", hmm)

f = Forwarder$new()
f$readSeq(seqFilename = "../example.seq", alphabetSize = 3, minNoEvals = 10)

hmm = readHMM("../example.hmm")

ll = f$forward(hmm)
stopifnot(all.equal(ll, -10457.52, tolerance=0.001))

ptll = f$ptforward(hmm)
stopifnot(all.equal(ptll, -10457.52, tolerance=0.001))

mrll = f$mrforward(hmm)
stopifnot(all.equal(mrll, -10457.52, tolerance=0.001))

l = f$getOrigSeqLength()
stopifnot(l == 10000)

m = f$getOrigAlphabetSize()
stopifnot(m == 3)

lprime = f$getSeqLength(4)
stopifnot(lprime == 3426)

mprime = f$getAlphabetSize(4)
stopifnot(mprime == 34)

p = f$getPair(7)
stopifnot(p[1] == 1 && p[2] == 1)


hmm = readHMM("../test_data/multiple_seqs.hmm")

f = Forwarder$new()
f$readSeqDirectory(dirname = "../test_data/seqs_directory", alphabetSize = 3, minNoEvals = 10)

ll = f$forward(hmm)
stopifnot(all.equal(ll, -52287.60, tolerance=0.001))

ptll = f$ptforward(hmm)
stopifnot(all.equal(ptll, -52287.60, tolerance=0.001))

mrll = f$mrforward(hmm)
stopifnot(all.equal(mrll, -52287.60, tolerance=0.001))
