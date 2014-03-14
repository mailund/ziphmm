\name{rzipHMM-package}
\alias{rZipHMM-package}
\alias{rZipHMM}
\docType{package}
\title{
rZipHMM
}
\description{
zipHMM is a library for hidden Markov models that exploits repetitions
in strings to greatly speed up the calculations of the likelihood of a
sequence. The library is released under the LGPL license.

The library analyses the input string and finds repetitive patterns
and then reduces the string - similar to how compression algorithms
compress strings - by replacing substrings with new symbols. The new
symbols correspond to often seen substrings, and we can precompute the
probability for an HMM scanning over such strings. The full likelihood
can then be computed similar to the traditional forward algorithm, but
much faster since the algorithm can skip over often seen substrings. 

A typical use case of zipHMM is for training HMM parameters, where the
likelihood function is often computed for several hundred different
HMM parameters while the input sequences are kept stable. zipHMM
allows you to preprocess the input sequences once such that you
significantly reduce the time used on every likelihood
computation. The preprocessing can be saved to disk in a data
structure for later use.
}
\details{
\tabular{ll}{
Package: \tab rZipHMM\cr
Type: \tab Package\cr
Version: \tab 1.1.0\cr
Date: \tab 2014-01-09\cr
License: \tab LGPL\cr
}
readHMMspec(filename) reads the specification (number of states and alphabet size) of an HMM from disk. The specifications are returned as a list: list("noStates" = ..., "alphabetSize" = ...).
readHMM(filename) reads an entire HMM (number of states, alphabet size, initial state probabilities (pi), transition probabilities (A) and emission probabilities (B)) from disk. The HMM is returned as a list: list("pi" = ..., "A" = ..., "B" = ..., "noStates" = ..., "alphabetSize" = ...)
writeHMM(hmm, filename) writes an HMM to disk. hmm should be a list on the form list("pi" = ..., "A" = ..., "B" = ...)
viterbi(seqFilename, hmm) computes the Viterbi path and the log-likelihood of this path of the sequence in seqFilename. hmm should be a list on the form list("pi" = ..., "A" = ..., "B" = ...). The path and its log-loglikelihood are returned as a list: list("loglik" = ..., "path" = ...)
posteriorDecoding(seqFilename, hmm) computes the posterior decoding table the posterior decoding (sequence of states) of the sequence in seqFilename. hmm should be a list on the form list("pi" = ..., "A" = ..., "B" = ...). The table and the decoding are returned in a list: list("path" = ..., "table" = ...)
Forwarder$new() creates a new non-initialized Forwarder object.
Forwarder$readSeq(seqFilename, alphabetSize, nStatesSave = NULL, minNoEvals = 1) preprocesses the sequence in seqFilename and saves the corresponding data structure in the Forwarder object. alphabetSize should be the size of the alphabet used in the input sequence. nStatesSave should contain the model sizes that the data structure will likely used on in a vector. If it is empty or NULL, the data structure will be optimized for 2 states (longer preprocessing time). minNoEvals should be the expected number of time the data structure will be reused.
Forwarder$readSeqDirectory(dirname, alphabetSize, nStatesSave = NULL, minNoEvals = 1) preprocesses the sequences in the directory dirname and saves the corresponding data structure in the Forwarder object. alphabetSize should be the size of the alphabet used in the input sequences. nStatesSave should contain the model sizes that the data structure will likely be used on in a vector. If it is empty or NULL, the data structure will be optimized for 2 states (longer preprocessing time). minNoEvals should be the expected number of time the data structure will be reused.
Forwarder$readFromDirectory(directory, nStates = NULL) reads the saved data structure of an earlier preprocessing from disk. The entire data structure is read if nStates is NULL, while only the sub data structure corresponding to nStates is read if a modes size is given.
Forwarder$writeToDirectory(directory) writes the data structure of a preprocessed sequence to disk.
Forwarder$forward(hmm) computes the likelihood of a preprocessed sequence in a specific HMM. hmm should be a list on the form list("pi" = ..., "A" = ..., "B" = ..., ...)
Forwarder$ptforward(hmm, deviceFilename = NULL) is a parallelized version of Forwarder$forward. It computes the likelihood of a preprocessed sequence in a specific HMM. hmm should be a list on the form list("pi" = ..., "A" = ..., "B" = ..., ...). If several input sequences were used, the parallelization works across the lengths of the sequences. You want to use this method if you have few very long sequences.
Forwarder$ptforward(hmm, deviceFilename = NULL) is a parallelized version of Forwarder$forward. It computes the likelihood of a preprocessed sequence in a specific HMM. hmm should be a list on the form list("pi" = ..., "A" = ..., "B" = ..., ...). If several input sequences were used, the parallelization works across set of sequences. You want to use this function if you have a lot of sequences regardless of their lengths.
Forwarder$getOrigSeqLength() returns the length of the original sequence.
Forwarder$getOrigAlphabetSize() returns the alphabet size of the original sequence.
Forwarder$getSeqLength(noStates) returns the length of the preprocessed sequence corresponding to a model with noStates states.
Forwarder$getAlphabetSize(noStates) returns the alphabet size of the preprocessed sequence corresponding to a model with noStates states.
Forwarder$getPair(symbol) returns the two symbols composing 'symbol' in the extended alphabet.
}
\author{
Andreas Sand

Maintainer: Andreas Sand <asand@birc.au.dk>
}
\references{
http://birc.au.dk/software/zipHMM
}
\keyword{ package }
\keyword{ rZipHMM }
\keyword{ forward algorithm }
\keyword{ Viterbi algorithm }
\keyword{ posterior decoding }
\keyword{ optimization }
\keyword{ string compression }
