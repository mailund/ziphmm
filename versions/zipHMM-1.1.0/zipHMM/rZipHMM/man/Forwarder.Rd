\name{Forwarder-class}
\docType{class}
\alias{Forwarder-class}
\alias{Forwarder}

\title{Class \code{"Forwarder"}}
\description{
The Forwarder class is used for preprocessing an input sequence and for computing the log-likelihood of this sequence very fast given a specific HMM. 

To preprocess a sequence, the method readSeq is used as illustrated below.

The data structure resulting from the preprocessing can be saved to disk, using the method writeToDirectory.

A saved data structure can be read from disk using the method readFromDirectory.

Finally, the forward method (or the parallelized ptforward) is used to compute the log-likelihood of the sequence given a specific HMM.
}
\section{Extends}{

All reference classes extend and inherit methods from \code{"\linkS4class{envRefClass}"}.

}
\references{
http://birc.au.dk/software/zipHMM
}
\author{
Andreas Sand <asand@birc.au.dk>
}

\seealso{
readHMMspec
readHMM
writeHMM
}
\usage{
f1 = Forwarder$new()
f1$readSeqDirectory(dirname = "sequences/", alphabetSize = 3, minNoEvals = 10)

f1$writeToDirectory("example_out")

## Do something crazy

f2 = Forwarder$new()
f2$readFromDirectory(directory = "example_out")
hmm = readHMM("example.hmm")

ll = f2$forward(hmm) # sequential implementation
#ll = f2$ptforward(hmm) # parallelization across the length of the sequences
#ll = f2$mrforward(hmm) # parallelization across the set of sequences
cat("loglikelihood:", ll, "\n")

l = f2$getOrigSeqLength()
m = f2$getOrigAlphabetSize()

lprime = f2$getSeqLength(4)
mprime = f2$getAlphabetSize(4)

p = f$getPair(7)

}
\keyword{ Forwarder }
\keyword{ the forward algorithm }

\section{Fields}{
  \describe{
    \item{\code{ptr}:}{Address of the C++ object that is wrapped by this object.}
  }
}
\section{Methods}{
  \describe{
    \item{\code{getPair(symbol)}:}{ returns the two symbols composing 'symbol' in the extended alphabet. }
    \item{\code{getAlphabetSize(noStates)}:}{ returns the alphabet size of the preprocessed sequence corresponding to a model with noStates states. }
    \item{\code{getSeqLength(noStates)}:}{ returns the length of the preprocessed sequence corresponding to a model with noStates states. }
    \item{\code{getOrigAlphabetSize()}:}{ returns the alphabet size of the original sequence. }
    \item{\code{getOrigSeqLength()}:}{ returns the length of the original sequence. }
    \item{\code{ptforward(hmm, deviceFilename)}:}{ is a parallelized version of Forwarder$forward. It computes the likelihood of a preprocessed sequence in a specific HMM. hmm should be a list on the form list("pi" = ..., "A" = ..., "B" = ..., ...). }
    \item{\code{forward(hmm)}:}{ computes the likelihood of a preprocessed sequence in a specific HMM. hmm should be a list on the form list("pi" = ..., "A" = ..., "B" = ..., ...) }
    \item{\code{writeToDirectory(directory)}:}{ writes the data structure of a preprocessed sequence to disk. }
    \item{\code{readFromDirectory(directory, nStates)}:}{ reads the saved data structure of an earlier preprocessing from disk. The entire data structure is read if nStates is NULL, while only the sub data structure corresponding to nStates is read if a modes size is given. }
    \item{\code{readSeq(seqFilename, alphabetSize, nStatesSave, minNoEvals)}:}{ preprocesses the sequence in seqFilename and saves the corresponding data structure in the Forwarder object. alphabetSize should be the size of the alphabet used in the input sequence. nStatesSave should contain the model sizes that the data structure will likely used on in a vector. If it is empty or NULL, the data structure will be optimized for 2 states (longer preprocessing time). minNoEvals should be the expected number of time the data structure will be reused. }
    \item{\code{initialize()}:}{ Initialize object. }
  }
}
