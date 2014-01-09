\name{posteriorDecoding}
\alias{posteriorDecoding}
\title{
posteriorDecoding
}
\description{
Computes the posterior decoding table the posterior decoding (sequence of states) of the sequence in seqFilename. hmm should be a list on the form list("pi" = ..., "A" = ..., "B" = ...).
}
\usage{
pd = posteriorDecoding(seqFilename, hmm)
path = pd$path
table = pd$table
}
\value{
The table and the decoding are returned in a list: list("path" = ..., "table" = ...).
}
\references{
http://birc.au.dk/software/zipHMM
}
\author{
Andreas Sand <asand@birc.au.dk>
}
\seealso{ viterbi }
\keyword{ hidden Markov model }
\keyword{ posterior decoding }
