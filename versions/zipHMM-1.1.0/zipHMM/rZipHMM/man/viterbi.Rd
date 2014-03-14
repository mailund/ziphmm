\name{viterbi}
\alias{viterbi}
\title{
viterbi
}
\description{
Computes the Viterbi path and the log-likelihood of this path of the sequence in seqFilename. hmm should be a list on the form list("pi" = ..., "A" = ..., "B" = ...).
}
\usage{
v = viterbi(seqFilename, hmm)
path = v$path
ll = v$loglik
}
\value{
The Viterbi path and its log-loglikelihood are returned as a list: list("loglik" = ..., "path" = ...)
}
\references{
http://birc.au.dk/software/zipHMM
}
\author{
Andreas Sand <asand@birc.au.dk>
}
\seealso{}
\keyword{ hidden Markov model }
\keyword{ the Viterbi algorithm }
