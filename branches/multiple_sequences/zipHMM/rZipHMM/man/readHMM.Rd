\name{readHMM}
\alias{readHMM}
\title{
readHMM
}
\description{
Reads an entire HMM (number of states, alphabet size, initial state probabilities (pi), transition probabilities (A) and emission probabilities (B)) from disk. 
}
\usage{
hmm = readHMM(filename)
pi = hmm$pi
A = hmm$A
B = hmm$B
noStates = hmm$noStates
alphabetSize = hmm$alphabetSize
}
\details{
For details on the file format see http://birc.au.dk/software/zipHMM.
}
\value{
The HMM is returned as a list: list("pi" = ..., "A" = ..., "B" = ..., "noStates" = ..., "alphabetSize" = ...).
}
\references{
http://birc.au.dk/software/zipHMM
}
\author{
Andreas Sand <asand@birc.au.dk>
}
\seealso{}
\keyword{ hidden Markov model }
