\name{readHMMspec}
\alias{readHMMspec}
\title{
readHMMspec
}
\description{
Reads the specification (number of states and alphabet size) of an HMM from disk.
}
\usage{
spec = readHMMspec(filename)
noStates = spec$noStates
alphabetSize = spec$alphabetSize
}
\details{
For details on the file format see http://birc.au.dk/software/zipHMM.
}
\value{
The specifications are returned as a list: list("noStates" = ..., "alphabetSize" = ...).
}
\references{
http://birc.au.dk/software/zipHMM
}
\author{
Andreas Sand <asand@birc.au.dk>
}
\seealso{}
\keyword{ hidden Markov model }
