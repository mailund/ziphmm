dyn.load("librZipHMM.so")

readHMMspec <- function(filename) {
  if(missing(filename))
    stop("filename is missing")
  noStates = as.integer(0)
  alphabetSize = as.integer(0)
  .Call("c_read_HMM_spec", noStates, alphabetSize, filename);
  return(list("noStates" = noStates, "alphabetSize" = alphabetSize))
}

readHMM <- function(filename) {
  if(missing(filename))
    stop("filename is missing")
  spec = readHMMspec(filename)
  noStates = spec$noStates
  alphabetSize = spec$alphabetSize

  pi = vector(mode = 'numeric', length = noStates)
  A  = matrix(0, nrow = noStates, ncol = noStates)
  B  = matrix(0, nrow = noStates, ncol = alphabetSize)

  .Call("c_read_HMM", pi, A, B, filename)
  
  return(list("pi" = pi, "A" = A, "B" = B, "noStates" = noStates, "alphabetSize" = alphabetSize))
}

writeHMM <- function(hmm, filename) {
  if(missing(hmm))
    stop("hmm is missing")
  if(missing(filename))
    stop("filename is missing")
  
  if(!is.list(hmm) || !("pi" %in% names(hmm) && "A" %in% names(hmm) && "B" %in% names(hmm)))
    stop("hmm should be a list of the form list(\"pi\" = ..., \"A\" = ..., \"B\" = ...)")
  
  pi = hmm$pi
  A = hmm$A
  B = hmm$B

  if(!is.vector(pi))
    stop("hmm$pi should be a vector.")
  if(!is.matrix(A))
    stop("hmm$A should be a matrix.")
  if(!is.matrix(B))
    stop("hmm$A should be a matrix.")
  
  noStates = dim(B)[1]
  alphabetSize = dim(B)[2]
  
  .Call("c_write_HMM", pi, A, B, noStates, alphabetSize, filename)

  cat()
}


Forwarder <- setRefClass(
  "Forwarder",

  fields = list(
    ptr = "ANY"
  ),

  methods = list(

    initialize = function() {
      ptr <<- .Call("Forwarder_new");
    },

    readSeq = function(seqFilename, alphabetSize, nStatesSave = NULL, minNoEvals = 1) {
      if(missing(seqFilename))
        stop("seqFilename is missing")
      if(missing(alphabetSize))
        stop("alphabetSize is missing")
      
      if(!is.numeric(alphabetSize))
        stop("alphabetSize should be a number")
      alphabetSize = as.integer(alphabetSize)

      if(missing(minNoEvals))
        minNoEvals = 1
      else if(!is.numeric(minNoEvals))
        stop("minNoEvals should be a number")
      minNoEvals = as.integer(minNoEvals)      

      if(!missing(nStatesSave))
        .Call("Forwarder_read_seq", ptr, seqFilename, alphabetSize, nStatesSave, minNoEvals)
      else
        .Call("Forwarder_read_seq", ptr, seqFilename, alphabetSize, NULL, minNoEvals)

      cat()
    },

    readSeqDirectory = function(dirname, alphabetSize, nStatesSave = NULL, minNoEvals = 1) {
      if(missing(dirname))
        stop("dirname is missing")
      if(missing(alphabetSize))
        stop("alphabetSize is missing")
      
      if(!is.numeric(alphabetSize))
        stop("alphabetSize should be a number")
      alphabetSize = as.integer(alphabetSize)

      if(missing(minNoEvals))
        minNoEvals = 1
      else if(!is.numeric(minNoEvals))
        stop("minNoEvals should be a number")
      minNoEvals = as.integer(minNoEvals)      

      if(!missing(nStatesSave))
        .Call("Forwarder_read_seq_directory", ptr, dirname, alphabetSize, nStatesSave, minNoEvals)
      else
        .Call("Forwarder_read_seq_directory", ptr, dirname, alphabetSize, NULL, minNoEvals)

      cat()
    },

    readFromDirectory = function(directory, nStates = NULL) {
      if(missing(directory))
        stop("directory is missing")
      
      if(!missing(nStates)) {
        .Call("Forwarder_read_from_directory_and_no_states", ptr, directory, nStates)
      } else {
        .Call("Forwarder_read_from_directory", ptr, directory)
      }

      cat()
    },

    writeToDirectory = function(directory) {
      .Call("Forwarder_write_to_directory", ptr, directory)

      cat()
    },

    forward = function(hmm) {
      if(missing(hmm))
        stop("hmm is missing")

      if(!is.list(hmm) || !("pi" %in% names(hmm) && "A" %in% names(hmm) && "B" %in% names(hmm)))
        stop("hmm should be a list of the form list(\"pi\" = ..., \"A\" = ..., \"B\" = ...)")
      
      pi = hmm$pi
      A = hmm$A
      B = hmm$B
      
      if(!is.vector(pi))
        stop("hmm$pi should be a vector.")
      if(!is.matrix(A))
        stop("hmm$A should be a matrix.")
      if(!is.matrix(B))
        stop("hmm$A should be a matrix.")
      
      noStates = dim(B)[1]
      alphabetSize = dim(B)[2]
      
      ll = .Call("Forwarder_forward", ptr, pi, A, B, noStates, alphabetSize)

      return(ll)
    },

    ptforward = function(hmm, deviceFilename = NULL) {
      if(missing(hmm))
        stop("hmm is missing")

      if(!is.list(hmm) || !("pi" %in% names(hmm) && "A" %in% names(hmm) && "B" %in% names(hmm)))
        stop("hmm should be a list of the form list(\"pi\" = ..., \"A\" = ..., \"B\" = ...)")
      
      pi = hmm$pi
      A = hmm$A
      B = hmm$B
      
      if(!is.vector(pi))
        stop("hmm$pi should be a vector.")
      if(!is.matrix(A))
        stop("hmm$A should be a matrix.")
      if(!is.matrix(B))
        stop("hmm$A should be a matrix.")
            
      noStates = dim(B)[1]
      alphabetSize = dim(B)[2]

      if(missing(deviceFilename)) {
        .Call("Forwarder_ptforward", ptr, pi, A, B, noStates, alphabetSize, "-")
      } else {
        .Call("Forwarder_ptforward", ptr, pi, A, B, noStates, alphabetSize, deviceFilename)
      }
      
    },

    mrforward = function(hmm, deviceFilename = NULL) {
      if(missing(hmm))
        stop("hmm is missing")

      if(!is.list(hmm) || !("pi" %in% names(hmm) && "A" %in% names(hmm) && "B" %in% names(hmm)))
        stop("hmm should be a list of the form list(\"pi\" = ..., \"A\" = ..., \"B\" = ...)")
      
      pi = hmm$pi
      A = hmm$A
      B = hmm$B
      
      if(!is.vector(pi))
        stop("hmm$pi should be a vector.")
      if(!is.matrix(A))
        stop("hmm$A should be a matrix.")
      if(!is.matrix(B))
        stop("hmm$A should be a matrix.")
            
      noStates = dim(B)[1]
      alphabetSize = dim(B)[2]

      if(missing(deviceFilename)) {
        .Call("Forwarder_mr_ptforward", ptr, pi, A, B, noStates, alphabetSize, "-")
      } else {
        .Call("Forwarder_mr_ptforward", ptr, pi, A, B, noStates, alphabetSize, deviceFilename)
      }
      
    },

    getOrigSeqLength = function() {
      return(.Call("Forwarder_get_orig_seq_length", ptr))
    },

    getOrigAlphabetSize = function() {
      return(.Call("Forwarder_get_orig_alphabet_size", ptr))
    },

    getSeqLength = function(noStates) {
      if(missing(noStates))
        stop("noStates is missing.")
      if(!is.numeric(noStates))
        stop("noStates should be a number")
      noStates = as.integer(noStates)
      
      return(.Call("Forwarder_get_seq_length", ptr, noStates))
    },

    getAlphabetSize = function(noStates) {
      if(missing(noStates))
        stop("noStates is missing.")
      if(!is.numeric(noStates))
        stop("noStates should be a number")
      noStates = as.integer(noStates)
      
      return(.Call("Forwarder_get_alphabet_size", ptr, noStates))
    },

    getPair = function(symbol) {
      if(missing(symbol))
        stop("symbol is missing.")
      if(!is.numeric(symbol))
        stop("symbol should be a number")
      symbol = as.integer(symbol)
      
      return(.Call("Forwarder_get_pair", ptr, symbol))
    }
  )

)

viterbi <- function(seqFilename, hmm) {
  if(missing(seqFilename))
    stop("seqFilename is missing.")

  if(missing(hmm))
    stop("hmm is missing")
  
  if(!is.list(hmm) || !("pi" %in% names(hmm) && "A" %in% names(hmm) && "B" %in% names(hmm)))
    stop("hmm should be a list of the form list(\"pi\" = ..., \"A\" = ..., \"B\" = ...)")
  
  pi = hmm$pi
  A = hmm$A
  B = hmm$B
  
  if(!is.vector(pi))
    stop("hmm$pi should be a vector.")
  if(!is.matrix(A))
    stop("hmm$A should be a matrix.")
  if(!is.matrix(B))
    stop("hmm$A should be a matrix.")
    
  noStates = dim(B)[1]
  alphabetSize = dim(B)[2]

  res = .Call("c_viterbi", seqFilename, pi, A, B, noStates, alphabetSize)

  return(list("loglik" = res[[1]], "path" = res[[2]]))
}

posteriorDecoding <- function(seqFilename, hmm) {
  if(missing(seqFilename))
    stop("seqFilename is missing.")

  if(missing(hmm))
    stop("hmm is missing")
  
  if(!is.list(hmm) || !("pi" %in% names(hmm) && "A" %in% names(hmm) && "B" %in% names(hmm)))
    stop("hmm should be a list of the form list(\"pi\" = ..., \"A\" = ..., \"B\" = ...)")
  
  pi = hmm$pi
  A = hmm$A
  B = hmm$B
  
  if(!is.vector(pi))
    stop("hmm$pi should be a vector.")
  if(!is.matrix(A))
    stop("hmm$A should be a matrix.")
  if(!is.matrix(B))
    stop("hmm$A should be a matrix.")
  
  noStates = dim(B)[1]
  alphabetSize = dim(B)[2]

  res = .Call("c_posterior_decoding", seqFilename, pi, A, B, noStates, alphabetSize)

  return(list("path" = res[[1]], "table" = res[[2]]))
}

calibrate <- function(deviceFilename) {
  if(missing(deviceFilename))
    .Call("c_calibrate", deviceFilename)
  else
    .Call("c_calibrate", "-")

  cat()
}
