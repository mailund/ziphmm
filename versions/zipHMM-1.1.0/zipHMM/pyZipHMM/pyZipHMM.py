from ctypes import *
from distutils.sysconfig import get_python_lib
from os import path

try:
    d = path.dirname(__file__)
    lib = cdll.LoadLibrary("%s/libpyZipHMM.so" % (d))
    library_location = "%s/libpyZipHMM.so" % (d)
except OSError:
    python_lib = get_python_lib()
    lib = cdll.LoadLibrary(python_lib + "/libpyZipHMM.so")
    library_location = python_lib + "/libpyZipHMM.so"
except OSError as e:
    print "Error: pyZipHMM not found:"
    print "\t libpyZipHMM.so missing"
    print "Looked at:", python_lib, '/libpyZipHMM.so and ./libpyZipHMM.so'
    print "{0}: {1}".format(e.errno, e.strerror)
    exit(-1)
    

## HMM IO
def readHMMspec(filename):
    nStates = c_uint()
    nObservables = c_uint()
    lib.c_read_HMM_spec(byref(nStates), byref(nObservables), c_char_p(filename))
    return (nStates, nObservables)

def readHMM(filename):
    pi = Matrix()
    A = Matrix()
    B = Matrix()
    lib.c_read_HMM(pi.obj, A.obj, B.obj, c_char_p(filename))
    return (pi, A, B)

def writeHMM(pi, A, B, filename):
    lib.c_write_HMM(pi.obj, A.obj, B.obj, c_char_p(filename))

lib.c_read_seq.restype = py_object
def readSeq(filename):
    return lib.c_read_seq(filename)
    
## Forwarder
lib.Forwarder_new.restype = c_void_p
lib.Forwarder_forward.restype = c_double
lib.Forwarder_pthread_forward.restype = c_double
lib.Forwarder_mr_pthread_forward.restype = c_double
lib.Forwarder_get_orig_seq_length.restype = c_uint
lib.Forwarder_get_orig_alphabet_size.restype = c_uint
lib.Forwarder_get_seq_length.restype = c_uint
lib.Forwarder_get_alphabet_size.restype = c_uint
lib.Forwarder_get_pair.restype = py_object


class Forwarder(object):
    def __init__(self):
        self.obj = c_void_p(lib.Forwarder_new())

    @staticmethod
    def fromSequence(seqFilename, alphabetSize, nStatesSave = None, minNoEvals = 1):
        forwarder = Forwarder()

        if nStatesSave != None:
            arr = ( c_uint * len(nStatesSave) )()
            arr[:] = nStatesSave
            lib.Forwarder_read_seq(forwarder.obj, c_char_p(seqFilename), alphabetSize, arr, len(nStatesSave), minNoEvals)
        else:
            arr = ( c_uint * 0 )()
            lib.Forwarder_read_seq(forwarder.obj, c_char_p(seqFilename), alphabetSize, arr, 0, minNoEvals)
        
        return forwarder

    @staticmethod
    def fromSequenceDirectory(dirname, alphabetSize, nStatesSave = None, minNoEvals = 1):
        forwarder = Forwarder()

        if nStatesSave != None:
            arr = ( c_uint * len(nStatesSave) )()
            arr[:] = nStatesSave
            lib.Forwarder_read_seq_directory(forwarder.obj, c_char_p(dirname), alphabetSize, arr, len(nStatesSave), minNoEvals)
        else:
            arr = ( c_uint * 0 )()
            lib.Forwarder_read_seq_directory(forwarder.obj, c_char_p(dirname), alphabetSize, arr, 0, minNoEvals)
        
        return forwarder

    @staticmethod
    def fromDirectory(directory, nStates = None):
        forwarder = Forwarder()
        if nStates == None:
            lib.Forwarder_read_from_directory(forwarder.obj, c_char_p(directory))
        else:
            lib.Forwarder_read_from_directory(forwarder.obj, c_char_p(directory), nStates)
        return forwarder

    def __del__(self):
        from ctypes import cdll
        lib = cdll.LoadLibrary(library_location)
        lib.Forwarder_destructor(self.obj)

    def forward(self, pi, A, B):
        return lib.Forwarder_forward(self.obj, pi.obj, A.obj, B.obj)

    def ptforward(self, pi, A, B, device_filename = None):
        if device_filename == None:
            return lib.Forwarder_pthread_forward(self.obj, pi.obj, A.obj, B.obj, "-")
        else :
            return lib.Forwarder_pthread_forward(self.obj, pi.obj, A.obj, B.obj, device_filename)

    def mrforward(self, pi, A, B, device_filename = None):
        if device_filename == None:
            return lib.Forwarder_mr_pthread_forward(self.obj, pi.obj, A.obj, B.obj, "-")
        else :
            return lib.Forwarder_mr_pthread_forward(self.obj, pi.obj, A.obj, B.obj, device_filename)

    def getOrigSeqLength(self):
        return lib.Forwarder_get_orig_seq_length(self.obj)

    def getOrigAlphabetSize(self):
        return lib.Forwarder_get_orig_alphabet_size(self.obj)

    def getSeqLength(self, no_states):
        return lib.Forwarder_get_seq_length(self.obj, no_states)

    def getAlphabetSize(self, no_states):
        return lib.Forwarder_get_alphabet_size(self.obj, no_states)

    def getPair(self, symbol):
        return lib.Forwarder_get_pair(self.obj, symbol)
         
    def writeToDirectory(self, directory):
        lib.Forwarder_write_to_directory(self.obj, c_char_p(directory))


## SimpleForwarder
lib.SimpleForwarder_new.restype = c_void_p
lib.SimpleForwarder_forward.restype = c_double

class SimpleForwarder(object):
    def __init__(self, seqFilename):
        self.obj = c_void_p(lib.SimpleForwarder_new(seqFilename))

    def forward(self, pi, A, B):
        return lib.SimpleForwarder_forward(self.obj, pi.obj, A.obj, B.obj)
 
## Sequence
lib.Sequence_new.restype = c_void_p
lib.Sequence_destructor.restype = c_int
lib.Sequence_get.restype = c_uint
lib.Sequence_len.restype = c_uint

class Sequence(object):

    def __init__(self):
        self.obj = c_void_p(lib.Sequence_new())

    def __del__(self):
        from ctypes import cdll
        lib = cdll.LoadLibrary(library_location)
        lib.Sequence_destructor(self.obj)

    def __len__(self):
        return lib.Sequence_len(self.obj)

    def __getitem__(self, key):
        if isinstance(key, slice) :
            return [ self[ii] for ii in xrange(*key.indices(len(self))) ]
        elif isinstance( key, int ) :
            if key < 0 :
                key += len( self )
            if key >= len( self ) :
                raise IndexError, "The index (%d) is out of range." % key
            return int(lib.Sequence_get(self.obj, key))
        else:
            raise TypeError, "Invalid argument type."

    def __str__(self):
        string = ""
        for i in xrange(len(self)):
            string = string + ("%d" % self[i]) + " "
        return string.strip()
        
## Matrix
lib.Matrix_new_empty.restype = c_void_p
lib.Matrix_new_height_width.restype = c_void_p
lib.Matrix_get_width.restype = c_uint
lib.Matrix_get_height.restype = c_uint
lib.Matrix_get.restype = c_double

class Matrix(object):

    def __init__(self, height = 0, width = 0):
        if height == 0 or width == 0:
            self.obj = c_void_p(lib.Matrix_new_empty())
        else:
            self.obj = c_void_p(lib.Matrix_new_height_width(height, width))

    def __del__(self):
        from ctypes import cdll
        lib = cdll.LoadLibrary(library_location)
        lib.Matrix_destructor(self.obj)

    def getWidth(self):
        return lib.Matrix_get_width(self.obj)

    def getHeight(self):
        return lib.Matrix_get_height(self.obj)

    def reset(self, height, width):
        lib.Matrix_reset(self.obj, c_uint(height), c_uint(width))

    def __setitem__(self, (row, column), value):
        lib.Matrix_set(self.obj, c_uint(row), c_uint(column), c_double(value))

    def __getitem__(self, (row, column)):
        return lib.Matrix_get(self.obj, row, column)

    @staticmethod
    def transpose(f, t):
        lib.Matrix_transpose(f.obj, t.obj)

    def p(self):
        lib.Matrix_print(self.obj)

## posterior decoding
def posteriorDecoding(seqFilename, pi, A, B):
    pdTable = Matrix()
    pdPath = Sequence()
    lib.c_posterior_decoding(pdPath.obj, pdTable.obj, pi.obj, A.obj, B.obj, c_char_p(seqFilename))
    return pdPath, pdTable

## Viterbi
lib.c_viterbi.restype = c_double

def viterbi(seqFilename, pi, A, B):
    viterbiPath = Sequence()
    viterbi_ll = lib.c_viterbi(viterbiPath.obj, pi.obj, A.obj, B.obj, c_char_p(seqFilename))
    return viterbiPath, viterbi_ll
        
## calibrate
def calibrate(deviceFilename = None):
    if deviceFilename == None:
        lib.c_calibrate("-")
    else:
        lib.c_calibrate(deviceFilename)


if __name__ == "__main__":
    print "Constructing Matrix(3,7)"
    m = Matrix(3, 7)
    print "Calling getHeight()"
    assert m.getHeight() == 3
    print "Calling getWidth()"
    assert m.getWidth() == 7
    print "Calling setitem method"
    m[1,2] = 0.5
    print "Calling getitem method"
    assert m[1, 2] == 0.5
    print "Calling reset method"
    m.reset(7,3)
    assert m.getHeight() == 7
    assert m.getWidth() == 3

    print "Calling readHMM method"
    (pi, A, B) = readHMM("test_data/test1.hmm")
    assert pi.getHeight() == 2
    assert pi.getWidth()  == 1
    assert A.getHeight()  == 2
    assert A.getWidth()   == 2
    assert B.getHeight()  == 2
    assert B.getWidth()   == 2

    print "Creating Forwarder object from files"
    f = Forwarder(newSeqFilename = "../new_seq.tmp", dataStructureFilename = "../data_structure.tmp")
    assert f.getOrigAlphabetSize() == 2
    assert f.getOrigSeqLength()    == 18
    assert f.getNewAlphabetSize()  == 4
    print "Calling forward on Forwarder object"
    assert abs(f.forward(pi, A, B)  - -12.5671022728) < 0.001

    print "Calling readHMMspec method"
    (nStates, nObservables) = readHMMspec("test_data/test1.hmm")
    assert nStates.value == 2
    assert nObservables.value == 2

    print "Creating Forwarder from sequence and hmm spec"
    f = Forwarder(seqFilename = "test_data/test1.seq", nStates = nStates, nObservables = nObservables)
    assert f.getOrigAlphabetSize() == 2
    assert f.getOrigSeqLength()    == 18
    assert f.getNewAlphabetSize()  == 4
    print "Calling forward"
    assert abs(f.forward(pi, A, B) - -12.5671022728) < 0.001
