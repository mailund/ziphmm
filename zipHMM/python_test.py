from pyZipHMM import *
import unittest
import logging
import sys

class TestSequenceFunctions(unittest.TestCase):

    def test_forwarder_io(self):
        log = logging.getLogger( "python_test")
        log.debug("Testing test_forwarder_io: ")

        directory = "../test_data/data_structure_test"
        
        forwarder1a = Forwarder.fromDirectory(directory = directory)
        forwarder1b = Forwarder.fromDirectory(directory = directory, nStates = 16)

        self.assertEqual(forwarder1a.getOrigSeqLength(), 10000);
        self.assertEqual(forwarder1a.getOrigAlphabetSize(), 3);
        self.assertEqual(forwarder1a.getAlphabetSize(16), 34);
        self.assertEqual(forwarder1a.getAlphabetSize(32), 33);
        # self.assertEqual(forwarder1a.getPair(20).first, 0);
        # self.assertEqual(forwarder1a.getPair(20).second, 4);
        self.assertEqual(forwarder1a.getSeqLength(16), 3457);
        self.assertEqual(forwarder1a.getSeqLength(32), 3489);
        
        self.assertEqual(forwarder1b.getOrigSeqLength(), 10000);
        self.assertEqual(forwarder1b.getOrigAlphabetSize(), 3);
        self.assertEqual(forwarder1b.getAlphabetSize(16), 34);
        self.assertEqual(forwarder1b.getAlphabetSize(32), 33);
        # self.assertEqual(forwarder1b.getPair(20).first, 0);
        # self.assertEqual(forwarder1b.getPair(20).second, 4);
        self.assertEqual(forwarder1b.getSeqLength(16), 3457);
        self.assertEqual(forwarder1b.getSeqLength(32), 3489);

        forwarder1a.writeToDirectory(directory + "_post");
        
        with open(directory + "/data_structure") as f1, open(directory + "_post/data_structure") as f2:
            self.assertEqual(f1.read(), f2.read())
            
        forwarder2a = Forwarder.fromDirectory(directory = directory + "_post")
        forwarder2b = Forwarder.fromDirectory(directory = directory + "_post", nStates = 16)
        
        pi, A, B = readHMM("../test_data/test6.hmm");
        
        self.assertAlmostEqual(forwarder1a.forward(pi, A, B), -10457.5211, delta=0.0001);
        self.assertAlmostEqual(forwarder1b.forward(pi, A, B), -10457.5211, delta=0.0001);
        self.assertAlmostEqual(forwarder2a.forward(pi, A, B), -10457.5211, delta=0.0001);
        self.assertAlmostEqual(forwarder2b.forward(pi, A, B), -10457.5211, delta=0.0001);

        log.debug("ok.\n")

    def test_forwarder(self):
        log = logging.getLogger( "python_test")
        log.debug("Testing test_forwarder: ")

        f = Forwarder.fromSequence(seqFilename = "../test_data/test6_10000.seq", alphabetSize = 3, nStatesSave = [2, 4, 8, 16, 32, 64], minNoEvals = 10)

        pi, A, B = readHMM("../test_data/test6.hmm")
    
        self.assertAlmostEqual(f.forward(pi, A, B), -10457.5210537, delta = 0.0001)
        self.assertAlmostEqual(f.ptforward(pi, A, B), -10457.5210537, delta = 0.0001)
        self.assertAlmostEqual(f.mrforward(pi, A, B), -10457.5210537, delta = 0.0001)
    
        log.debug("ok.\n")

    def test_forwarder_multiple_seqs(self):
        log = logging.getLogger( "python_test")
        log.debug("Testing test_forwarder: ")

        f = Forwarder.fromSequenceDirectory(dirname = "../test_data/seqs_directory" , alphabetSize = 3, nStatesSave = [2, 4, 8, 16, 32, 64], minNoEvals = 10)

        pi, A, B = readHMM("../test_data/multiple_seqs.hmm")

        self.assertAlmostEqual(f.forward(pi, A, B),   -52287.605268, delta = 0.0001)
        self.assertAlmostEqual(f.ptforward(pi, A, B), -52287.605268, delta = 0.0001)
        self.assertAlmostEqual(f.mrforward(pi, A, B), -52287.605268, delta = 0.0001)

        log.debug("ok.\n")
        

    def test_SimpleForwarder(self):
        log = logging.getLogger( "python_test")
        log.debug("Testing test_SimpleForwarder: ")
       
        f = SimpleForwarder(seqFilename = "../test_data/test6_10000.seq")

        pi, A, B = readHMM("../test_data/test6.hmm")
    
        self.assertAlmostEqual(f.forward(pi, A, B), -10457.5210537, delta = 0.0001)
    
        log.debug("ok.\n")

    def test_nStatesSave_constructor(self):
        log = logging.getLogger( "python_test")
        log.debug("Testing test_nStatesSave_constructor: ")
       
        f = Forwarder.fromSequence(seqFilename = "../test_data/test6_10000.seq", alphabetSize = 3, nStatesSave = [2, 4, 8, 16, 32, 64])
        pi, A, B = readHMM("../test_data/test6.hmm")
    
        self.assertAlmostEqual(f.forward(pi, A, B), -10457.5210537, delta = 0.0001)
        self.assertAlmostEqual(f.ptforward(pi, A, B), -10457.5210537, delta = 0.0001)
    
        log.debug("ok.\n")

    def test_viterbi(self):
        log = logging.getLogger( "python_test")
        log.debug("Testing test_viterbi: ")


        pi, A, B = readHMM("../test_data/test1.hmm")
        viterbiPath, viterbi_ll = viterbi("../test_data/test1.seq", pi, A, B)

        self.assertAlmostEqual(viterbi_ll,  -15.4332, delta = 0.0001)
        self.assertEqual(viterbiPath[1], 1)
        self.assertEqual(viterbiPath[-1], 1)

        pi, A, B = readHMM("../test_data/test4.hmm")
        viterbiPath, viterbi_ll = viterbi("../test_data/test4.seq", pi, A, B)

        self.assertAlmostEqual(viterbi_ll,  -24.7258, delta = 0.0001)
        self.assertEqual(viterbiPath[1], 1)
        self.assertEqual(viterbiPath[-1], 3)

        log.debug("ok.\n")
        

if __name__ == "__main__":
    logging.basicConfig( stream=sys.stderr )
    logging.getLogger( "python_test").setLevel( logging.DEBUG )
    unittest.main()
