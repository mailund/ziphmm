#include "hmm_io.hpp"
#include "seq_io.hpp"
#include "simple_stop_forwarder.hpp"
#include "forwarder.hpp"
#include "simple_forwarder.hpp"
#include "prob_spaces.hpp"
#include "timer.hpp"
#include "test_utils.hpp"


using namespace zipHMM;

void test_forwarder(const std::string &hmmFilename, const std::string &seqFilename, size_t minNumEvals, double expected_ll) {
  Matrix pi, A, B;
  size_t no_states, alphabet_size;

  std::cout << "Testing " << hmmFilename << " " << seqFilename << ": ";
  std::cout.flush();

  read_HMM(pi, A, B, hmmFilename);
  no_states = pi.get_height();
  alphabet_size = B.get_width();

  Forwarder forwarder;
  forwarder.read_seq(seqFilename, alphabet_size, minNumEvals);
  SimpleForwarder s_forwarder(seqFilename);
  
  assertClose(forwarder.forward(pi, A, B), expected_ll, "Forwarder: Difference in likelihood log likelihood:", 0.0001);
  assertClose(forwarder.pthread_forward(pi, A, B), expected_ll, "pthread_forward: Difference in likelihood log likelihood:", 0.0001);
  assertClose(forwarder.pthread_forward_par_stage1(pi, A, B), expected_ll, "pthread_forward_orig: Difference in likelihood log likelihood:", 0.0001);
  assertClose(s_forwarder.forward(pi, A, B), expected_ll, "s_forward: Difference in log likelihood:", 0.0001);

  std::cout << "ok." << std::endl;
}

void test_forwarder_long_sequence(const std::string &hmmFilename, const std::string &seqFilename, double expected_ll) {
 Matrix pi, A, B;
  size_t no_states, alphabet_size;

  std::cout << "Testing " << hmmFilename << " " << seqFilename << ": ";
  std::cout.flush();

  read_HMM(pi, A, B, hmmFilename);
  no_states = pi.get_height();
  alphabet_size = B.get_width();

  Forwarder forwarder;
  forwarder.read_seq(seqFilename, alphabet_size, 500);
  SimpleForwarder s_forwarder(seqFilename);
  SimpleStopForwarder ss_forwarder;
  ss_forwarder.read_seq(seqFilename, alphabet_size, no_states);
      
  assertClose(forwarder.forward(pi, A, B), expected_ll, "Forwarder: Difference in likelihood log likelihood:", 0.0001);
  assertClose(forwarder.pthread_forward(pi, A, B), expected_ll, "pthread_forward: Difference in likelihood log likelihood:", 0.0001);
  assertClose(forwarder.pthread_forward_par_stage1(pi, A, B), expected_ll, "pthread_forward_orig: Difference in likelihood log likelihood:", 0.0001);
  assertClose(ss_forwarder.forward(pi, A, B), expected_ll, "t_forward: Difference in log likelihood:", 0.0001);
  assertClose(s_forwarder.forward(pi, A, B), expected_ll, "s_forward: Difference in log likelihood:", 0.0001);

  std::cout << "ok." << std::endl;
}

void test_forwarder_io() {
  const std::string &directory = "test_data/data_structure_test";
  std::cout << "Testing " << directory << ": ";

  Forwarder forwarder1a;
  forwarder1a.read_from_directory(directory);
  Forwarder forwarder1b;
  forwarder1b.read_from_directory(directory, 16);

  assertEqual<size_t>(forwarder1a.get_orig_seq_length(), 10000);
  assertEqual<size_t>(forwarder1a.get_orig_alphabet_size(), 3);
  assertEqual<size_t>(forwarder1a.get_alphabet_size(16), 34);
  assertEqual<size_t>(forwarder1a.get_alphabet_size(32), 33);
  assertEqual<unsigned>(forwarder1a.get_pair(20).first, 0);
  assertEqual<unsigned>(forwarder1a.get_pair(20).second, 4);
  assertEqual<size_t>(forwarder1a.get_seq_length(16), 3457);
  assertEqual<size_t>(forwarder1a.get_seq_length(32), 3489);

  assertEqual<size_t>(forwarder1b.get_orig_seq_length(), 10000);
  assertEqual<size_t>(forwarder1b.get_orig_alphabet_size(), 3);
  assertEqual<size_t>(forwarder1b.get_alphabet_size(16), 34);
  assertEqual<unsigned>(forwarder1b.get_pair(20).first, 0);
  assertEqual<unsigned>(forwarder1b.get_pair(20).second, 4);
  assertEqual<size_t>(forwarder1b.get_seq_length(16), 3457);

  forwarder1a.write_to_directory(directory + "_post");  
  assertFilesEqual(directory + "/data_structure", directory + "_post/data_structure");
  
  Forwarder forwarder2a;
  forwarder2a.read_from_directory(directory + "_post");
  Forwarder forwarder2b;
  forwarder2b.read_from_directory(directory + "_post", 16);
  Matrix pi, A, B;
  read_HMM(pi, A, B, "test_data/test6.hmm");

  assertClose(forwarder1a.forward(pi, A, B), -10457.5211, "forwarder1 likelihood", 0.0001);
  assertClose(forwarder1b.forward(pi, A, B), -10457.5211, "forwarder1b likelihood", 0.0001);
  assertClose(forwarder2b.forward(pi, A, B), -10457.5211, "forwarder2b likelihood", 0.0001);
  
  std::cout << "ok." << std::endl;
}


int main(int argc, char **argv) {
  std::cout << "Testing library. Expect it to take 1-2 minutes" << std::endl;
  std::cout << "-- Testing forwarder.forward:" << std::endl;
  test_forwarder("test_data/test0.hmm",             "test_data/test0.seq",             1,   -12.4767);
  test_forwarder("test_data/test1.hmm",             "test_data/test1.seq",             1,   -12.5671);
  test_forwarder("test_data/test2.hmm",             "test_data/test2.seq",             1,   -27.531);
  test_forwarder("test_data/test2.hmm",             "test_data/test2_10000.seq",       1,   -13742.6069);
  test_forwarder("test_data/test3.hmm",             "test_data/test3.seq",             1,   LogSpace::ONE);
  test_forwarder("test_data/test4.hmm",             "test_data/test4.seq",             1,   -24.0953);
  test_forwarder("test_data/test4.hmm",             "test_data/test4_10000.seq",       1,   -44239.73129);
  test_forwarder("test_data/test5.hmm",             "test_data/test5.seq",             1,   -5.78270);
  test_forwarder("test_data/test6.hmm",             "test_data/test6.seq",             1,   -18.8829);
  test_forwarder("test_data/test6.hmm",             "test_data/test6_10000.seq",       500,   -10457.5211);
  test_forwarder("test_data/test6.hmm",             "test_data/test6_100000.seq",      500,   -104126.2410);
  test_forwarder("test_data/test8States.hmm",       "test_data/test8States.seq",       1,   LogSpace::ONE);
  test_forwarder("test_data/testOffByOneEm.hmm",    "test_data/testOffByOneEm.seq",    1,   -16.6355);
  test_forwarder("test_data/testOffByOneTrans.hmm", "test_data/testOffByOneTrans.seq", 1,   LogSpace::ONE);
  test_forwarder("test_data/testOneState.hmm",      "test_data/testOneState.seq",      1,   -12.4766);
  test_forwarder("test_data/testVectorReserve.hmm", "test_data/testVectorReserve.seq", 1,   -2.41869);
  test_forwarder("test_data/testZeroAdd.hmm",       "test_data/testZeroAdd.seq",       1,   LogSpace::ZERO);
  test_forwarder("test_data/testNoSolution.hmm",    "test_data/testNoSolution.seq",    1,   LogSpace::ZERO);

  std::cout << "-- Testing io:" << std::endl;
  test_forwarder_io();
}
