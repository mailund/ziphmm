#include "zipHMM/hmm_io.hpp"
#include "zipHMM/forwarder.hpp"
#include "zipHMM/matrix.hpp"

#include <iostream>

using namespace zipHMM;

int main(int argc, char **args) {
  Forwarder f;
  size_t alphabet_size = 3;
  f.read_seq_directory("sequences", alphabet_size);
  f.write_to_directory("example_out");

  Matrix pi, A, B;
  read_HMM(pi, A, B, "example.hmm");

  std::cout << "log-likelihood: " << f.forward(pi, A, B) << std::endl;
  // std::cout << "log-likelihood: " << f.pthread_forward(pi, A, B)  << std::endl; // parallelized on the lengths of the sequences
  // std::cout << "log-likelihood: " << f.mr_pthread_forward(pi, A, B)  << std::endl; // parallelized version on the number of sequences
  return 0;
}
