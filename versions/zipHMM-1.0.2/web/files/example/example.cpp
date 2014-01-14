#include "zipHMM/hmm_io.hpp"
#include "zipHMM/forwarder.hpp"
#include "zipHMM/matrix.hpp"
#include "zipHMM/viterbi.hpp"
#include "zipHMM/posterior_decoding.hpp"

#include <iostream>
#include <vector>

using namespace zipHMM;

int main(int argc, char **args) {
  
  Forwarder f;
  size_t alphabet_size = 3;
  f.read_seq("example.seq", alphabet_size);
  f.write_to_directory("example_out");

  Matrix pi, A, B;
  read_HMM(pi, A, B, "example.hmm");

  std::cout << "log-likelihood: " << f.forward(pi, A, B) << std::endl;
  // std::cout << "log-likelihood: " << f.pthread_forward(pi, A, B)  << std::endl; # parallelized version


  std::vector<unsigned> pd_path;
  Matrix pd_table;
  posterior_decoding("example.seq", pi, A, B, pd_path, pd_table);
  std::cout << "posterior path[0:10]: ";
  for(size_t i = 0; i < 10; ++i)
    std::cout << pd_path[i] << " ";
  std::cout << std::endl;
  
  std::cout << "posterior table column 0 - 9: " << std::endl;
  for(size_t r = 0; r < 2; ++r) {
    for(size_t c = 0; c < 10; ++c) {
      std::cout << pd_table(r, c) << "\t";
    }
    std::cout << std::endl;
  }

  std::vector<unsigned> viterbi_path;
  double viterbi_ll = viterbi("example.seq", pi, A, B, viterbi_path);
  std::cout << "viterbi loglikelihood: " << viterbi_ll << std::endl;
  std::cout << "viterbi path[0:10]: ";
  for(size_t i = 0; i < 10; ++i)
    std::cout << viterbi_path[i] << " ";
  std::cout << std::endl;
  
  return 0;
}
