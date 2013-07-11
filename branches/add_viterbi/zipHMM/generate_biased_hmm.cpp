#include "hmm_io.hpp"

#include <cstdlib>
#include <sstream>
#include <fstream>
#include <iostream>
#include <string>
#include <cmath>
#include <time.h>

namespace {
    void make_random_prob_column_vector(zipHMM::Matrix &vector) {
    double sum = 0.0;
    for(unsigned i = 0; i < vector.get_height(); ++i) {
      sum += (vector(i,0) = rand() / (double)RAND_MAX);
    }
    for(unsigned i = 0; i < vector.get_height(); ++i) {
      vector(i,0) /= sum;
    }  
  }

  void make_random_row_prob_matrix(zipHMM::Matrix &m) {
    for(unsigned r = 0; r < m.get_height(); ++r) {
      double sum = 0.0;
      for(unsigned c = 0; c < m.get_width(); ++c) {
	sum += (m(r,c) = rand() / (double)RAND_MAX);
      }
      for(unsigned c = 0; c < m.get_width(); ++c) {
	m(r,c) /= sum;
      }
    }
  }

  void make_biased_row_prob_matrix(zipHMM::Matrix &m) {
    for(unsigned r = 0; r < m.get_height(); ++r) {
      double sum = 0.0;
      for(unsigned c = 0; c < m.get_width(); ++c) {
	sum += m(r,c) = (rand() / (double)RAND_MAX) * double(c+1);
      }
      for(unsigned c = 0; c < m.get_width(); ++c) {
	m(r,c) /= sum;
      }
    }
  }

  void generate_biased_hmm(zipHMM::Matrix &pi, zipHMM::Matrix &A, zipHMM::Matrix &B) {
    make_random_prob_column_vector(pi);
    make_random_row_prob_matrix(A);
    make_biased_row_prob_matrix(B);
  }

  void help(const std::string cmd) {
    std::cout << "Usage: " << cmd << " "
	      << "no_states alphabet_size hmm_filename"
	      << std::endl;
    exit(-1);
  }

}



int main(int argc, char **argv) {
  std::string hmm_filename;
  std::ofstream hmm_file_stream;
  unsigned no_states, alphabet_size;
  zipHMM::Matrix pi, A, B;

  if(argc != 4)
    help(argv[0]);

  no_states = unsigned(atoi(argv[1]));
  alphabet_size = unsigned(atoi(argv[2]));
  hmm_filename = argv[3];
  
  pi.reset(no_states, 1);
  A.reset(no_states, no_states);
  B.reset(no_states, alphabet_size);

  srand ( unsigned(time(NULL)) );
  generate_biased_hmm(pi, A, B);
  
  write_HMM(pi, A, B, hmm_filename);
  
  return 0;
}
