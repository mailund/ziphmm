#include "matrix.hpp"
#include "hmm_io.hpp"
#include "seq_io.hpp"

#include <cstdlib>
#include <sstream>
#include <fstream>
#include <iostream>
#include <string>
#include <cmath>
#include <ctime>

namespace {

  double random_prob() {
    return rand()/(RAND_MAX + 1.0);
  }

  unsigned sample_initial_state(zipHMM::Matrix &pi) {
    size_t no_states = pi.get_height();
    double prob = random_prob();
    double sum = 0.0;
    for(size_t i = 0; i < no_states; ++i) {
      sum += pi(0,i);
      if(sum > prob)
	return unsigned(i);
    }
    return unsigned(no_states-1);
  }

  unsigned sample_emission_from_state(zipHMM::Matrix &B, unsigned state) {
    size_t alphabet_size = B.get_width();
    double prob = random_prob();
    double sum = 0.0;
    for(size_t i = 0; i < alphabet_size; ++i) {
      sum += B(state, i);
      if(sum > prob)
	return unsigned(i);
    }
    return unsigned(alphabet_size-1);
  }

  unsigned sample_transition_from_state(zipHMM::Matrix &A, unsigned state) {
    size_t no_states = A.get_height();
    double prob = random_prob();
    double sum = 0.0;
    for(size_t i = 0; i < no_states; ++i) {
      sum += A(state, i);
      if(sum > prob)
	return unsigned(i);
    }
    return unsigned(no_states-1);
  }  
  
  double generate_seq(size_t length,
		      zipHMM::Matrix &pi, zipHMM::Matrix &A, zipHMM::Matrix &B, 
		      std::vector<unsigned> &obsseq, std::vector<unsigned> &hiddenseq) {

    double loglikelihood = 0.0;
  
    unsigned state;
    unsigned prev_state;
    unsigned emission;
  
    hiddenseq[0] = (state = sample_initial_state(pi));
    obsseq[0]    = (emission = sample_emission_from_state(B, state));
  
    loglikelihood += std::log( pi(state, 0) );
    loglikelihood += std::log( B(state, emission) );
  
    for(size_t i = 1; i < length; ++i) {
      prev_state   = state;
      hiddenseq[i] = (state = sample_transition_from_state(A, state));
      obsseq[i]    = (emission = sample_emission_from_state(B, state));
    
      loglikelihood += std::log( A(prev_state, state) );
      loglikelihood += std::log( B(state, emission) );
    }
  
    return loglikelihood;
  }

  double generate_seq(size_t length, zipHMM::Matrix &pi, zipHMM::Matrix &A, zipHMM::Matrix &B, std::ostream &obsseq_out, std::ostream &hiddenseq_out) {
  
    std::vector<unsigned> obsseq(length);
    std::vector<unsigned> hiddenseq(length);
  
    double loglikelihood = generate_seq(length, pi, A, B, obsseq, hiddenseq);
  
    for(size_t i = 0; i < length; ++i) {
      obsseq_out << obsseq[i] << " ";
      hiddenseq_out << hiddenseq[i] << " ";
    }
  
    return loglikelihood;
  }

  void help(const std::string cmd) {
    std::cout << "Usage: " << cmd << " "
	      << "hmmfilename length outputObsSeqFilename outputHiddenSeqFilename" << std::endl;
    exit(-1);
  }

}

int main(int argc, char **argv) {
  std::string hmm_filename;
  size_t length;
  std::string output_obsseq_filename;
  std::string output_hiddenseq_filename;
  
  zipHMM::Matrix pi, A, B;

  if(argc != 5)
    help(argv[0]);

  hmm_filename = argv[1];
  length = unsigned(atoi(argv[2]));
  output_obsseq_filename = argv[3];
  output_hiddenseq_filename = argv[4];
  
  read_HMM(pi, A, B, hmm_filename);
  
  std::ofstream output_obsseq_stream(output_obsseq_filename.c_str());
  std::ofstream output_hiddenseq_stream(output_hiddenseq_filename.c_str());
  
  if(!output_obsseq_stream.is_open()) {
    std::cerr << "Unable to open \"" << output_obsseq_filename << "\"" << std::endl;
    exit(-1);
  }
  if(!output_hiddenseq_stream.is_open()) {
    std::cerr << "Unable to open \"" << output_hiddenseq_filename << "\"" << std::endl;
    exit(-1);
  }

  srand ( unsigned(time(NULL)) );
  double likelihood = generate_seq(length, pi, A, B, output_obsseq_stream, output_hiddenseq_stream);
  
  output_obsseq_stream.close();
  output_hiddenseq_stream.close();
  
  std::cout << likelihood << std::endl;
  
  return 0;
}
