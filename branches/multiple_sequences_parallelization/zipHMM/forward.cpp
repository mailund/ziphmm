#include "forwarder.hpp"
#include "hmm_io.hpp"

#include <iostream>
#include <unistd.h>
#include <getopt.h>
#include <vector>
#include <algorithm>

static const struct option opts[] = {
  { "sequence", required_argument, NULL, 's' },
  { "model", required_argument, NULL, 'm'},
  { "output", required_argument, NULL, 'o' },
  { "directory", required_argument, NULL, 'd'},
  { "pthread", no_argument, NULL, 'p'},
  { "help", no_argument, NULL, 'h' },
  { NULL, no_argument, NULL, 0 }
};

static const char *opt_string = "s:?m:?o:?d:?p?h?";

void print_usage(char *command) {
  std::cout << "Usage: " << command << " (-s <sequence filename> -m <HMM filename> [-e #expected forward calls] [-o <output directory>] ) | (-d <preprocessing directory> -m <HMM filename>) [-p]";
  std::cout << std::endl;
}

int main(int argc, char **args) {
  std::string sequence_filename;
  std::string preproc_dir;
  std::string hmm_filename;
  size_t minNoEvals = 0;
  std::string output_directory;
  bool parallel = false;

  if(argc < 3) {
    print_usage(args[0]);
    exit(-1);
  }

  int index = 0;
  int opt = getopt_long(argc, args, opt_string, opts, &index);
  while(opt != -1) {
    switch( opt ) {
    case 's':
      sequence_filename = optarg;
      break;
    case 'm':
      hmm_filename = optarg;
      break;
    case 'e':
      minNoEvals = size_t(atoi(optarg));
    case 'o':
      output_directory = optarg;
      break;
    case 'd':
      preproc_dir = optarg;
      break;
    case 'p':
      parallel = true;
      break;
    case 'h':
      print_usage(args[0]);
      exit(-1);
      break;
    default:
      std::cout << "Unrecognized command line parameter." << std::endl;
      print_usage(args[0]);
      exit(-1);
      break;
    }
    opt = getopt_long(argc, args, opt_string, opts, &index);
  }

  if(hmm_filename.empty() || (sequence_filename.empty() && preproc_dir.empty()) || (!sequence_filename.empty() && !preproc_dir.empty())) {
    print_usage(args[0]);
    exit(-1);
  }

  zipHMM::Matrix pi, A, B;
  zipHMM::read_HMM(pi, A, B, hmm_filename);
  size_t no_states = pi.get_height();

  double loglikelihood = 0.0;
  if(!sequence_filename.empty()) {
    zipHMM::Forwarder forwarder;
    if(minNoEvals != 0)
      forwarder.read_seq(sequence_filename, no_states, minNoEvals);
    else 
      forwarder.read_seq(sequence_filename, no_states);

    if(parallel)
      loglikelihood = forwarder.pthread_forward(pi, A, B);
    else
      loglikelihood = forwarder.forward(pi, A, B);

    if(!output_directory.empty())
      forwarder.write_to_directory(output_directory);

  } else {
    zipHMM::Forwarder forwarder;
    forwarder.read_from_directory(preproc_dir, no_states);
    
    if(parallel)
      loglikelihood = forwarder.pthread_forward(pi, A, B);
    else
      loglikelihood = forwarder.forward(pi, A, B);
  }
  
  std::cout.precision(12);
  std::cout << loglikelihood << std::endl;
}
