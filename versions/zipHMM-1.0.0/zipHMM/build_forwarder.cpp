#include "forwarder.hpp"

#include <iostream>
#include <unistd.h>
#include <getopt.h>
#include <vector>
#include <algorithm>


static const struct option opts[] = {
  { "sequence", required_argument, NULL, 's' },
  { "alphabet-size", required_argument, NULL, 'M'},
  { "no-states", required_argument, NULL, 'N' },
  { "output", required_argument, NULL, 'o' },
  { "help", no_argument, NULL, 'h' },
  { NULL, no_argument, NULL, 0 }
};

static const char *opt_string = "s:M:N:?o:h?";

void print_usage(char *command) {
  std::cout << "Usage: " << command << " -s <sequence filename> -M <alphabet size> -o <output directory> [-N <number of states>]*" << std::endl;
  std::cout << "See http://birc.au.dk/software/zipHMM for examples." << std::endl;
}

int main(int argc, char **args) {
  char *sequence_filename = NULL;   /* -s option */
  size_t alphabet_size = 0;         /* -M option */
  char *output_dirName = NULL;      /* -o option */
  std::vector<size_t> nStatesSave;  /* -N option */

  if(argc < 4) {
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
    case 'M':
      alphabet_size = (size_t) atol(optarg);
      break;
    case 'o':
      output_dirName = optarg;
      break;
    case 'N':
      nStatesSave.push_back( (size_t) atol(optarg) );
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

  if(sequence_filename == NULL)
    std::cout << "no sequence filename given." << std::endl;
  if(output_dirName == NULL)
    std::cout << "no output directory given." << std::endl;
  if(alphabet_size == 0)
    std::cout << "no alphabet size given." << std::endl;
  if(sequence_filename == NULL || alphabet_size == 0 || output_dirName == NULL) {
    print_usage(args[0]);
    exit(-1);
  }

  std::sort(nStatesSave.begin(), nStatesSave.end());

  std::cout << "Building forwarder:" << std::endl;
  std::cout << "-- sequence filename:\t" << sequence_filename << std::endl;
  if(!nStatesSave.empty()) {
    std::cout << "-- number of states:\t";
    for(std::vector<size_t>::iterator it = nStatesSave.begin(); it != nStatesSave.end(); ++it)
      std::cout << (*it) << ", ";
    std::cout << "\b\b" << std::endl;
  }
  std::cout << "-- alphabet size:\t" << alphabet_size << std::endl << std::endl;

  zipHMM::Forwarder f;
  f.read_seq(sequence_filename, alphabet_size, nStatesSave);
  f.write_to_directory(output_dirName);

  std::cout << "Forwarder saved to " << output_dirName << "/" << std::endl;


  return EXIT_SUCCESS;
}
