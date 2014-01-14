#include <cstdlib>
#include <sstream>
#include <fstream>
#include <iostream>
#include <string>
#include <cmath>
#include <time.h>

namespace {
  double random_prob() {
    return rand()/(RAND_MAX + 1.0);
  }

  void generate_01_seq(unsigned length, double one_freq, std::ofstream &output_stream) {
    for(unsigned i = 0; i < length; ++i) {
      if(random_prob() < one_freq)
	output_stream << 1 << " ";      
      else
	output_stream << 0 << " ";
    }
  }

  void help(const std::string cmd) {
    std::cout << "Usage: " << cmd << " "
	      << "length one_freq output_filename"
	      << std::endl;
    exit(-1);
  }

}


int main(int argc, char **argv) {
  unsigned length;
  double one_freq;
  std::string output_filename;
  
  if(argc != 4)
    help(argv[0]);

  length = unsigned(atoi(argv[1]));
  one_freq = atof(argv[2]);
  output_filename = argv[3];
  
  std::ofstream output_stream(output_filename.c_str());
  
  if(!output_stream.is_open()) {
    std::cerr << "Unable to open \"" << output_filename << "\"" << std::endl;
    exit(-1);
  }

  srand ( unsigned(time(NULL)) );
  generate_01_seq(length, one_freq, output_stream);
  
  output_stream.close();
  
  return 0;
}
