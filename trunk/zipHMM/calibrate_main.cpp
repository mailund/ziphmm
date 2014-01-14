#include "calibrate.hpp"

#include <string>

void help(const std::string &execName) {
        std::cout << "Usage: " << execName << " [filename]" << std::endl;
        exit(-1);
}

int main(int argc, char **argv) {
  std::string filename = zipHMM::DEFAULT_DEVICE_FILENAME;
  switch(argc) {
  case 1:
    break;
  case 2:
    filename = argv[1];
    break;
  default:
    help(argv[0]);
  }
  
  zipHMM::calibrate(filename);
}
