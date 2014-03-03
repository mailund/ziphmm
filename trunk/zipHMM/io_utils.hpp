#ifndef IO_UTILS_HPP
#define IO_UTILS_HPP

#include <cstdlib>
#include <iostream>

namespace zipHMM {

  const size_t MAX_PATH_LENGTH = 256;

  template<typename T>
  T read_or_exit(std::istream &in, const std::string &expected_desc) {
    T result;
    if(!(in >> result)) {
      std::cerr << "Error reading " << expected_desc << "." << std::endl;
      exit(-1);
    }
    return result;
  }
  
  bool read_token(std::istream &in, const std::string &token);

  void read_token_or_exit(std::istream &in, const std::string &expected_token);
  
  bool read_token_or_tell(std::istream &in, const std::string &expected_token);

  std::string get_user_dir();

  void mk_dir(const std::string &path);
  
  std::string get_working_directory();
}

#endif
