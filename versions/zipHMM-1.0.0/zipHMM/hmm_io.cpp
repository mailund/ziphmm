#include "hmm_io.hpp"
#include "io_utils.hpp"

#include <fstream>
#include <cstdlib>

namespace zipHMM {

  void read_HMM_spec(size_t &no_states, size_t &alphabet_size, const std::string &filename) {
    std::ifstream in(filename.c_str());
      
    if(!in) {
      std::cerr << "Unable to open \"" << filename << "\"" << std::endl;
      exit(-1);
    }

    read_HMM_spec(no_states, alphabet_size, in);

    in.close();
  }

  void read_HMM_spec(size_t &no_states, size_t &alphabet_size, std::istream &in) {
    read_token_or_exit(in, "no_states");
    no_states = read_or_exit<size_t>(in, "number of states");
    
    read_token_or_exit(in, "alphabet_size");
    alphabet_size = read_or_exit<size_t>(in, "number of observables");
  }

  void read_HMM(Matrix &result_pi, Matrix &result_A, Matrix &result_B, const std::string &filename) {
    std::ifstream in(filename.c_str());
      
    if(!in) {
      std::cerr << "Unable to open \"" << filename << "\"" << std::endl;
      exit(-1);
    }

    read_HMM(result_pi, result_A, result_B, in);

    in.close();
  }

  void read_HMM(Matrix &result_pi, Matrix &result_A, Matrix &result_B, std::istream &in) {
    size_t no_states;
    size_t alphabet_size;
    read_HMM_spec(no_states, alphabet_size, in);
      
    read_token_or_exit(in, "pi");
    result_pi.reset(no_states, 1);
    for(size_t i = 0; i < no_states; ++i)
      result_pi(i, 0) = read_or_exit<double>(in, "initial probability");

    read_token_or_exit(in, "A");
    result_A.reset(no_states, no_states);
    for(size_t from = 0; from < no_states; ++from) {
      for(size_t to = 0; to < no_states; ++to)
	result_A(from, to) = read_or_exit<double>(in, "transition probability");
    }
      
    read_token_or_exit(in, "B");
    result_B.reset(no_states, alphabet_size);
    for(size_t i = 0; i < no_states; ++i) {
      for(size_t ob = 0; ob < alphabet_size; ++ob)
	result_B(i, ob) = read_or_exit<double>(in, "emission probability");
    }
  }

  void write_HMM(const Matrix &pi, const Matrix &A, const Matrix &B, const std::string &filename) {
    std::ofstream out(filename.c_str());
      
    if(!out) {
      std::cerr << "Unable to open \"" << filename << "\"" << std::endl;
      exit(-1);
    }
      
    write_HMM(pi, A, B, out);
      
    out.close();
  }

  void write_HMM(const Matrix &pi, const Matrix &A, const Matrix &B, std::ostream &out) {
    size_t no_states = pi.get_height();
    size_t alphabet_size = B.get_width();
      
    out << "no_states" << std::endl;
    out << no_states << std::endl;
      
    out << "alphabet_size" << std::endl;
    out << alphabet_size << std::endl;
      
    out << "pi" << std::endl;
    for(size_t i = 0; i < no_states; ++i)
      out << pi(i, 0) << std::endl;
      
    out << "A" << std::endl;
    for(size_t from = 0; from < no_states; ++from) {
      std::string sep = "";
      for(size_t to = 0; to < no_states; ++to) {
	out << sep << A(from, to);
	sep = " ";
      }
      out << std::endl;
    }

    out << "B" << std::endl;
    for(size_t i = 0; i < no_states; ++i) {
      std::string sep = "";
      for(size_t ob = 0; ob < alphabet_size; ++ob) {
	out << sep << B(i, ob);
	sep = " ";
      }
      out << std::endl;
    }
  }
}
