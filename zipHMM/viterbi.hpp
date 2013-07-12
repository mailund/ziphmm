#ifndef VITERBI_HPP
#define VITERBI_HPP

#include "matrix.hpp"

#include <vector>

namespace zipHMM {
  
  double viterbi(const std::vector<unsigned> &seq,
		 const Matrix &initProbs,
		 const Matrix &transProbs,
		 const Matrix &emProbs,
		 std::vector<unsigned> &viterbi_path);

  double viterbi(const std::string &seq_filename,
		 const Matrix &pi,
		 const Matrix &A,
		 const Matrix &B,
		 std::vector<unsigned> &viterbi_path);

}

#endif
