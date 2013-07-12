#ifndef VITERBI_HPP
#define VITERBI_HPP

#include "matrix.hpp"
#include "seq_io.hpp"

#include <vector>
#include <algorithm>
#include <limits>
#include <cmath>

namespace zipHMM {
  
  double viterbi(const std::vector<unsigned> &seq,
		 const Matrix &pi,
		 const Matrix &A,
		 const Matrix &B,
		 std::vector<unsigned> &viterbi_path) {

    size_t no_states = A.get_height();
    size_t length = seq.size();

    Matrix viterbi_table(no_states, length);

    //init
    for(size_t r = 0; r < no_states; ++r) {
      viterbi_table(r, 0) = std::log(pi(r, 0) * B(r, seq[0]));
    }

    // recursion
    for(size_t c = 1; c < length; ++c) {
      for(size_t r = 0; r < no_states; ++r) {
	double max_value = -std::numeric_limits<double>::max();
	for(size_t prev_state = 0; prev_state < no_states; ++prev_state) {
	  max_value = std::max(max_value, viterbi_table(prev_state, c - 1) + std::log(A(prev_state, r) * B(r, seq[c])));
	}
	viterbi_table(r,c) = max_value;
      }
    }

    double path_ll = viterbi_table(0, length - 1);
    size_t end_point = 0;
    for(size_t r = 1; r < no_states; ++r) {
      if(viterbi_table(r, length - 1) > path_ll) {
	path_ll = viterbi_table(r, length - 1);
	end_point = r;
      }
    }


    // backtrack
    viterbi_path.resize(length);
    viterbi_path[length - 1] = unsigned(end_point);

    size_t current_state = end_point;
    for(size_t c = length - 1; c > 0; --c) {
      double max_value = viterbi_table(0, c - 1) + std::log(A(0, current_state) * B(current_state, seq[c]));
      size_t max_state = 0;
      for(size_t prev_state = 1; prev_state < no_states; ++prev_state) {
	double prev_state_ll = viterbi_table(prev_state, c - 1) + std::log(A(prev_state, current_state) * B(current_state, seq[c]));
	if(prev_state_ll > max_value) {
	  max_value = prev_state_ll;
	  max_state = prev_state;
	}
      }
      viterbi_path[c - 1] = unsigned(max_state);
      current_state = max_state;
    }

    return path_ll;
  }

  double viterbi(const std::string &seq_filename,
		 const Matrix &pi,
		 const Matrix &A,
		 const Matrix &B,
		 std::vector<unsigned> &viterbi_path) {
    std::vector<unsigned> seq;
    readSeq(seq, seq_filename);
    return viterbi(seq, pi, A, B, viterbi_path);
  }

}

#endif
