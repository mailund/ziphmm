#include "posterior_decoding.hpp"
#include "matrix.hpp"
#include "seq_io.hpp"

#include <vector>

namespace zipHMM {

  void forward(Matrix &forward_table, 
	       std::vector<double> &scales,
	       const Matrix &pi, const Matrix &A, const Matrix&B, 
	       const std::vector<unsigned> &seq) {

    size_t no_states = A.get_height();
    size_t length = seq.size();

    scales.resize(length);
    forward_table.reset(no_states, length);

    //init
    double sum = 0.0;
    for(size_t r = 0; r < no_states; ++r) {
      forward_table(r, 0) = pi(r, 0) * B(r, seq[0]);
      sum += forward_table(r, 0);
    }
    scales[0] = sum;
    for(size_t r = 0; r < no_states; ++r)
      forward_table(r, 0) /= sum;

    // recursion
    for(size_t c = 1; c < length; ++c) {
      sum = 0.0;
      for(size_t r = 0; r < no_states; ++r) {
	forward_table(r,c) = 0.0;
	for(size_t prev_state = 0; prev_state < no_states; ++prev_state) {
	  forward_table(r,c) += forward_table(prev_state, c - 1) * A(prev_state, r);
	}
	forward_table(r,c) *= B(r, seq[c]);
	sum += forward_table(r,c);
      }
      
      scales[c] = sum;
      // scale
      for(size_t r = 0; r < no_states; ++r)
	forward_table(r,c) /= sum;
    }
  }

  void backward(Matrix &backward_table, 
		const Matrix &pi, const Matrix &A, const Matrix&B, 
		const std::vector<unsigned> &seq, 
		const std::vector<double> &scales,
		const Matrix &forward_table) {
    size_t no_states = A.get_height();
    size_t length = seq.size();

    backward_table.reset(no_states, length);

    // init
    for(size_t r = 0; r < no_states; ++r) {
      backward_table(r, length - 1) = 1.0;
    }

    // recursion
    for(size_t c = length - 2; c <= length - 2; --c) { // weird because of unsigned wrap around when negative.
      for(size_t r = 0; r < no_states; ++r) {
	backward_table(r,c) = 0.0;
	for(size_t next_state = 0; next_state < no_states; ++next_state) {
	  backward_table(r, c) += backward_table(next_state, c + 1) * A(r, next_state) * B(next_state, seq[c + 1]);
	}
	
	backward_table(r,c) /= scales[c+1];
      }
    }
  }

  void posterior_decoding(std::vector<unsigned> &pd_path,
			  Matrix &pd_table,
			  const Matrix &pi,
			  const Matrix &A,
			  const Matrix &B,
			  const std::vector<unsigned> &seq) {
    
    size_t no_states = A.get_height();
    size_t length = seq.size();
    
    posterior_decoding(pd_table, pi, A, B, seq);

    pd_path.resize(length);
    for(size_t c = 0; c < length; ++c) {
      size_t max_state = 0;
      for(size_t r = 1; r < no_states; ++r) {
	if(pd_table(r, c) > pd_table(max_state, c))
	  max_state = r;
      }

      pd_path[c] = unsigned(max_state);
    }
  }

  void posterior_decoding(std::vector<unsigned> &pd_path,
			  Matrix &pd_table,
			  const Matrix &pi,
			  const Matrix &A,
			  const Matrix &B,
			  const std::string &seq_filename) {

    std::vector<unsigned> seq;
    readSeq(seq, seq_filename);
    posterior_decoding(pd_path, pd_table, pi, A, B, seq);
  }

  void posterior_decoding(Matrix &pd_table,
			  const Matrix &pi,
			  const Matrix &A,
			  const Matrix &B,
			  const std::vector<unsigned> &seq) {
    
    Matrix forward_table, backward_table;
    std::vector<double> scales;

    size_t no_states = A.get_height();
    size_t length = seq.size();

    forward(forward_table, scales, pi, A, B, seq);
    backward(backward_table, pi, A, B, seq, scales, forward_table);

    pd_table.reset(no_states, length);
    for(size_t r = 0; r < no_states; ++r) {
      for(size_t c = 0; c < length; ++c) {
	pd_table(r, c) = forward_table(r, c) * backward_table(r, c);
      }
    }
  }

  void posterior_decoding(Matrix &pd_table,
			  const Matrix &pi,
			  const Matrix &A,
			  const Matrix &B,
			  const std::string &seq_filename) {

    std::vector<unsigned> seq;
    readSeq(seq, seq_filename);
    posterior_decoding(pd_table, pi, A, B, seq);
  }

} // namespace
