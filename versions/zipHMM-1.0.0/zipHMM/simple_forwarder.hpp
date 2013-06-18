#ifndef SIMPLE_FORWARD_HPP
#define SIMPLE_FORWARD_HPP

#include "matrix.hpp"
#include "hmm_utils.hpp"
#include "seq_io.hpp"
#include "hmm_io.hpp"

#include <vector>
#include <deque>
#include <string>

namespace zipHMM {
  
  class SimpleForwarder {
  private:
    std::vector<unsigned> seq;
    
  public:
    SimpleForwarder(const std::string &seq_filename) {
      seq = std::vector<unsigned>();
      readSeq(seq, seq_filename);
    }

    double forward(const Matrix &pi, const Matrix &A, const Matrix &B) const {
      std::deque<double> scales;
      double *symbol2scale = new double[B.get_width()];
      Matrix *symbol2matrix = new Matrix[B.get_width()];
      Matrix res;
      Matrix tmp;
      double loglikelihood;

      make_em_trans_probs_array(symbol2scale, symbol2matrix, A, B);

      // compute C_1 and push corresponding scale
      scales.push_back( std::log( init_apply_em_prob(res, pi, B, seq[0]) ) );
    
      // multiply matrices across the sequence
      for(size_t i = 1; i < seq.size(); ++i) {
	Matrix::blas_mult(symbol2matrix[seq[i]], res, tmp);
	Matrix::copy(tmp, res);
	scales.push_back( std::log( res.normalize() ) );
	scales.push_back( symbol2scale[seq[i]] );
      }
    
      // compute loglikelihood by summing log of scales
      loglikelihood = 0.0;
      for(std::deque<double>::iterator it = scales.begin(); it != scales.end(); ++it) {
	loglikelihood += (*it);
      }
    
      delete[] symbol2scale;
      delete[] symbol2matrix;
    
      return loglikelihood;
    }
  };

}

#endif
