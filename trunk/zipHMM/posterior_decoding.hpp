#ifndef POSTERIOR_DECODING_HPP
#define POSTERIOR_DECODING_HPP

#include "matrix.hpp"

#include <vector>

namespace zipHMM {
  
  void posterior_decoding(std::vector<unsigned> &pd_path,
			  Matrix &pd_table,
			  const Matrix &pi,
			  const Matrix &A,
			  const Matrix &B,
			  const std::vector<unsigned> &seq);

  void posterior_decoding(std::vector<unsigned> &pd_path,
			  Matrix &pd_table,
			  const Matrix &pi,
			  const Matrix &A,
			  const Matrix &B,
			  const std::string &seq_filename);

  void posterior_decoding(Matrix &pd_table,
			    const Matrix &initProbs,
			    const Matrix &transProbs,
			    const Matrix &emProbs,
			    const std::vector<unsigned> &seq);

  void posterior_decoding(Matrix &pd_table,
			  const Matrix &pi,
			  const Matrix &A,
			  const Matrix &B,
			  const std::string &seq_filename);

}

#endif
