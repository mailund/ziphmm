#ifndef POSTERIOR_DECODING_HPP
#define POSTERIOR_DECODING_HPP

#include "matrix.hpp"

#include <vector>

namespace zipHMM {
  
  void posterior_decoding(const std::vector<unsigned> &seq,
			  const Matrix &pi,
			  const Matrix &A,
			  const Matrix &B,
			  std::vector<unsigned> &pd_path,
			  Matrix &pd_table);

  void posterior_decoding(const std::string &seq_filename,
			  const Matrix &pi,
			  const Matrix &A,
			  const Matrix &B,
			  std::vector<unsigned> &pd_path,
			  Matrix &pd_table);

  void posterior_decoding(const std::vector<unsigned> &seq,
			  const Matrix &initProbs,
			  const Matrix &transProbs,
			  const Matrix &emProbs,
			  Matrix &pd_table);

  void posterior_decoding(const std::string &seq_filename,
			  const Matrix &pi,
			  const Matrix &A,
			  const Matrix &B,
			  Matrix &pd_table);

}

#endif
