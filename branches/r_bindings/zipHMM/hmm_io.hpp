#ifndef HMM_IO_HPP
#define HMM_IO_HPP

#include "matrix.hpp"

#include <string>
#include <iostream>

namespace zipHMM { 
  void read_HMM_spec(size_t &nStates, size_t &nObservables, const std::string &filename);
  void read_HMM_spec(size_t &nStates, size_t &nObservables, std::istream &in);

  void read_HMM(Matrix &resultInitProbs,
	       Matrix &resultTransProbs,
	       Matrix &resultEmProbs,
	       const std::string &filename);

  void read_HMM(Matrix &resultInitProbs,
	       Matrix &resultTransProbs,
	       Matrix &resultEmProbs,
	       std::istream &in);
  
  void write_HMM(const Matrix &initProbs,
		const Matrix &transProbs,
		const Matrix &emProbs,
		const std::string &filename);

  void write_HMM(const Matrix &initProbs,
		const Matrix &transProbs,
		const Matrix &emProbs,
		std::ostream &out);
}

#endif
