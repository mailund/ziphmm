#ifndef MATRIX_HPP
#define MATRIX_HPP

#include "prob_spaces.hpp"
#include "timer.hpp"

#include <iostream>
#include <cmath>
#include <cstring>
#include <vector>
#include <algorithm>
#include <cstdlib>

#if defined __APPLE__
  #include <Accelerate/Accelerate.h>
#else

extern "C" {
  #include "cblas.h"
}

#endif

namespace zipHMM {
  class Matrix {

    double *content;
    size_t height, width;

    size_t idx(size_t i, size_t j) const
    {
      return width * i + j;
    }

    //Not implemented. Do not copy matrices.
    Matrix(const Matrix&);
    Matrix &operator=(const Matrix&);

  public:

    Matrix()
      : content(new double[1]), // allocating array of size 1 (not 0) to avoid warning.
	height(0), 
	width(0)
    {}

    Matrix(size_t height, size_t width)
      : content(new double[width * height]), height(height), width(width) { }

    Matrix(const std::vector<double> &column)
      : content(new double[column.size()]), height(column.size()), width(1) {
      for(size_t i = 0; i < column.size(); ++i)
	content[i] = column[i];
    }
    
    Matrix(const std::vector<std::vector<double> > &matrix)
      : content(0), height(0), width(0) {
      if(matrix.size() == 0) {
	height = 0;
	width = 0;
	content = new double[1]; // allocating array of length 1 (not 0) to avoid warnings.
	return;
      }
      
      height = matrix.size();
      width = matrix[0].size();
      content = new double[width * height];
      
      for(size_t i = 0; i < height; ++i) {
	for(size_t j = 0; j < width; ++j)
	  content[idx(i, j)] = matrix[i][j];
      }
    }

    ~Matrix() {
      delete[] content;
    }

    size_t get_width()  const { return width; }
    size_t get_height() const { return height; }

    void reset(size_t height, size_t width) {
      if(height != this->height || width != this->width) {
	delete[] content;
	this->height = height;
	this->width = width;
	content = new double[width * height];
      }
    }
    
    //Accessors
    double &at(size_t i, size_t j)       { return content[idx(i, j)]; }
    double  at(size_t i, size_t j) const { return content[idx(i, j)]; }

    //For convenience in C++
    double &operator()(size_t i, size_t j)       { return at(i, j); }
    double  operator()(size_t i, size_t j) const { return at(i, j); }

    //For raw access in C++
    const double *get_raw_data() const { return content; }

    //For exporting to scripting languages
    void   set(size_t i, size_t j, double value) { at(i, j) = value; }
    double get(size_t i, size_t j) const         { return at(i, j); }

    void rowSwap(size_t r1, size_t r2) { 
      for(size_t i = 0; i < width; ++i) 
	std::swap(at(r1, i), at(r2, i)); 
    }

    double normalize() {
      double sum = 0.0;
      
      for(size_t i = 0; i < width*height; ++i)
	sum += content[i];
      for(size_t i = 0; i < width*height; ++i)
	content[i] /= sum;
      
      return sum;
    }

    void log() {
      for(unsigned i = 0; i < width*height; ++i)
	content[i] = std::log(content[i]);
    }

    unsigned colArgMax(unsigned j) {
      double bestVal = LogSpace::ZERO;
      unsigned bestIdx = static_cast<unsigned>(-1);
      
      for(unsigned i = 0; i < height; ++i) {
	if(at(i, j) >= bestVal) {
	  bestVal = at(i, j);
	  bestIdx = i;
	}
      }
      
      return bestIdx;
    }

    void print() const {
      for(size_t r = 0; r < height; ++r) {
	for(size_t c = 0; c < width; ++c) 
	  std::cout << content[idx(r,c)] << "\t";
	std::cout << std::endl;
      }
    }

    static inline
    void mult(const Matrix &lhs, const Matrix &rhs, Matrix &res) {
      const size_t depth = lhs.width;
      const size_t width = rhs.width;
      
      res.reset(lhs.height, width);
      
      double *outPtr = res.content;
      const double *lhsPtr = lhs.content;
      const double *rhsPtr = rhs.content;
      
      const double * const lhsStop = lhsPtr + lhs.width * lhs.height;
      const double * const rhsStop = rhsPtr + width;
      
      for(; lhsPtr < lhsStop; lhsPtr += depth, rhsPtr -= width) {
	for(; rhsPtr < rhsStop; lhsPtr -= depth, rhsPtr -= (depth-1)*width - 1) {
	  double val = (*lhsPtr++) * (*rhsPtr);
	  for(size_t k = 1; k < depth; ++k) {
	    rhsPtr += width;
	    val = val + (*lhsPtr++) * (*rhsPtr);
	  }
	  
	  *outPtr++ = val;
	}
      }
    }

    static inline
    void blas_mult(const Matrix &lhs, const Matrix &rhs, Matrix &res) {
      const size_t height = lhs.height;
      const size_t width = rhs.width;
      const size_t depth = lhs.width;
      
      res.reset(height, width);
      
      cblas_dgemm(CblasRowMajor,
		  CblasNoTrans,
		  CblasNoTrans,
		  static_cast<int>(height), 
		  static_cast<int>(width), 
		  static_cast<int>(depth),
		  1.0,
		  lhs.content, 
		  static_cast<int>(lhs.width), // changed
		  rhs.content, 
		  static_cast<int>(rhs.width), // changed
		  0.0,
		  res.content, 
		  static_cast<int>(res.width)); // changed
    }

    static inline
    void blas_matrix_vector_mult(const Matrix &lhs, const Matrix &rhs, Matrix &res) {
      const size_t nRows = lhs.height;
      const size_t nCols = lhs.width;
      
      res.reset(nRows, 1);
      
      cblas_dgemv(CblasRowMajor,
		  CblasNoTrans,
		  static_cast<int>(nRows), 
		  static_cast<int>(nCols),
		  1.0,
		  lhs.content,
		  static_cast<int>(nRows),
		  rhs.content, 
		  1,
		  0.0,
		  res.content, 
		  1);
    }

    template<typename Space>
    static inline
    void maxMult(const Matrix &lhs, const Matrix &rhs, Matrix &res) {
      const size_t depth = lhs.width;
      
      res.reset(lhs.height, rhs.width);
      
      double *outPtr = res.content;
      const double *lhsPtr = lhs.content;
      
      for(size_t i = 0; i < lhs.height; ++i) {
	for(size_t j = 0; j < rhs.width; ++j) {
	  double val = Space::mult(*lhsPtr++, rhs(0, j));
	  for(size_t k = 1; k < depth; ++k) {
	    double newVal = Space::mult(*lhsPtr++, rhs(k, j));
	    if(newVal > val)
	      val = newVal;
	  }
	  lhsPtr -= depth;
	  *outPtr++ = val;
	}
	lhsPtr += depth;
      }
    }
    
    template<typename Space>
    static inline
    size_t argMaxMult(const Matrix &lhs, size_t leftRow, const Matrix &rhs, size_t rightCol) {
      const size_t depth = lhs.width;
      
      double bestVal = Space::mult(lhs(leftRow, 0), rhs(0, rightCol));
      size_t bestIdx = 0;
      
      for(size_t i = 1; i < depth; ++i) {
	double val = Space::mult(lhs(leftRow, i), rhs(i, rightCol));
	
	if(val >= bestVal) {
	  bestVal = val;
	  bestIdx = i;
	}
      }
      
      return bestIdx;
    }
    
    
    static inline
    void copy(const Matrix &from, Matrix &to) {
      to.reset(from.get_height(), from.get_width());
      memcpy(to.content, from.content, sizeof(double) * from.get_height() * from.get_width());
    }

    static inline
    void transpose(const Matrix &from, Matrix &to) {
      to.reset(from.get_height(), from.get_width());
      for(size_t r = 0; r < from.get_height(); ++r) 
	for(size_t c = 0; c < from.get_width(); ++c)
	  to(c, r) = from(r, c);
    }

    static inline
    void copy(const double *from, size_t height, size_t width, Matrix &to) {
      to.reset(height, width);
      memcpy(to.content, from, sizeof(double) * height * width);
    }

    #define TIME_BLAS_MULT_TIMES 200

    static inline
    double time_blas_mult(const size_t N) {
      Matrix M(N,N);
      Matrix tmp;
      Matrix res;

      for(size_t r = 0; r < N; ++r)
	for(size_t c = 0; c < N; ++c)
	  M(r, c) = ((double) std::rand() / (RAND_MAX));
      copy(M, tmp);

      Timer timer;
      timer.start();

      for(unsigned i = 0; i < TIME_BLAS_MULT_TIMES; ++i) {
	blas_mult(tmp, M, res);
	res.normalize();
	copy(res, tmp);
      }

      timer.stop();
      
      return timer.timeElapsed() / TIME_BLAS_MULT_TIMES;      
    }

    #define TIME_BLAS_MATRIX_VECTOR_MULT_TIMES 200
    
    static inline 
    double time_blas_matrix_vector_mult(const size_t N) {
      Matrix M(N, N);
      Matrix V(N, 0);
      Matrix res;

      for(size_t r = 0; r < N; ++r) {
	for(size_t c = 0; c < N; ++c)
	  M(r, c) = ((double) std::rand() / (RAND_MAX));
	V(r, 0) = ((double) std::rand() / (RAND_MAX));
      }

      Timer timer;
      timer.start();

      for(unsigned i = 0; i < TIME_BLAS_MATRIX_VECTOR_MULT_TIMES; ++i) {
	blas_matrix_vector_mult(M, V, res);
	res.normalize();
	copy(res, M);
      }

      timer.stop();
      
      return timer.timeElapsed() / TIME_BLAS_MATRIX_VECTOR_MULT_TIMES;      
    }

  };
} // namespace

std::ostream &operator<<(std::ostream &out, const zipHMM::Matrix &mat);

#endif
