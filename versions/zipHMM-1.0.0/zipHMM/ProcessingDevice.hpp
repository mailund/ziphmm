#ifndef PROCESSINGDEVICE_HPP
#define PROCESSINGDEVICE_HPP

#include "matrix.hpp"
#include "hmm_utils.hpp"
#include "Stage1JobControl.hpp"
#include "Stage2JobControl.hpp"
#include "ViterbiJobControl.hpp"
#include "debug.hpp"

#include <vector>

namespace zipHMM {

  class ProcessingDevice {
  
  public:
    virtual ~ProcessingDevice() { }

    virtual void computeSymbol2ScaleAndSymbol2Matrix(Stage1JobControl &control) = 0;
    
    virtual void likelihoodVector(Stage2JobControl &control) = 0;
    virtual void likelihoodMatrix(Stage2JobControl &control) = 0;

    virtual void join() = 0;

    virtual void setParameters(const Matrix *pi, const Matrix *A, const Matrix *B, 
			       const std::map<unsigned, s_pair> *symbol2pair, double *symbol2scale, Matrix *symbol2matrix, 
			       const std::vector<unsigned> *seq) = 0;

    virtual void setHMM(const Matrix *pi, const Matrix *A, const Matrix *B) = 0;
    virtual void setSeq(const std::vector<unsigned> *seq) = 0;
    virtual void setSymbol2scale(double *symbol2scale) = 0;
    virtual void setSymbol2matrix(Matrix *symbol2matrix) = 0;

    virtual const Matrix &getInitProbs() const = 0;
    virtual const Matrix &getEmProbs() const = 0;

    virtual void clear() = 0;
  };
}

#endif
