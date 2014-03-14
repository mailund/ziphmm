#include "PThreadProcessingDevice.hpp"

#include "prob_spaces.hpp"
#include "hmm_utils.hpp"
#include "Stage2JobControl.hpp"
#include "MapReduceJobControl.hpp"
#include "debug.hpp"

#include <pthread.h>
#include <map>
#include <deque>
#include <vector>

namespace zipHMM {
  
  namespace {

    struct InternalStage2JobControl {
      PThreadProcessingDevice &self;
      Stage2JobControl &control;
      
      InternalStage2JobControl(PThreadProcessingDevice &self, Stage2JobControl &control)
	: self(self), control(control)
      {}
    };

    struct InternalMapReduceJobControl {
      PThreadProcessingDevice &self;
      MapReduceJobControl &control;
      
      InternalMapReduceJobControl(PThreadProcessingDevice &self, MapReduceJobControl &control)
	: self(self), control(control)
      {}
    };
  }

  void processLikelihood(Matrix &resultMatrix, double &resultLogLikelihood, PThreadProcessingDevice &device, size_t seqBegin, size_t seqEnd) {
    const std::vector<unsigned> &seq = device.getSeq();
            
    Matrix *result = new Matrix();
    Matrix *temp = new Matrix();

    Matrix::copy(resultMatrix, *result);
    resultLogLikelihood = LinearSpace::toLogSpace(result->normalize());
    // std::cout << "resultLogLikelihood: " << resultLogLikelihood << std::endl;

    for(size_t i = seqBegin; i < seqEnd; ++i) {
      Matrix::mult(device.getMatrix(seq[i]), *result, *temp);
      std::swap(result, temp);
      double scale = LinearSpace::toLogSpace(result->normalize());
      // std::cout << "scale: " << scale << std::endl;
      resultLogLikelihood += scale + device.getScale(seq[i]);
    }
      
    Matrix::copy(*result, resultMatrix);
      
    delete result;
    delete temp;
  }

  void *processLikelihoodVector(void *vInternalStage2JobControl) {
    InternalStage2JobControl *internalControl = static_cast<InternalStage2JobControl *>(vInternalStage2JobControl);
    PThreadProcessingDevice &device = internalControl->self;
    Stage2JobControl &control = internalControl->control;
      
    device.clearResult();

    Matrix resultMatrix;
    double resultLogLikelihood = LogSpace::ONE;
    double tempLogLikelihood = LogSpace::ONE;

    resultLogLikelihood += std::log( init_apply_em_prob(resultMatrix, device.getInitProbs(), device.getEmProbs(), device.getSeq()[0]) );

    pthread_mutex_lock(&control.mutex);
    while(control.headBlock < control.tailBlock) {
      size_t block = control.headBlock++;
      // std::cout << "Vector took block " << block << std::endl;
      pthread_mutex_unlock(&control.mutex);
	
      processLikelihood(resultMatrix, tempLogLikelihood, device, control.getBlockBegin(block), control.getBlockEnd(block));

      resultLogLikelihood += tempLogLikelihood;
	
      Matrix::copy(resultMatrix, *(control.resultMatrices[block]));
      control.resultLogLikelihoods[block] = resultLogLikelihood;
	
      pthread_mutex_lock(&control.mutex);
    }
    pthread_mutex_unlock(&control.mutex);

    delete internalControl;
    return 0;
  }

  void *processLikelihoodMatrix(void *vInternalStage2JobControl) {
    InternalStage2JobControl *internalControl = static_cast<InternalStage2JobControl *>(vInternalStage2JobControl);
    PThreadProcessingDevice &device = internalControl->self;
    Stage2JobControl &control = internalControl->control;

    device.clearResult();

    const size_t N = device.getInitProbs().get_height();

    pthread_mutex_lock(&control.mutex);
    while((control.tailBlock - control.headBlock) > N && control.tailBlock > 1) {
      size_t block = --control.tailBlock;
      // std::cout << "Matrix took block " << block << std::endl;
      pthread_mutex_unlock(&control.mutex);
	
      Matrix resultMatrix;
      double resultLogLikelihood = LogSpace::ONE;
	
      Matrix::copy( device.getMatrix(device.getSeq()[control.getBlockBegin(block)]) , resultMatrix );
      resultLogLikelihood += device.getScale( device.getSeq()[control.getBlockBegin(block)] );
      // std::cout << "resultLogLikelihood: " << resultLogLikelihood << std::endl;

      double tmpLoglikelihood = LogSpace::ONE;
      processLikelihood(resultMatrix, tmpLoglikelihood, device, control.getBlockBegin(block) + 1, control.getBlockEnd(block));
      resultLogLikelihood += tmpLoglikelihood;

      Matrix::copy(resultMatrix, *(control.resultMatrices[block]));
      control.resultLogLikelihoods[block] = resultLogLikelihood;
	
      pthread_mutex_lock(&control.mutex);
    }
    pthread_mutex_unlock(&control.mutex);

    delete internalControl;
    return 0;
  }

  void PThreadProcessingDevice::likelihoodVector(Stage2JobControl &control) {
    pthread_create(&thread, 0, processLikelihoodVector, new InternalStage2JobControl(*this, control));
  }

  void PThreadProcessingDevice::likelihoodMatrix(Stage2JobControl &control) {
    pthread_create(&thread, 0, processLikelihoodMatrix, new InternalStage2JobControl(*this, control));
  }

  void *processMapReduceLoglikelihood(void *vInternalMapReduceJobControl) {
    InternalMapReduceJobControl *internalControl = static_cast<InternalMapReduceJobControl *>(vInternalMapReduceJobControl);
    PThreadProcessingDevice &device = internalControl->self;
    MapReduceJobControl &control = internalControl->control;

    std::deque<double> scales;
    Matrix resultMatrix;
    Matrix tmp;
    
    device.clearResult();

    pthread_mutex_lock(&control.mutex);
    while(control.headBlock < control.nBlocks) {
      size_t block = control.headBlock++;
      // std::cout << "mr took block " << block << std::endl;
      pthread_mutex_unlock(&control.mutex);

      for(size_t i = control.getBlockBegin(block); i < control.getBlockEnd(block); ++i) {
	// std::cout << "block " << block << " processing sequence " << i << std::endl;
	
	const std::vector<unsigned> &sequence = device.getSeq(i);
	
	scales.clear();

	// compute C_1 and push corresponding scale
	scales.push_back(std::log(init_apply_em_prob(resultMatrix, device.getInitProbs(), device.getEmProbs(), sequence[0])));
	
	// multiply matrices across the sequence
	for(size_t i = 1; i < sequence.size(); ++i) {
	  Matrix::blas_mult(device.getMatrix(sequence[i]), resultMatrix, tmp);
	  Matrix::copy(tmp, resultMatrix);
	  scales.push_back( std::log( resultMatrix.normalize() ) );
	  scales.push_back( device.getScale(sequence[i]) );
	}
	
	// compute loglikelihood by summing log of scales
	double loglikelihood = 0.0;
	for(std::deque<double>::iterator it = scales.begin(); it != scales.end(); ++it) {
	  loglikelihood += (*it);
	}
	

	control.resultLogLikelihoods[block] += loglikelihood;
      }

      pthread_mutex_lock(&control.mutex);
    }
    pthread_mutex_unlock(&control.mutex);
    
    delete internalControl;
    return 0;
  }

  void PThreadProcessingDevice::mapReduceLoglikelihood(MapReduceJobControl &control) {
    pthread_create(&thread, 0, processMapReduceLoglikelihood, new InternalMapReduceJobControl(*this, control));
  }

 
  void PThreadProcessingDevice::setParameters(const Matrix *pi, const Matrix *A, const Matrix *B, 
					      const std::map<unsigned, s_pair> *symbol2pair, double *symbol2scale, Matrix *symbol2matrix) {
    setHMM(pi, A, B);
    setSymbol2pair(symbol2pair);
    setSymbol2scale(symbol2scale);
    setSymbol2matrix(symbol2matrix);
  }

  void PThreadProcessingDevice::setHMM(const Matrix *pi, const Matrix *A, const Matrix *B) {
    this->pi = pi;
    this->A = A;
    this->B = B;
  }

  void PThreadProcessingDevice::join() {
    pthread_join(thread, 0);
  }

  void PThreadProcessingDevice::computeSymbol2logMatrix() {
    if(symbol2logMatrix == 0) {
      Matrix *symbol2logMatrix = new Matrix[B->get_width()];
      const size_t nStates = A->get_height();

      for(size_t symbol = 0; symbol < B->get_width(); ++symbol) {
	symbol2logMatrix[symbol].reset(nStates, nStates);
	for(size_t r = 0; r < nStates; ++r) {
	  for(size_t c = 0; c < nStates; ++c)
	    symbol2logMatrix[symbol](r, c) = A->at(c, r) * B->at(r, symbol);
	}
	symbol2logMatrix[symbol].log();
      }
    }
  }

  void PThreadProcessingDevice::clearCache() { // delete?
    if(symbol2logMatrix != 0)
      delete[] symbol2logMatrix;
  }

  void PThreadProcessingDevice::clearResult() { // delete?
    
  }

  void PThreadProcessingDevice::clear() { // delete?
    clearCache();
    clearResult();
  }

} // namespace zipHMM
