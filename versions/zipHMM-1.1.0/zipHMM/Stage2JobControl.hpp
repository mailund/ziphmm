#ifndef STAGE_2_JOB_CONTROL_HPP
#define STAGE_2_JOB_CONTROL_HPP

#include "debug.hpp"

#include <vector>

namespace zipHMM {

  struct Stage2JobControl {
    const size_t seqLength, nBlocks;

    pthread_mutex_t mutex;
    size_t headBlock, tailBlock;
    std::vector<Matrix *> resultMatrices;
    std::vector<double> resultLogLikelihoods;
    
    Stage2JobControl(size_t seqLength, size_t nBlocks)
      : seqLength(seqLength), nBlocks(nBlocks),
	headBlock(0), tailBlock(nBlocks),
	resultMatrices(nBlocks, 0),
	resultLogLikelihoods(nBlocks, 0.0) {
      pthread_mutex_init(&mutex, 0);
      for(size_t i = 0; i < nBlocks; ++i)
	resultMatrices[i] = new Matrix();
    }

    ~Stage2JobControl() {
      pthread_mutex_destroy(&mutex);
      for(size_t i = 0; i < nBlocks; ++i)
	delete resultMatrices[i];
    }

    size_t getBlockBegin(size_t block) const {
      return static_cast<size_t>(static_cast<uint_least64_t>(block) * (seqLength-1) / nBlocks + 1);
    }
    
    size_t getBlockEnd(size_t block) const { return getBlockBegin(block + 1); }

    size_t getMaxBlockLength() const {
      return ((seqLength - 1) + (nBlocks - 1)) / nBlocks;
    }
  }; // struct

} // namespace

#endif
