#ifndef VITERBI_JOB_CONTROL
#define VITERBI_JOB_CONTROL

#include "matrix.hpp"

#include <vector>

namespace zipHMM {

  struct ViterbiJobControl {
    const size_t seqLength, nBlocks;
    
    pthread_mutex_t mutex;
    size_t headBlock, tailBlock;
    std::vector<Matrix *> resultMatrices;
    std::vector<unsigned> *resultSequence;
    
    ViterbiJobControl(size_t seqLength, size_t nBlocks, std::vector<unsigned> *resultSequence = 0)
      : seqLength(seqLength), nBlocks(nBlocks),
	headBlock(0), tailBlock(nBlocks),
	resultMatrices(nBlocks, 0),
	resultSequence(resultSequence) {
      
      pthread_mutex_init(&mutex, 0);
      for(size_t i = 0; i < nBlocks; ++i)
	resultMatrices[i] = new Matrix();
    }

    ~ViterbiJobControl() {
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
  };

}

#endif
