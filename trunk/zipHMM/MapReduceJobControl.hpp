#ifndef MAP_REDUCE_JOB_CONTROL_HPP
#define MAP_REDUCE_JOB_CONTROL_HPP

#include "debug.hpp"

#include <vector>

namespace zipHMM {

  struct MapReduceJobControl {
    const size_t vectorLength, nBlocks;

    pthread_mutex_t mutex;
    size_t headBlock;
    std::vector<double> resultLogLikelihoods;
    
    MapReduceJobControl(size_t vectorLength, size_t nBlocks)
      : vectorLength(vectorLength), nBlocks(nBlocks),
	headBlock(0),
	resultLogLikelihoods(nBlocks, 0.0) {
      pthread_mutex_init(&mutex, 0);
    }

    ~MapReduceJobControl() {
      pthread_mutex_destroy(&mutex);
    }

    size_t getBlockBegin(size_t block) const {
      return static_cast<size_t>(static_cast<uint_least64_t>(block) * vectorLength / nBlocks);
    }
    
    size_t getBlockEnd(size_t block) const { 
      return getBlockBegin(block + 1);
    }

    size_t getMaxBlockLength() const {
      return ((vectorLength - 1) + (nBlocks - 1)) / nBlocks;
    }
  }; // struct

} // namespace

#endif
