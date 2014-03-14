#include "calibrate.hpp"

#include "timer.hpp"
#include "forwarder.hpp"
#include "seq_io.hpp"
#include "performance_description.hpp"

#include <iostream>
#include <limits>
#include <cmath>
#include <vector>

namespace zipHMM {

  const size_t NREPS = 10;
  const size_t NSTATES = 4;
  const size_t NOBSERVABLES = 10;
  const size_t SEQLENGTH = 1000000;

  namespace {

    std::vector<double> makeProbVector(size_t size) {
      std::vector<double> result(size, 0.0);
      double sum = 0.0;
      
      for(size_t i = 0; i < size; ++i)
	sum += result[i] = rand();
      for(size_t i = 0; i < size; ++i)
	result[i] /= sum;
      
      return result;
    }

    void makeProbMatrix(Matrix &result) {
      for(size_t i = 0; i < result.get_height(); ++i) {
	std::vector<double> vals = makeProbVector(result.get_width());
	for(size_t j = 0; j < result.get_width(); ++j)
	  result(i, j) = vals[j];
      }
    }

    double testNDevices(const Forwarder &forwarder, const DeviceDescriptor &descriptor) {
      Matrix pi(makeProbVector(NSTATES));
      Matrix A(NSTATES, NSTATES);
      makeProbMatrix(A);
      Matrix B(NSTATES, NOBSERVABLES);
      makeProbMatrix(B);

      Timer timer;
      timer.start();
      forwarder.pthread_forward(pi, A, B, descriptor);
      timer.stop();

      return timer.timeElapsed();
    }
                           
    size_t findBestNDevices(const DeviceFactory &factory) {
      if(!factory.supportsMultipleDevices())
	return 1;


      std::cout << "\tCreating sequence" << std::endl;
      std::vector<unsigned> seq;
      for(size_t i = 0; i < SEQLENGTH; ++i)
	seq.push_back( (unsigned) rand() % NOBSERVABLES );
      writeSeq(seq, "tmp.seq");

      std::cout << "\tCreating forwarder" << std::endl;
      std::vector<size_t> nStatesSave;
      zipHMM::Forwarder forwarder;
      forwarder.read_seq("tmp.seq", NOBSERVABLES, nStatesSave);

      size_t bestNDevices = 1;
      double bestTime = std::numeric_limits<double>::infinity();

      for(size_t nDevices = 1; nDevices <= 2*bestNDevices + 3; ++nDevices) {
	std::cout << "    Trying " << nDevices << " devices." << std::endl;

	DeviceDescriptor descriptor(factory, nDevices);
	double timeSum = 0;
	for(size_t i = 0; i < NREPS; ++i)
	  timeSum += testNDevices(forwarder, descriptor);

	std::cout << "      " << timeSum << "s." << std::endl;
	if(timeSum < bestTime) {
	  std::cout << "      New best." << std::endl;
	  bestNDevices = nDevices;
	  bestTime = timeSum;
	}
      }

      return bestNDevices;
    }
    
    void calibrateDevice(DeviceDescriptor &result, const DeviceFactory &factory) {
      std::cout << "Calibrating device type " << factory.getTypeId() << std::endl;
      
      srand( (unsigned) time(0) );
      
      size_t nDevices = findBestNDevices(factory);
      
      result = DeviceDescriptor(factory, nDevices);
    }
  }

  void calibrateDevice(DeviceDescriptor &result, const DeviceTypeId &typeId) {
    calibrateDevice(result, deviceIdToFactory(typeId));
  }
  
  void calibrateDevices(std::vector<DeviceDescriptor> &result) {
    result.resize(N_DEVICE_FACTORIES);
    for(size_t i = 0; i < N_DEVICE_FACTORIES; ++i)
      calibrateDevice(result[i], *(ALL_DEVICE_FACTORIES[i]));
  }

  void calibrate(const std::string &filename) {
    std::string f = filename;
    if(std::strcmp(filename.c_str(), "-") == 0)
      f = zipHMM::DEFAULT_DEVICE_FILENAME;
    std::vector<zipHMM::DeviceDescriptor> descriptors;
    zipHMM::calibrateDevices(descriptors);
    writeDescriptors(descriptors, filename);
  }
}
